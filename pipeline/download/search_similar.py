"""Sequence-similarity completeness backstop for the Filo3D structure set.

The taxonomy search (``search_rcsb``) only finds structures RCSB has annotated
with a Filoviridae source organism. This stage runs an RCSB **sequence** search
(mmseqs2) against curated, tag-free filovirus reference sequences (one per
protein, in ``references.fa``) to catch structures the taxonomy query misses:
mis-annotated organisms, synthetic constructs, antibody-antigen depositions and
designed immunogens.

Every hit is annotated with its best alignment (which reference, %identity,
e-value, reference coverage). Hits **not** already in the taxonomy set are the
review queue: they are flagged ``in_taxonomy=False`` and classified with the
same protein/species rules, but are NOT auto-promoted into the trusted set --
you vet ``candidates.csv`` and add the good ones. This makes the manual
"search similar sequences in the PDB" step reproducible and re-runnable.
"""

import sys
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src._config import load_json, load_params  # noqa: E402
from pipeline.download.annotate import classify_protein, classify_species  # noqa: E402
from pipeline.download.search_rcsb import GRAPHQL_URL, fetch_metadata, flatten  # noqa: E402

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def read_fasta(path):
    """Minimal FASTA reader -> list of (protein_label, sequence)."""
    records, name, seq = [], None, []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name = line[1:].split()[0]  # first token = protein label (e.g. GP)
                seq = []
            elif line:
                seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def sequence_search(seq, identity, evalue, retries=3, timeout=60):
    """Verbose RCSB sequence search -> list of (entity_id, match_context)."""
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": evalue,
                "identity_cutoff": identity,
                "sequence_type": "protein",
                "value": seq,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "return_all_hits": True,
            "results_content_type": ["experimental"],
            "results_verbosity": "verbose",
        },
    }
    for attempt in range(retries):
        try:
            resp = requests.post(SEARCH_URL, json=query, timeout=timeout)
            if resp.status_code == 204:
                return []
            resp.raise_for_status()
            out = []
            for hit in resp.json().get("result_set", []):
                ctx = {}
                for svc in hit.get("services", []):
                    for node in svc.get("nodes", []):
                        mc = node.get("match_context") or []
                        if mc:
                            ctx = mc[0]
                            break
                out.append((hit["identifier"], ctx))
            return out
        except requests.RequestException as exc:
            if attempt == retries - 1:
                raise
            print(f"  seq-search retry {attempt + 1}/{retries}: {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    return []


def parse_args():
    parser = ArgumentParser(
        description="Sequence-similarity completeness backstop (review queue).",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("params")
    parser.add_argument("taxonomy_csv", help="raw.csv from search_rcsb (the trusted set)")
    parser.add_argument("output_csv", help="Output candidates.csv (review queue)")
    return parser.parse_args()


def main():
    args = parse_args()
    params = load_params(args.params)
    s = params["similar"]
    identity, evalue = s["identity_cutoff"], s["evalue_cutoff"]
    min_cov = s.get("min_ref_coverage", 0.0)

    tax_df = pd.read_csv(args.taxonomy_csv)
    tax_ids = set(tax_df["entity_id"])

    refs = read_fasta(s["references"])
    print(f"# References: {len(refs)}  ({', '.join(p for p, _ in refs)})")
    print(f"# Sequence search @ identity>={identity}, e<={evalue}, ref_coverage>={min_cov}")

    # Best match per entity across all references.
    best = {}  # entity_id -> dict(ref, identity, evalue, ref_cov, hit_cov)
    for label, seq in refs:
        hits = sequence_search(seq, identity, evalue)
        qlen = len(seq)
        kept = 0
        for eid, ctx in hits:
            aln = ctx.get("alignment_length") or 0
            ref_cov = aln / qlen if qlen else 0.0
            if ref_cov < min_cov:
                continue
            kept += 1
            ident = ctx.get("sequence_identity")
            rec = {
                "matched_ref": label,
                "pct_identity": round(ident * 100, 1) if ident is not None else None,
                "evalue": ctx.get("evalue"),
                "ref_coverage": round(ref_cov, 2),
                "hit_coverage": round(aln / ctx["subject_length"], 2)
                if ctx.get("subject_length") else None,
            }
            if eid not in best or (rec["pct_identity"] or 0) > (best[eid]["pct_identity"] or 0):
                best[eid] = rec
        print(f"  {label:5s} hits={len(hits):4d} kept={kept:4d}")
        time.sleep(0.2)

    all_ids = sorted(best)
    new_ids = [e for e in all_ids if e not in tax_ids]
    print(f"\n# Union hits: {len(all_ids)}  |  in taxonomy set: {len(all_ids) - len(new_ids)}"
          f"  |  NEW candidates: {len(new_ids)}")

    # Fetch metadata + classify the NEW candidates.
    protein_rules = load_json(params["reference"]["proteins"])["proteins"]
    taxonomy = load_json(params["reference"]["taxonomy"])
    rows = []
    if new_ids:
        entities = fetch_metadata(new_ids, batch_size=params["search"].get("batch_size", 200))
        meta = {e["rcsb_id"]: flatten(e) for e in entities if e}
        for eid in new_ids:
            m = meta.get(eid, {})
            b = best[eid]
            rows.append({
                "entity_id": eid,
                "PDB": m.get("PDB", eid.split("_")[0].lower()),
                "in_taxonomy": False,
                **b,
                "Protein_guess": classify_protein(
                    m.get("description"), m.get("polymer_type"), protein_rules),
                "Code_guess": classify_species(
                    m.get("taxid", ""), m.get("organism", ""), taxonomy)["code"],
                "description": m.get("description"),
                "organism": m.get("organism"),
                "title": m.get("title"),
                "Review": "",  # fill: keep / drop
            })

    cand = pd.DataFrame(rows)
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    cand.to_csv(args.output_csv, index=False)
    print(f"# Wrote {len(cand)} new candidates -> {args.output_csv}")
    if not cand.empty:
        print("\n=== NEW candidates by reference / guessed protein ===")
        print(cand.groupby(["matched_ref", "Protein_guess"]).size().to_string())
        print("\n=== preview ===")
        cols = ["entity_id", "matched_ref", "pct_identity", "ref_coverage",
                "Code_guess", "organism", "description"]
        print(cand[cols].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
