"""Seed the Filo3D structure set from the RCSB PDB.

Two-step, fully reproducible acquisition:

1. **Search API** (search.rcsb.org): find every experimental polymer entity whose
   source-organism lineage contains Filoviridae (NCBI taxid 11266). Using the
   lineage id captures all genera/species below the family in one query.
2. **Data API / GraphQL** (data.rcsb.org): batch-fetch metadata for each hit
   (description, organism, resolution, method, release date, entry-level partner
   descriptions and ligands).

Output: a per-polymer-entity ``raw.csv`` (one row per protein chain-group). The
downstream ``annotate`` stage classifies each row into a canonical filovirus
protein + species; ``download`` fetches the CIF files and quality metrics.
"""

import sys
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

import pandas as pd
import requests

# Make ``src`` importable when run as a module from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src._config import load_params  # noqa: E402

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

CONTENT_TYPE_MAP = {"experimental": "experimental", "computational": "computational"}

GRAPHQL_QUERY = """
query($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    rcsb_polymer_entity { pdbx_description }
    entity_poly {
      rcsb_sample_sequence_length
      rcsb_entity_polymer_type
      pdbx_seq_one_letter_code_can
    }
    rcsb_entity_source_organism { ncbi_scientific_name ncbi_taxonomy_id }
    rcsb_polymer_entity_container_identifiers { entry_id auth_asym_ids }
    entry {
      rcsb_id
      struct { title }
      exptl { method }
      rcsb_accession_info { initial_release_date }
      rcsb_entry_info { resolution_combined }
      polymer_entities {
        rcsb_polymer_entity { pdbx_description }
        rcsb_entity_source_organism { ncbi_scientific_name }
      }
      nonpolymer_entities { nonpolymer_comp { chem_comp { id name } } }
    }
  }
}
"""


def build_search_query(taxonomy_id, return_type, content_types):
    """Search for polymer entities whose source-organism lineage includes taxid."""
    return {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                "operator": "exact_match",
                "value": str(taxonomy_id),
            },
        },
        "return_type": return_type,
        "request_options": {
            "return_all_hits": True,
            "results_content_type": [
                CONTENT_TYPE_MAP.get(c, c) for c in content_types
            ],
        },
    }


def run_search(query, retries=3, timeout=60):
    for attempt in range(retries):
        try:
            resp = requests.post(SEARCH_URL, json=query, timeout=timeout)
            if resp.status_code == 204:  # no content = zero hits
                return []
            resp.raise_for_status()
            data = resp.json()
            return [hit["identifier"] for hit in data.get("result_set", [])]
        except requests.RequestException as exc:
            if attempt == retries - 1:
                raise
            print(f"  search retry {attempt + 1}/{retries}: {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    return []


def fetch_metadata(entity_ids, batch_size=200, retries=3, timeout=120):
    """Batch-fetch polymer-entity metadata via the GraphQL Data API."""
    from tqdm import tqdm

    entities = []
    for i in tqdm(range(0, len(entity_ids), batch_size), desc="GraphQL batches"):
        batch = entity_ids[i : i + batch_size]
        for attempt in range(retries):
            try:
                resp = requests.post(
                    GRAPHQL_URL,
                    json={"query": GRAPHQL_QUERY, "variables": {"ids": batch}},
                    timeout=timeout,
                )
                resp.raise_for_status()
                payload = resp.json()
                if "errors" in payload:
                    raise RuntimeError(payload["errors"])
                entities.extend(payload["data"]["polymer_entities"])
                break
            except (requests.RequestException, RuntimeError) as exc:
                if attempt == retries - 1:
                    raise
                print(f"  batch retry {attempt + 1}/{retries}: {exc}", file=sys.stderr)
                time.sleep(2 * (attempt + 1))
    return entities


def _first(lst):
    return lst[0] if lst else None


def flatten(entity):
    """Flatten one GraphQL polymer-entity record into a flat metadata row."""
    if entity is None:
        return None
    ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
    poly = entity.get("entity_poly") or {}
    desc = (entity.get("rcsb_polymer_entity") or {}).get("pdbx_description")

    organisms = entity.get("rcsb_entity_source_organism") or []
    org_names = [o.get("ncbi_scientific_name") for o in organisms if o.get("ncbi_scientific_name")]
    org_taxids = [o.get("ncbi_taxonomy_id") for o in organisms if o.get("ncbi_taxonomy_id") is not None]

    entry = entity.get("entry") or {}
    entry_info = entry.get("rcsb_entry_info") or {}
    resolution = _first(entry_info.get("resolution_combined") or [])
    methods = [m.get("method") for m in (entry.get("exptl") or []) if m.get("method")]
    release = (entry.get("rcsb_accession_info") or {}).get("initial_release_date")
    title = (entry.get("struct") or {}).get("title")

    # Entry-level sibling descriptions (for complex / partner detection downstream).
    sibling_descs = []
    for pe in entry.get("polymer_entities") or []:
        d = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description")
        if d:
            sibling_descs.append(d)

    ligands = []
    for ne in entry.get("nonpolymer_entities") or []:
        comp = ((ne.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        if comp.get("id"):
            ligands.append(comp["id"])

    return {
        "entity_id": entity.get("rcsb_id"),
        "PDB": (ids.get("entry_id") or "").lower(),
        "description": desc,
        "organism": " | ".join(org_names),
        "taxid": ";".join(str(t) for t in org_taxids),
        "seq_length": poly.get("rcsb_sample_sequence_length"),
        "polymer_type": poly.get("rcsb_entity_polymer_type"),
        "sequence": (poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", ""),
        "auth_chains": ",".join(ids.get("auth_asym_ids") or []),
        "title": title,
        "method": ";".join(methods),
        "resolution": resolution,
        "release_date": (release or "")[:10],
        "entry_descriptions": " | ".join(sibling_descs),
        "ligands": ",".join(sorted(set(ligands))),
    }


def parse_args():
    parser = ArgumentParser(
        description="Search RCSB for filovirus polymer entities and fetch metadata.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("params", help="Path to params.yaml")
    parser.add_argument("output_csv", help="Output path for the raw per-entity metadata CSV")
    return parser.parse_args()


def main():
    args = parse_args()
    params = load_params(args.params)
    s = params["search"]

    query = build_search_query(s["taxonomy_id"], s["return_type"], s["content_types"])
    print(f"# Searching RCSB for taxonomy lineage {s['taxonomy_id']} (Filoviridae)...")
    ids = run_search(query)
    print(f"# Search hits ({s['return_type']}): {len(ids)}")
    if not ids:
        print("No hits returned; writing empty CSV.", file=sys.stderr)

    entities = fetch_metadata(ids, batch_size=s.get("batch_size", 200)) if ids else []
    rows = [r for r in (flatten(e) for e in entities) if r is not None]

    df = pd.DataFrame(rows).sort_values(["PDB", "entity_id"]).reset_index(drop=True)
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    print(f"# Wrote {len(df)} entity rows across {df['PDB'].nunique()} structures -> {args.output_csv}")


if __name__ == "__main__":
    main()
