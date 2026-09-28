"""Classify raw RCSB entities into canonical filovirus proteins and species.

Reads the per-entity ``raw.csv`` produced by ``search_rcsb`` and adds:

- ``Protein``   canonical filovirus protein (NP, VP35, VP40, GP, VP30, VP24, L),
                or the polymer type (e.g. RNA) for non-protein entities, else Other.
- ``Code``/``Virus``/``Species``/``Genus``  filovirus taxonomy assignment.
- ``Chimeric``  the entity is a fusion of >1 distinct protein or carries a
                non-viral expression/solubility tag.
- ``HasAntibody``/``Receptor``/``Complex``  entry-level context flags.

Writes an augmented ``aug.csv`` (one row per polymer entity) and prints a
protein x species summary.
"""

import re
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from functools import lru_cache
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src._config import load_json, load_params  # noqa: E402

# Non-viral fusion / expression tags that mark a construct as chimeric.
FUSION_TAGS = [
    "maltose", "maltodextrin", "mbp", "glutathione", "gst", "sumo", "thioredoxin",
    "trx", "green fluorescent", "gfp", "superfolder", "t4 lysozyme", "lysozyme",
    "fc region", "immunoglobulin g", "bril", "apocytochrome", "designed",
]


@lru_cache(maxsize=4096)
def _kw_re(kw):
    """Whole-token matcher: alphanumeric neighbours block a match, so 'np'
    matches 'NP' / 'EBOV-NP' but not 'input', and 'gp' matches 'Shed GP'."""
    return re.compile(r"(?<![a-z0-9])" + re.escape(kw) + r"(?![a-z0-9])")


def _has_kw(text, keywords):
    return any(_kw_re(kw).search(text) for kw in keywords)


def classify_protein(desc, polymer_type, protein_rules):
    """First-match protein classification; non-protein entities keep their type."""
    if isinstance(polymer_type, str) and polymer_type.strip() and polymer_type != "Protein":
        return polymer_type  # e.g. RNA, DNA
    if not isinstance(desc, str) or not desc.strip():
        return "Other"
    low = desc.lower()
    for rule in protein_rules:
        if _has_kw(low, rule["keywords"]):
            return rule["name"]
    return "Other"


def matched_proteins(text, protein_rules):
    """Set of canonical proteins whose keywords appear anywhere in text."""
    low = text.lower()
    return {r["name"] for r in protein_rules if _has_kw(low, r["keywords"])}


def classify_species(taxid_field, organism_field, taxonomy):
    by_taxid = taxonomy["by_taxid"]
    for tid in str(taxid_field).split(";"):
        tid = tid.strip()
        if tid in by_taxid:
            return by_taxid[tid]
    low = str(organism_field).lower()
    for rule in taxonomy["by_name"]:
        if rule["match"] in low:
            return rule
    return {"code": "Unknown", "virus": "Unknown", "species": "Unknown", "genus": "Unknown"}


def is_chimeric(desc, protein_rules):
    if not isinstance(desc, str) or "," not in desc:
        return False
    low = desc.lower()
    if any(tag in low for tag in FUSION_TAGS):
        return True
    # Comma-joined parts that resolve to >1 distinct filovirus protein = true fusion.
    parts = [p.strip() for p in desc.split(",")]
    proteins = set()
    for part in parts:
        proteins |= matched_proteins(part, protein_rules)
    return len(proteins) > 1


def parse_args():
    parser = ArgumentParser(
        description="Annotate raw filovirus entities with protein + species.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("params")
    parser.add_argument("input_csv", help="raw.csv from search_rcsb")
    parser.add_argument("output_csv", help="Output aug.csv")
    return parser.parse_args()


def main():
    args = parse_args()
    params = load_params(args.params)
    ref = params["reference"]
    protein_rules = load_json(ref["proteins"])["proteins"]
    partner_rules = load_json(ref["proteins"])["partners"]
    taxonomy = load_json(ref["taxonomy"])

    df = pd.read_csv(args.input_csv)

    df["Protein"] = [
        classify_protein(d, pt, protein_rules)
        for d, pt in zip(df["description"], df["polymer_type"])
    ]

    tax = [classify_species(t, o, taxonomy) for t, o in zip(df["taxid"], df["organism"])]
    df["Code"] = [t["code"] for t in tax]
    df["Virus"] = [t["virus"] for t in tax]
    df["Species"] = [t["species"] for t in tax]
    df["Genus"] = [t["genus"] for t in tax]

    df["Chimeric"] = [is_chimeric(d, protein_rules) for d in df["description"]]

    # Entry-level partner / complex context (from sibling entity descriptions).
    entry_desc = df["entry_descriptions"].fillna("")
    for partner in partner_rules:
        col = "HasAntibody" if partner["name"] == "antibody" else partner["name"]
        df[col] = [_has_kw(txt.lower(), partner["keywords"]) for txt in entry_desc]
    df["Receptor"] = df.get("NPC1", False) | df.get("TIM-1", False)
    multi_protein = pd.Series(
        [len(matched_proteins(txt, protein_rules)) > 1 for txt in entry_desc],
        index=df.index,
    )
    df["Complex"] = multi_protein | df["HasAntibody"] | df["Receptor"]

    df.to_csv(args.output_csv, index=False)

    print(f"# Annotated {len(df)} entities -> {args.output_csv}")
    print(f"# Unassigned protein: {(df.Protein == 'Other').sum()}   "
          f"Unknown species: {(df.Code == 'Unknown').sum()}   "
          f"Chimeric: {df.Chimeric.sum()}")
    print("\n=== Protein x Species (entities) ===")
    prot = df[df.Protein.isin([r["name"] for r in protein_rules])]
    print(pd.crosstab(prot.Protein, prot.Code, margins=True).to_string())
    print("\n=== Structures with antibody / receptor complexes ===")
    print(f"antibody: {df[df.HasAntibody].PDB.nunique()} structures   "
          f"receptor(NPC1/TIM-1): {df[df.Receptor].PDB.nunique()} structures")


if __name__ == "__main__":
    main()
