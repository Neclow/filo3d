"""Write per-entity amino-acid FASTA for the filovirus protein set.

The canonical one-letter sequence is carried in the metadata (fetched from the
RCSB Data API), so this stage just filters to protein entities and writes a
FASTA whose headers embed the useful database fields:

    >{entity_id}|{Protein}|{Code}|res={resolution}|chains={auth_chains}

One record per polymer entity (a homo-oligomer is one entity / one sequence).
Non-protein entities (RNA) and rows without a sequence are skipped.
"""

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src._config import load_params  # noqa: E402

PROTEINS = {"NP", "VP35", "VP40", "GP", "VP30", "VP24", "L"}


def parse_args():
    parser = ArgumentParser(
        description="Write AA FASTA for filovirus protein entities.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("params")
    parser.add_argument("input_csv", help="Master table clean.csv (or aug.csv)")
    parser.add_argument("output_fasta")
    parser.add_argument(
        "--all", action="store_true",
        help="Keep every polymer entity, not just canonical filovirus proteins.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    load_params(args.params)  # validates params file is present/consistent
    df = pd.read_csv(args.input_csv)

    keep = df["polymer_type"] == "Protein"
    if not args.all:
        keep &= df["Protein"].isin(PROTEINS)
    keep &= df["sequence"].fillna("").str.len() > 0
    sub = df[keep]

    out = Path(args.output_fasta)
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out, "w", encoding="utf-8") as fh:
        for _, r in sub.iterrows():
            res = "" if pd.isna(r.get("resolution")) else r["resolution"]
            header = (
                f">{r['entity_id']}|{r['Protein']}|{r['Code']}"
                f"|res={res}|chains={r.get('auth_chains', '')}"
            )
            seq = str(r["sequence"])
            fh.write(header + "\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
            n += 1

    print(f"# Wrote {n} protein sequences -> {out}")
    by_prot = sub["Protein"].value_counts()
    print("# By protein: " + ", ".join(f"{k}={v}" for k, v in by_prot.items()))


if __name__ == "__main__":
    main()
