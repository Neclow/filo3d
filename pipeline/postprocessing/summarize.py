"""Compute Filo3D summary tables from the master table.

Produces:
- a **protein x species** crosstab of polymer entities, and
- a one-line **release summary** (structures, entities, complexes, methods),

both as GitHub-flavoured Markdown. Prints to stdout, and can write a standalone
Markdown file (--out) and/or splice the tables into the README between marker
comments (--readme), so the published counts stay in sync with the data.

    pixi run summarize                      # print to stdout
    pixi run summarize --readme README.md   # update README in place
"""

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src._config import load_params  # noqa: E402

PROTEINS = ["GP", "NP", "VP35", "VP40", "VP30", "VP24", "L"]
CROSSTAB_START = "<!-- filo3d:crosstab:start -->"
CROSSTAB_END = "<!-- filo3d:crosstab:end -->"


def crosstab_markdown(df):
    prot = df[df["Protein"].isin(PROTEINS)]
    ct = pd.crosstab(prot["Protein"], prot["Code"], margins=True, margins_name="All")
    ct = ct.reindex([p for p in PROTEINS if p in ct.index] + ["All"])
    cols = list(ct.columns)
    lines = [
        "| Protein | " + " | ".join(cols) + " |",
        "|" + "---|" * (len(cols) + 1),
    ]
    for name, row in ct.iterrows():
        lines.append(f"| **{name}** | " + " | ".join(str(int(v)) for v in row) + " |")
    return "\n".join(lines)


def release_summary(df):
    prot = int((df["polymer_type"] == "Protein").sum())
    rna = int((df["polymer_type"] == "RNA").sum())
    species = [c for c in df["Code"].dropna().unique() if c != "Unknown"]
    ab = df.loc[df["HasAntibody"], "PDB"].nunique()
    rec = df.loc[df["Receptor"], "PDB"].nunique()
    methods = ", ".join(
        f"{k} {v}" for k, v in df.drop_duplicates("PDB")["method"].value_counts().items()
    )
    return "\n".join([
        "| | |",
        "|---|---|",
        f"| Structures | **{df['PDB'].nunique()}** |",
        f"| Polymer entities | **{len(df)}** ({prot} protein, {rna} RNA) |",
        f"| Species | {', '.join(sorted(species))} |",
        f"| Antibody complexes | {ab} structures |",
        f"| Receptor (NPC1/TIM-1) complexes | {rec} structures |",
        f"| Methods | {methods} |",
    ])


def update_readme(path, block):
    text = Path(path).read_text(encoding="utf-8")
    if CROSSTAB_START not in text or CROSSTAB_END not in text:
        raise SystemExit(
            f"README markers not found; add these where the table should go:\n"
            f"  {CROSSTAB_START}\n  {CROSSTAB_END}"
        )
    pre = text.split(CROSSTAB_START)[0]
    post = text.split(CROSSTAB_END)[1]
    new = f"{pre}{CROSSTAB_START}\n{block}\n{CROSSTAB_END}{post}"
    Path(path).write_text(new, encoding="utf-8")
    print(f"# Updated crosstab in {path}")


def parse_args():
    parser = ArgumentParser(
        description="Compute Filo3D summary tables.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("params")
    parser.add_argument("input_csv", help="Master table clean.csv")
    parser.add_argument("--out", help="Write the Markdown tables to this file")
    parser.add_argument("--readme", help="Splice the crosstab into this README (between markers)")
    return parser.parse_args()


def main():
    args = parse_args()
    load_params(args.params)
    df = pd.read_csv(args.input_csv)

    ct = crosstab_markdown(df)
    summary = release_summary(df)
    combined = (
        "### Protein × species (polymer entities)\n\n" + ct + "\n"
    )

    print(summary + "\n\n" + combined)

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(
            "## Release summary\n\n" + summary + "\n\n" + combined, encoding="utf-8"
        )
        print(f"# Wrote {args.out}")
    if args.readme:
        update_readme(args.readme, combined.rstrip())


if __name__ == "__main__":
    main()
