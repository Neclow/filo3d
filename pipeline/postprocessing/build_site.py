"""Build the data file for the Filo3D GitHub Pages browser.

Reads the master table (clean.csv) and writes docs/data.json — a compact,
column-oriented payload (meta + column names + row arrays) consumed by the
static searchable table in docs/index.html. One row per polymer entity.
"""

import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from datetime import date
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src._config import load_params  # noqa: E402

METHOD = {
    "X-RAY DIFFRACTION": "X-ray",
    "ELECTRON MICROSCOPY": "cryo-EM",
    "SOLUTION NMR": "NMR",
}

# Displayed column -> source field (order matters; must match docs/index.html)
COLUMNS = [
    ("PDB", "PDB"),
    ("Protein", "Protein"),
    ("Species", "Code"),
    ("Virus", "Virus"),
    ("Description", "description"),
    ("Method", "method"),
    ("Res (Å)", "resolution"),
    ("Len", "seq_length"),
    ("Chains", "auth_chains"),
    ("Ab", "HasAntibody"),
    ("Receptor", "Receptor"),
    ("Complex", "Complex"),
    ("Chimeric", "Chimeric"),
    ("Clash %ile", "percentile_clashscore"),
    ("Released", "release_date"),
]


def cell(field, value):
    if pd.isna(value):
        return ""
    if field in {"HasAntibody", "Receptor", "Complex", "Chimeric"}:
        return "Yes" if bool(value) else ""
    if field == "method":
        return METHOD.get(str(value), str(value).split(";")[0])
    if field == "resolution":
        try:
            return round(float(value), 2)
        except (TypeError, ValueError):
            return ""
    if field in {"seq_length", "percentile_clashscore"}:
        try:
            return int(round(float(value)))
        except (TypeError, ValueError):
            return ""
    return value


def parse_args():
    parser = ArgumentParser(
        description="Build docs/data.json for the Filo3D Pages browser.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("params")
    parser.add_argument("input_csv", help="Master table clean.csv")
    parser.add_argument("output_json", help="e.g. docs/data.json")
    return parser.parse_args()


def main():
    args = parse_args()
    load_params(args.params)
    df = pd.read_csv(args.input_csv).sort_values(["Protein", "Code", "PDB"])

    fields = [f for _, f in COLUMNS]
    rows = [[cell(f, r[f]) for f in fields] for _, r in df.iterrows()]

    prot = df[df["polymer_type"] == "Protein"]
    payload = {
        "meta": {
            "generated": date.today().isoformat(),
            "n_structures": int(df["PDB"].nunique()),
            "n_entities": int(len(df)),
            "n_protein": int(len(prot)),
            "species": sorted(c for c in df["Code"].dropna().unique() if c != "Unknown"),
        },
        "columns": [name for name, _ in COLUMNS],
        "rows": rows,
    }

    out = Path(args.output_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"# Wrote {len(rows)} rows -> {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
