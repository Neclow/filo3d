"""Download CIF structures and attach wwPDB validation quality metrics.

Input is the annotated per-entity ``aug.csv``. This stage:

1. downloads one mmCIF file per unique PDB entry (via biotite / RCSB), and
2. fetches wwPDB validation-report percentiles per entry (clashscore, Ramachandran
   and rotamer outliers, and their global percentiles),

then writes ``clean.csv`` = the per-entity metadata joined with per-structure
quality metrics and a ``Downloaded`` flag. This is the master table of the
Filo3D data backbone.

The validation fetch is ported from the BCov3D pipeline (Neclow).
"""

import gzip
import multiprocessing
import sys
import urllib.request
import xml.etree.ElementTree as ET
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from functools import partial
from pathlib import Path

import pandas as pd
from biotite.database import rcsb
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src._config import load_params  # noqa: E402

VALIDATION_URL = (
    "https://files.wwpdb.org/pub/pdb/validation_reports/{mid}/{pdb}/{pdb}_validation.xml.gz"
)

VALIDATION_FIELDS = {
    "absolute-percentile-clashscore": "percentile_clashscore",
    "relative-percentile-clashscore": "percentile_clashscore_relative",
    "absolute-percentile-percent-rama-outliers": "percentile_rama",
    "relative-percentile-percent-rama-outliers": "percentile_rama_relative",
    "absolute-percentile-percent-rota-outliers": "percentile_rota",
    "relative-percentile-percent-rota-outliers": "percentile_rota_relative",
    "clashscore": "clashscore",
    "percent-rama-outliers": "rama_outliers",
    "percent-rota-outliers": "rota_outliers",
}


def fetch_validation(pdb_id):
    """Fetch wwPDB validation percentiles for a single PDB entry."""
    pdb = pdb_id.lower()
    mid = pdb[1:3]
    url = VALIDATION_URL.format(mid=mid, pdb=pdb)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            xml_bytes = gzip.decompress(resp.read())
        root = ET.fromstring(xml_bytes)
        entry = root.find("Entry")
        if entry is None:
            return pdb_id, {}
        return pdb_id, {
            col: float(entry.attrib[attr])
            for attr, col in VALIDATION_FIELDS.items()
            if attr in entry.attrib
        }
    except Exception:
        return pdb_id, {}


def parse_args():
    parser = ArgumentParser(
        description="Download filovirus CIFs and attach wwPDB validation metrics.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("params")
    parser.add_argument("input_csv", help="Annotated aug.csv (per-entity)")
    parser.add_argument("output_dir", help="Directory for downloaded CIF files")
    parser.add_argument("output_csv", help="Output master table clean.csv")
    parser.add_argument("--fmt", default="cif", help="Download format")
    return parser.parse_args()


def main():
    args = parse_args()
    params = load_params(args.params)
    threads = params.get("threads", 16)

    df = pd.read_csv(args.input_csv)
    pdb_ids = sorted(df["PDB"].dropna().astype(str).str.lower().unique())
    print(f"# Structures to download: {len(pdb_ids)} ({len(df)} entities)")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Download CIFs in parallel.
    with multiprocessing.Pool(processes=threads) as pool:
        fetch_fn = partial(rcsb.fetch, format=args.fmt, target_path=str(out_dir), verbose=False)
        for _ in tqdm(
            pool.imap_unordered(fetch_fn, pdb_ids), total=len(pdb_ids), desc="Downloading CIF"
        ):
            pass

    downloaded = {p for p in pdb_ids if (out_dir / f"{p}.cif").exists()}
    print(f"# Downloaded: {len(downloaded)}/{len(pdb_ids)}")

    # Fetch wwPDB validation metrics in parallel.
    print("# Fetching wwPDB validation reports...")
    with multiprocessing.Pool(processes=threads) as pool:
        val_results = list(
            tqdm(
                pool.imap_unordered(fetch_validation, pdb_ids),
                total=len(pdb_ids),
                desc="Validation",
            )
        )
    val_df = pd.DataFrame({pid: m for pid, m in val_results if m}).T
    n_val = len(val_df)
    print(f"# Validation metrics for {n_val}/{len(pdb_ids)} structures")

    df["Downloaded"] = df["PDB"].str.lower().isin(downloaded)
    if not val_df.empty:
        val_df.index.name = "PDB"
        df = df.merge(val_df, left_on="PDB", right_index=True, how="left")

    df.to_csv(args.output_csv, index=False)
    print(f"# Wrote master table ({len(df)} rows) -> {args.output_csv}")


if __name__ == "__main__":
    main()
