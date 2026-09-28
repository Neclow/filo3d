# Filo3D

A curated, quality-scored database of high-resolution **Filoviridae** protein
structures from the [RCSB PDB](https://www.rcsb.org) — a filovirus counterpart to
[CoV3D](https://cov3d.ibbr.umd.edu)
([Gowthaman et al., *NAR* 2021](https://doi.org/10.1093/nar/gkaa731)).

Filo3D covers the full filovirus proteome — **NP, VP35, VP40, GP, VP30, VP24, L** —
across *Ebolavirus*, *Marburgvirus*, *Cuevavirus* (Lloviu) and *Dianlovirus*
(Měnglà), with amino-acid sequences, structure quality metrics, and antibody /
receptor complex annotation. The dataset is assembled by a reproducible pipeline,
not by hand, so it can be rebuilt or refreshed with one command.

## Current release (2026-09-28)

| | |
|---|---|
| Structures | **237** |
| Polymer entities | **342** (333 protein, 8 RNA) |
| Species | EBOV, SUDV, BDBV, RESTV, TAFV, MARV, RAVV, LLOV, MLAV |
| Antibody complexes | 66 structures |
| Receptor (NPC1) complexes | 2 structures |
| Methods | X-ray 164, cryo-EM 64, NMR 8 |

### Protein × species (polymer entities)

| Protein | BDBV | EBOV | LLOV | MARV | MLAV | RAVV | RESTV | SUDV | TAFV | All |
|---|---|---|---|---|---|---|---|---|---|---|
| **GP** | 9 | 132 | 0 | 9 | 2 | 10 | 0 | 12 | 0 | 174 |
| **NP** | 2 | 16 | 2 | 8 | 7 | 0 | 1 | 3 | 1 | 40 |
| **VP35** | 0 | 34 | 0 | 10 | 0 | 0 | 6 | 2 | 0 | 52 |
| **VP40** | 1 | 13 | 0 | 1 | 0 | 0 | 0 | 18 | 0 | 33 |
| **VP30** | 0 | 9 | 1 | 1 | 2 | 0 | 1 | 0 | 0 | 14 |
| **VP24** | 0 | 6 | 0 | 1 | 0 | 0 | 1 | 2 | 0 | 10 |
| **L** | 0 | 7 | 0 | 2 | 0 | 0 | 0 | 1 | 0 | 10 |
| **All** | 12 | 217 | 3 | 32 | 11 | 10 | 9 | 38 | 1 | 333 |

Full provenance and the sequence-similarity completeness check: [`notes/provenance.md`](notes/provenance.md).

## Data

| File | Contents |
|---|---|
| `data/metadata/clean.csv` | Master table: one row per entity — protein, species, sequence, resolution, method, wwPDB validation percentiles, complex/chimeric flags |
| `data/metadata/aug.csv` | Annotated metadata (pre-download) |
| `data/aa/fa/filo3d.fa` | Amino-acid sequences (metadata-rich headers) |
| `data/cif/` | mmCIF structures — not in git; run `pixi run download` |

## Reproduce

```bash
pixi install      # dependencies
pixi run repro    # full pipeline (search → annotate → download → extract_aa)
```

Individual stages: `pixi run search | annotate | download | extract_aa | search_similar`.

## How it's built

Two acquisition routes feed the database:

1. **Taxonomy (trusted core)** — RCSB search on source-organism lineage =
   Filoviridae (NCBI taxid `11266`), then GraphQL metadata fetch, protein/species
   classification, CIF download + validation metrics, and AA extraction.
2. **Sequence similarity (completeness backstop)** — `search_similar` runs an RCSB
   mmseqs2 search against tag-free reference sequences (`data/metadata/references.fa`)
   to catch structures the taxonomy route misses (e.g. host-annotated or synthetic
   constructs). Hits are written to a review queue, never auto-added.

Classification rules (`proteins.json`, `taxonomy.json`) and all thresholds
(`params.yaml`) are editable. Taxonomy backbone follows Carroll et al. 2013.

## Citation

If you use Filo3D, please cite this repository. The design follows CoV3D
(Gowthaman et al., *Nucleic Acids Research* 2021).
