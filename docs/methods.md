# Data & methods

Filo3D is assembled by a reproducible pipeline (Pixi + DVC), not by hand. See the
[GitHub repository](https://github.com/Neclow/filo3d) for the code.

## Acquisition — two routes

1. **NCBI taxonomy (trusted core).** An RCSB search on source-organism lineage =
   Filoviridae (taxid `11266`) returns every experimental polymer entity below the
   family; metadata is fetched from the RCSB Data API. Entities are classified into
   canonical proteins and species, then CIFs are downloaded and annotated with wwPDB
   validation percentiles.
2. **Sequence-similarity backstop.** An RCSB mmseqs2 search against tag-free reference
   sequences catches structures the taxonomy route misses (e.g. host-annotated or
   synthetic constructs, such as the bat *Myotis* VP35 `6DKU`). These are written to a
   review queue, never auto-added.

Taxonomy backbone follows Carroll et al. 2013.

## Composition (polymer entities)

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

## Notes

- Species codes: EBOV (Zaire), SUDV (Sudan), BDBV (Bundibugyo), RESTV (Reston),
  TAFV (Taï Forest), MARV (Marburg), RAVV (Ravn), LLOV (Lloviu), MLAV (Měnglà).
- One row per polymer entity; a structure may contribute several (e.g. GP1 + GP2, or
  a bound Fab). Sequences carry expression tags as deposited.
- Amino-acid sequences and the full master table are in the repository
  (`data/aa/fa/filo3d.fa`, `data/metadata/clean.csv`).
