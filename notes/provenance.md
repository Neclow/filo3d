# Filo3D data provenance

## Acquisition method

- **Source:** RCSB PDB, via the Search API (`search.rcsb.org/rcsbsearch/v2`) +
  Data API GraphQL (`data.rcsb.org/graphql`).
- **Query:** all *experimental* polymer entities whose source-organism
  taxonomy lineage contains **Filoviridae** (`rcsb_entity_source_organism.taxonomy_lineage.id == 11266`).
  Using the lineage id captures every genus/species below the family in one query.
- **Return type:** `polymer_entity` (one hit per protein chain-group).
- **Classification:** rule-based, from RCSB entity descriptions + source organism
  (`data/metadata/proteins.json`, `data/metadata/taxonomy.json`). Whole-token
  keyword matching. Taxonomy backbone follows Carroll et al. 2013 (`papers/carroll2013.pdf`).

## Snapshot — 2026-09-28

- Structures: **237**  |  polymer entities: **342**
- Antibody-complex structures: **66**  |  receptor complexes (NPC1/TIM-1): **2**
- Chimeric entities: **11** (VP35–NP fusions, MBP-tagged L/VP35 constructs)
- wwPDB validation coverage: **340/342** entities
- Methods (by structure): X-ray = 164, cryo-EM = 64, solution NMR = 8

### Protein × species (polymer entities)

```
Code     BDBV  EBOV  LLOV  MARV  MLAV  RAVV  RESTV  SUDV  TAFV  All
Protein
GP          9   132     0     9     2    10      0    12     0  174
L           0     7     0     2     0     0      0     1     0   10
NP          2    16     2     8     7     0      1     3     1   40
VP24        0     6     0     1     0     0      1     2     0   10
VP30        0     9     1     1     2     0      1     0     0   14
VP35        0    34     0    10     0     0      6     2     0   52
VP40        1    13     0     1     0     0      0    18     0   33
All        12   217     3    32    11    10      9    38     1  333
```
(8 RNA entities and 1 poly-Ala fragment are excluded from the protein crosstab.)

## Completeness — sequence-similarity backstop (`search_similar`)

The taxonomy query only finds structures RCSB annotates with a Filoviridae source
organism. To measure what it misses, an RCSB **sequence** search (mmseqs2) is run
against curated, tag-free reference sequences (`data/metadata/references.fa`, one
per protein; GP reference from Neil, ectodomain, no fibritin/foldon tag). Hits
not already in the taxonomy set become a review queue (`candidates.csv`), not
auto-added.

**Critical caveat learned:** reference sequences MUST be tag-free. Filovirus GP is
deposited as fibritin/foldon fusions and L/VP35 as MBP fusions; searching with
tagged sequences returns every MBP / T4-fibritin / unrelated fusion in the PDB
(~1,900 false hits at 90% identity). Tag-free references + a reference-coverage
filter remove that noise.

Findings (2026-09-28):
- **identity ≥ 0.5, e ≤ 0.1, coverage ≥ 0.3:** 178 union hits, **0 new** — the
  taxonomy set is complete at the 50% level.
- **identity ≥ 0.3 (same filters):** exactly **1 new** candidate, no noise:
  - `6DKU_1` — VP35 interferon-inhibitory domain from *Myotis lucifugus*
    (little brown bat), 36.9% identity, full coverage. Annotated under the **bat
    host organism** (taxid 59463), so invisible to the taxonomy query. A genuine
    borderline filovirus-related entry for manual review.

## Known limitations / TODO

- `Orthoebolavirus sp.` entries default to EBOV via the name fallback — genus is
  correct, species is a best guess; revisit if strain-level resolution matters.
- Species is assigned from the source organism, not from sequence phylogeny.
- Antibody/receptor flags are inferred from sibling entity descriptions in the
  same PDB entry (entry-level context), not from interface analysis.
