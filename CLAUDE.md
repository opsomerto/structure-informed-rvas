# Structure-informed rare-variant association (SIR) — orientation

Structure-informed RVAS: testing whether rare variants in neuropsychiatric disease
are enriched at **protein–protein interface residues**, using PIONEER's
partner-specific interface predictions.

Collaboration: Sherif Gerges (Daly lab, Stanley Center/Broad) with Hilary Finucane's
group; `annotation_test.py` is Emily Nason's code.

---

## CRITICAL: use the right branch

This repo is on branch **`annotation-test-updated`** (commit `0d327f5`), which fixes an
off-by-one bug in the neighborhood expansion. `expand_annot_neighborhood` previously did:

```python
return np.where(is_neighbor > 0)[0]        # 0-based indices used as 1-based aa_pos
return np.where(is_neighbor > 0)[0] + 1    # FIXED
```

Any result generated before this fix is wrong for radius > 0 runs. It does **not** affect
the global/overall enrichment test (no neighborhood) or PyMOL structural distances.

Ignore `annotation_test.py.orig` and `annotation_test.py.orig.2` — stale copies.

---

## The three repos

| path | role | git |
|---|---|---|
| `GITHUB/sir-annotation` | per-interface Fisher test, results, PyMOL mapping | yes |
| `GITHUB/ppi_pioneer` | PIONEER download, annotation building, global enrichment test | no |
| `GITHUB/structure-informed-rvas` | upstream pipeline + reference structures | yes |

**Reference data** lives at `../structure-informed-rvas/sir-reference-data/`:
`pdb_files/` holds 21,943 gzipped AlphaFold monomer models (`AF-<uniprot>-F1-model_v6.pdb.gz`)
and `pae_files/` the matching PAE JSONs. **Do not bulk-read these.**

---

## Key files

**Analysis**
- `annotation_test.py` — entry point `annotation_test(df_rvas, annotation_file, annotation_id,
  reference_dir, neighborhood_radius, pae_cutoff, results_dir, df_filter, ignore_ac)`.
  Per-interface 2×2 Fisher (variants on the interface vs rest of the same protein), then BH-FDR
  pooled across all tested interfaces.
- `debug_annotation_test.py` — heavily-printed 7-step walkthrough of one protein. Best
  orientation for a new reader.
- `compute_asd_medium_r5_fixed.py` / `run_asd_medium_r5_fixed.py` — the corrected ASD re-run.
- `../ppi_pioneer/overall_enrichment_test.py` — pooled global test (radius 0, single Fisher,
  no FDR).

**Data**
- `ASD_mapped.tsv` — ASD variants: case = de novo in probands, control = untransmitted
  parental. Columns `uniprot_id, aa_pos, aa_ref, aa_alt, Variant ID, ac_case, ac_control`.
  753,447 rows / 19,150 proteins / 28,435 case alleles / 726,339 control alleles.
  NOTE `ppi_pioneer/` contains a *different* file of the same name — check which a script uses.
- `pioneer_human_{very_high,high,medium}_annotation.tsv` — `uniprot_id, aa_pos, annotation_id`
  where `annotation_id` is `<focal>_<partner>`. Tiers are cumulative (medium ⊃ high ⊃ very_high).
- `../ppi_pioneer/annotations/pioneer_human_medium_interfaces_long.tsv` — adds the **`source`**
  column: `PDB` (experimental co-structure) / `HM` (homology model) / `PIONEER` (pure prediction).
  Essential for judging how much to trust an interface.

**Results**
- `results/ASD_ppi_medium_r5_FIXED/asd_medium_r5_FIXED.results.tsv` — the current headline run.
  6,040 interfaces tested, 174 FDR-significant, 13 focal genes.
- `results/ASD_ppi_medium_r5_FIXED/sig_genes_all_partners.tsv` — tidy per-partner table for plotting.
- `results/PIONEER_PDB_structures.md` — which hits have experimental structures.
- `results/ASD_ppi_r5/pmls/*.py` — PyMOL mapping scripts, one per structure.

---

## Conventions

**PyMOL colour scheme** (keep consistent):
purple = interface residue only; orange = interface + case variant; red = case variant
off-interface; deepblue = control variant. Never rescale spheres (`sphere_scale` stays 1).

**Numbering.** Deposited PDB numbering often differs from UniProt (engineered constructs,
deleted loops, cleaved signal peptides). Every mapping script aligns the observed sequence
to the UniProt canonical with `difflib.SequenceMatcher` rather than assuming an offset.
Always verify with a spot-check before trusting a mapping.

**`ignore_ac=True`** binarises counts and drops variants shared between cases and controls
(i.e. case-only vs control-only).

---

## Findings worth knowing

- **The de novo class drives the signal, not the disorder.** Global enrichment: ASD de novo
  OR 1.28 (p=1.4e-33), SCZ de novo OR 1.30 (p=0.008), SCZ inherited/case-control OR 1.00
  (null). Not a hub artifact.
- **Most significant interfaces have no structural backing**: of the 174, only 12 (7%) are
  `PDB`-sourced, 16 `HM`, 146 `PIONEER`.
- **Partner multiplication inflates counts.** PIONEER predicts the same surface against many
  partners, so one finding appears as many interfaces (TCF4: 133 significant, only 1 PDB-backed;
  GABRA1's β2/β3 pair shares 135 of ~150 residues). Count genes, not interfaces.
- **The validated hits are known disease genes** — SMAD4's R496/I500 hotspot was published as
  the Myhre syndrome mechanism in 2012. Treat structural confirmations as positive controls
  for the method, not as new biology.
- **Watch direction.** H3C1 has 19 "significant" interfaces but 18 are *control*-enriched
  (OR<1) — filter `or > 1` before reporting counts.
- **Failed validation**: CNOT3–CREBBP (AF3 ipTM ≤ 0.36 across 4 runs; PIONEER's predicted site
  contradicts the only experimental mapping), DZIP1–NUS1 (interaction itself unvalidated),
  MECP2 (84% disordered, no rigid interface to model).

---

## Gotchas

- `fdr_reject` reads as the string `"True"`/`"False"` — coerce before filtering
  (in R: `fdr_reject %in% c("True","TRUE",TRUE)`).
- Rows with `p == NaN` are untested (no variants on one side) — filter them out.
- Proteins >2,700 aa are split into multiple AlphaFold fragments (F1, F2…) and need the
  guide-based stitching; ~56 giant proteins (titin etc.) are excluded outright.
- Background jobs do not survive between shell calls in this environment — use chunked,
  resumable foreground runs.
