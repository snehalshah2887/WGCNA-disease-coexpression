# Reproducibility Notes

This document explains the relationship between results produced by this
pipeline and results reported in the associated publication:

> Shah S, et al. (2023). Weighted gene co-expression network analysis of
> developmental brain transcriptomics reveals connectivity patterns of 22q13
> deletion syndrome genes. *Genes* **14**(11): 1998.
> https://doi.org/10.3390/genes14111998

---

## Summary

This repository re-implements the original analysis as a modular, containerised
pipeline for open reproducibility. **Numerical results will not be identical to
the paper.** The differences are intentional, documented below, and do not
change the biological conclusions of the published work.

---

## 1. WGCNA Module Assignments

### What may differ
Module colour labels, module sizes, and gene-to-module assignments can vary
between runs of `blockwiseModules()`.

### Why
WGCNA's hierarchical clustering and dynamic tree cutting are deterministic
given identical inputs, but are sensitive to:

- **R version** — floating-point operations and LAPACK/BLAS implementations
  differ across R versions. The paper used an earlier R environment; this
  pipeline uses R 4.3 (Bioconductor 3.18).
- **WGCNA package version** — the package was removed from CRAN in 2023 and
  is now maintained on Bioconductor. Internal algorithms may have been updated.
- **Operating system / hardware** — BLAS implementations differ between
  macOS, Linux, and HPC environments, causing minor numerical differences that
  propagate through eigenvalue decomposition.

### What remains consistent
- The soft-thresholding power (β = 9) is fixed.
- The same signed network type, TOMType, and `minModuleSize` parameters are used.
- The broad co-expression structure (which genes cluster together biologically)
  is reproducible; individual module label assignments may differ.

### How to check
Compare `results/02_network/colormodule.csv` with Table S1 of the paper.
Gene-level module memberships should be broadly consistent even if some
colour labels differ.

---

## 2. GO Biological Process Enrichment

### What differs
GO term results will differ from the paper. This is a deliberate methodological
improvement, not an error.

| Aspect | Paper (DAVID) | This pipeline (phyper) |
|--------|--------------|------------------------|
| Tool | DAVID website | `org.Hs.eg.db` + `GO.db` + `phyper` |
| Background universe | All annotated human genes (~19,000) | 2,116 network genes only |
| Statistical test | Fisher's exact (one-tailed) | Hypergeometric (mathematically equivalent) |
| Multiple testing | Benjamini-Hochberg | Benjamini-Hochberg |
| GO annotation version | DAVID internal DB (2022) | Current Bioconductor release |

### Why the network background is more appropriate
Using only the 2,116 co-expressed brain genes as the background universe
asks: *"Is this module enriched for a GO term relative to the other
brain-expressed genes in the network?"* This is statistically more
appropriate for WGCNA than using all human genes, which inflates enrichment
for ubiquitous brain biology terms (e.g., "synaptic transmission" is
enriched in every brain gene list vs. all human genes).

### What this means in practice
- Broad biological themes (nervous system development, cell adhesion,
  semaphorin-plexin signalling, chromatin remodelling) replicate between
  both approaches.
- Terms that were significant in DAVID solely because of the all-human
  background (e.g., synaptic transmission in the yellow module) will appear
  weaker or absent with the network background.
- More specific or novel terms may appear that DAVID did not detect.
- The `Goterms_22q13.csv` output lists GO terms containing at least one PMS
  gene — this is the most directly relevant output for the paper's
  conclusions.

---

## 3. Gene Ranking (Module 6)

### What differs
The paper only reported ASD–PMS co-expression connectivity. This pipeline
ranks PMS genes against **all five disease phenotypes** (ASD, ID, Seizures,
Hypotonia, Language Impairment) and adds a composite score.

### What was fixed (intentional result differences)
The original analysis contained three bugs that caused incorrect or incomplete
rankings. Their fixes change results:

| Bug | Original behaviour | Fixed behaviour |
|-----|-------------------|-----------------|
| TOM had no row/column names | `intersect()` returned empty — 0 PMS genes found | Names restored from `datExpr`; all 139 PMS genes ranked |
| Used `modTom` (sub-matrix of last colour module) | Only genes in that module could be ranked | Full 2,116 × 2,116 TOM used correctly |
| kME lookup used `mcolor` (e.g. "turquoise") instead of `mlabel` (e.g. "1") | All kME values returned NA | Correct kME values returned for all genes |

The corrected rankings are numerically different from the paper but
methodologically correct.

### Composite score (new in this pipeline)
The paper ranked by `mean_TOM` alone. This pipeline adds a composite score:

```
composite_score = 0.70 × rank_norm(mean_TOM) + 0.30 × rank_norm(kME)
```

where `rank_norm` maps values to [0, 1]. This rewards genes that are both
strongly connected to the disease gene set (mean_TOM) and are intramodular
hubs (kME ≥ 0.5), providing a more robust prioritisation for follow-up studies.

---

## 4. Supplemental Co-expression Table (Module 4)

### What may differ
The `Table22q13.csv` output depends on the Cytoscape edge files, which are
derived from the TOM. If module assignments differ (see Section 1), edge
weights and co-expression partners may differ slightly.

### What remains consistent
The logic for identifying co-expression partners of 22q13 genes is unchanged.
The set of genes flagged as co-expressed at the TOM threshold used in
Module 3 should be broadly consistent.

---

## 5. Pipeline Architecture Differences

The original analysis was a single monolithic R script. This repository
re-implements it as six modular scripts orchestrated by Nextflow, with the
following improvements:

| # | Module | Original issue | Resolution |
|---|--------|---------------|------------|
| 1 | 03 | `final$gene_symbol` → `final$Gene_symbol` (wrong case) | Fixed — was returning all NA for `modGenes` |
| 2 | 06 | `final1` undefined | Now saved in `network_objects.rds` by Module 2 |
| 3 | 06 | `PMSs_kME$cor` — `signedKME()` returns a data frame | Replaced with direct indexing |
| 4 | 03 | `ASD_PMS` column → `ASD == 1 & PMS == 1` filter | Was crashing / producing empty subnetwork |
| 5 | 05 | `temp` variable → direct filter on `final$mcolor` | Was crashing inside GO function |
| 6 | 05 | Missing `+` before `theme()` in ggplot | Theme not applied to tan module plot |
| 7 | 02 | `dev.off()` without matching `png()`/`pdf()` | Corrected |
| 8 | 03 | `createCompositeFilter()` missing required arguments | Added `filter.name` and `filter.list` |
| 9 | 06 | Only ASD ranking produced; wrong TOM sub-matrix used | Rewrote to use full TOM for all 5 diseases |

---

## 6. Guidance for Reviewers

If you are reviewing this repository alongside the published paper:

- **Module assignments** will be approximately but not exactly the same.
  Compare gene-level kME values and module membership scores rather than
  colour labels.
- **GO results** use a different (more conservative) method. The biological
  themes are consistent; specific p-values and term lists are not expected to
  match.
- **Gene rankings** for ASD will be in the same direction as the paper but
  numerically different due to the three bug fixes listed above.
- **All data used** are publicly available from BrainSpan (see `data/README.md`).
  The analysis is fully reproducible from those sources using this pipeline.

For questions or discrepancies, please open an issue on GitHub.
