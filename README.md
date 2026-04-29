# 22q13 / Phelan-McDermid Syndrome — WGCNA Co-expression Network Analysis

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A523.04.0-brightgreen)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-24.0%2B-blue)](https://www.docker.com/)
[![R](https://img.shields.io/badge/R-4.3.x-276DC3)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18-87B13F)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A reproducible, containerised Nextflow + Docker pipeline for Weighted Gene
Co-expression Network Analysis (WGCNA) applied to 22q13 deletion syndrome
(Phelan-McDermid Syndrome, PMS) using BrainSpan developmental brain
transcriptomics data.

> **Associated publication**
> Shah S, et al. (2023). Weighted gene co-expression network analysis of
> developmental brain transcriptomics reveals connectivity patterns of 22q13
> deletion syndrome genes. *Genes* **14**(11): 1998.
> https://doi.org/10.3390/genes14111998

> **Why this repository matters**
> The WGCNA R package was removed from CRAN in 2023. It is now maintained
> exclusively on Bioconductor. This pipeline provides a fully containerised
> environment so that the analysis can be reproduced without manually
> resolving package availability or environment issues.

> **Reproducibility disclaimer**
> This pipeline has been re-implemented and improved for open reproducibility.
> Numerical results (module assignments, GO terms, gene rankings) may differ
> from those reported in the associated publication. These differences are
> expected and are documented in detail in
> [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). The biological conclusions of
> the paper are not affected.

---

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [Input Data](#input-data)
- [Pipeline Architecture](#pipeline-architecture)
- [Running the Pipeline](#running-the-pipeline)
- [Output Structure](#output-structure)
- [Cytoscape Visualisation](#cytoscape-visualisation)
- [Running Individual Modules](#running-individual-modules)
- [HPC / SLURM](#hpc--slurm)
- [Reproducibility Notes](#reproducibility-notes)
- [Citation](#citation)
- [License](#license)

---

## Overview

This pipeline analyses gene co-expression patterns in the developing human
brain (BrainSpan dataset) for 22q13-region genes and five clinical phenotypes:
Autism Spectrum Disorder (ASD), Intellectual Disability (ID), Seizures,
Hypotonia, and Language Impairment. It:

1. Validates and passes through the pre-filtered expression matrix
2. Constructs a signed, weighted co-expression network (soft-threshold power = 9)
3. Assigns genes to colour-labelled modules and tests disease enrichment (Fisher exact, FDR)
4. Computes the full Topological Overlap Matrix (TOM) and exports Cytoscape network files
5. Builds a co-expression table for 22q13 genes across all modules
6. Runs GO Biological Process enrichment using `org.Hs.eg.db` + hypergeometric test
7. Ranks all 22q13 genes by composite co-expression connectivity to each phenotype's gene set

---

## Prerequisites

### Option A — Docker + Nextflow (recommended)

| Tool     | Minimum version | Install |
|----------|----------------|---------|
| Docker   | 24.0           | https://docs.docker.com/get-docker/ |
| Nextflow | 23.04.0        | `curl -s https://get.nextflow.io \| bash` |
| Java     | 11             | Required by Nextflow |

### Option B — Local R execution

| Package        | Source       | Version tested |
|----------------|--------------|----------------|
| R              | —            | 4.3.x          |
| WGCNA          | Bioconductor | 1.72           |
| org.Hs.eg.db   | Bioconductor | 3.17.0         |
| GO.db          | Bioconductor | 3.17.0         |
| AnnotationDbi  | Bioconductor | 1.62.x         |
| RCy3           | Bioconductor | 2.22           |
| dplyr          | CRAN         | 1.1.x          |
| ggplot2        | CRAN         | 3.4.x          |
| tidyr          | CRAN         | 1.3.x          |
| matrixStats    | CRAN         | 1.0.x          |
| forcats        | CRAN         | 1.0.x          |
| stringr        | CRAN         | 1.5.x          |
| openxlsx       | CRAN         | 4.2.x          |

### Optional — Cytoscape network visualisation

| Tool      | Minimum version | Install |
|-----------|----------------|---------|
| Cytoscape | 3.9.0          | https://cytoscape.org/download.html |

---

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/22q13_WGCNA.git
cd 22q13_WGCNA

# 2. Prepare input data (see data/README.md for full details)
#    Place gene_expression_filtered.csv in data/input/

# 3. Build the Docker image (once)
docker build -t 22q13-wgcna:latest -f docker/Dockerfile .

# 4. Run the full pipeline
nextflow run nextflow/main.nf -params-file conf/params.yaml -profile docker
```

To resume after a partial run (avoids re-running completed steps):

```bash
nextflow run nextflow/main.nf -params-file conf/params.yaml -profile docker -resume
```

---

## Input Data

The pipeline requires one input file:

```
data/input/gene_expression_filtered.csv
```

This is a combined gene annotation + BrainSpan RNA-seq expression matrix:
- Columns 1–10: gene annotation (Ensembl ID, gene symbol, disease flags, 22q13 flag)
- Columns 11–534: RPKM expression values across 524 BrainSpan developmental brain samples
- Pre-filtered to 2,116 genes with mean expression ≥ 0.3

See [`data/README.md`](data/README.md) for detailed format specification and
instructions on how to obtain and prepare the data from public sources.

> **Data included**: `gene_expression_filtered.csv` (11 MB) is distributed
> with this repository. It is derived from the publicly available BrainSpan
> atlas; the gene annotation list was published in the supplemental material
> of the associated paper. See `data/README.md` for full provenance and
> instructions to reproduce the file from raw BrainSpan downloads.

---

## Pipeline Architecture

```
 data/input/gene_expression_filtered.csv
              │
 ┌────────────▼─────────────────────────────────────────────┐
 │ Module 1: Preprocessing                                   │
 │   Validates input, extracts gene annotation block         │
 │   Script: scripts/01_preprocess.R                        │
 └────────────┬─────────────────────────────────────────────┘
              │ gene_expression_filtered.csv + gene_annotation.rds
 ┌────────────▼─────────────────────────────────────────────┐
 │ Module 2: Network Construction                   [16 GB] │
 │   Soft-threshold selection · blockwiseModules            │
 │   Module–disease Fisher enrichment · heatmap             │
 │   Script: scripts/02_build_network.R                     │
 └────────────┬─────────────────────────────────────────────┘
              │ network_objects.rds
    ┌─────────┼──────────────┬──────────────┐
    │         │              │              │
 ┌──▼──┐  ┌──▼──┐       ┌───▼───┐      ┌───▼───┐
 │ M3  │  │ M4  │       │  M5   │      │  M6   │
 │TOM +│  │Supp.│       │  GO   │      │ Gene  │
 │Cyto │  │Table│       │ Enrich│      │Ranking│
 │[16GB│  │     │       │(phyper│      │(5 dis.│
 │]    │  │     │       │)      │      │)      │
 └─────┘  └─────┘       └───────┘      └───────┘
```

| Module | Script | Label | Key output |
|--------|--------|-------|-----------|
| 1 | `01_preprocess.R` | small | `gene_annotation.rds` |
| 2 | `02_build_network.R` | **high_mem** | `network_objects.rds` |
| 3 | `03_export_network.R` | **high_mem** | `tom_objects.rds`, Cytoscape files |
| 4 | `04_supplement_table.R` | small | `Table22q13.csv` |
| 5 | `05_go_analysis.R` | small | GO enrichment CSVs + bar plots |
| 6 | `06_gene_ranking.R` | medium | Ranked PMS gene CSVs + top-5 table |

Modules 3, 4, 5, and 6 run in parallel after Module 2 completes.

---

## Running the Pipeline

### Full pipeline with Docker (recommended)

```bash
nextflow run nextflow/main.nf -params-file conf/params.yaml -profile docker
```

### Custom parameters on the command line

```bash
nextflow run nextflow/main.nf \
    --expression data/input/gene_expression_filtered.csv \
    --outdir     results/ \
    -profile docker
```

### With live Cytoscape import

Requires Cytoscape 3.9+ running on your local machine. Must use `-profile local`
because Docker containers cannot reach the host Cytoscape instance.

```bash
nextflow run nextflow/main.nf \
    --expression     data/input/gene_expression_filtered.csv \
    --skip_cytoscape false \
    -profile local
```

### Available profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run all processes inside the `22q13-wgcna` Docker container |
| `singularity` | Run with Singularity (for HPC environments) |
| `local` | Run without a container (all R packages must be installed locally) |
| `slurm` | Submit jobs to a SLURM cluster via Singularity |

### Resource allocation

Modules 2 and 3 (WGCNA network construction and TOM computation) are
memory-intensive. Default allocations:

| Label | CPUs | Memory | Wall time |
|-------|------|--------|-----------|
| `small` | 2 | 8 GB | 1 h |
| `medium` | 4 | 16 GB | 4 h |
| `high_mem` | 4 | 24 GB (× attempt) | 12 h (× attempt) |

On OOM exit (Linux kill signal 137 or SLURM exit 140), `high_mem` processes
automatically retry with doubled memory and extended wall time (up to 2
retries, max 48 GB / 36 h).

Adjust limits in `nextflow/nextflow.config` for your infrastructure.

---

## Output Structure

```
results/
├── 01_preprocessed/
│   ├── gene_expression_filtered.csv   ← Validated expression matrix
│   └── gene_annotation.rds            ← Gene annotation R object
│
├── 02_network/
│   ├── network_objects.rds            ← All WGCNA R objects
│   ├── colormodule.csv                ← Gene → module colour
│   ├── summary_colors.csv             ← Disease gene counts per module
│   ├── color_modules/
│   │   ├── <color>_module_V1.csv      ← Module gene lists (WGCNA labels)
│   │   └── <color>_module_V2.csv      ← Module gene lists (annotated)
│   ├── TOM/
│   │   └── datExpr-block.1.RData      ← Block TOM (large file)
│   └── plots/
│       ├── sampleClustering.png
│       ├── soft_threshold.png
│       ├── module_dendrogram.png
│       ├── disease_distribution_barplot.png
│       └── enrichment_heatmap.png
│
├── 03_network_export/
│   ├── tom_objects.rds                ← Full TOM + disease gene table
│   └── cytoscape_files/
│       ├── CytoscapeInput-edges-<color>.txt
│       ├── CytoscapeInput-nodes-<color>.txt
│       └── CytoscapeInput-edges-10_29_22.txt   ← Consolidated edge file
│
├── 04_supplement_table/
│   └── Table22q13.csv                 ← Co-expression partners for 22q13 genes
│
├── 05_go_analysis/
│   ├── Goterms_22q13.csv              ← GO BP terms containing ≥1 PMS gene
│   ├── top5_go_all.csv                ← Top 5 GO terms per module
│   ├── filtered_go/
│   │   └── <color>_filtered_module.csv
│   └── plots/
│       ├── go_barplot_<color>.png
│       └── go_barplot_all_modules.png
│
└── 06_gene_ranking/
    ├── ranked_ASD_PMS.csv             ← 139 PMS genes ranked vs ASD
    ├── ranked_ID_PMS.csv              ← 139 PMS genes ranked vs ID
    ├── ranked_Seizures_PMS.csv        ← 139 PMS genes ranked vs Seizures
    ├── ranked_Hypotonia_PMS.csv       ← 139 PMS genes ranked vs Hypotonia
    ├── ranked_LangImp_PMS.csv         ← 139 PMS genes ranked vs Language Impairment
    ├── ranked_all_diseases_PMS.csv    ← Wide-format summary (all 5 diseases)
    └── top5_PMS_disease_genes.csv     ← Top-5 hub genes per disease (publication table)
```

### Gene ranking output columns

Each `ranked_<disease>_PMS.csv` contains:

| Column | Description |
|--------|-------------|
| `rank` | Rank by composite score (1 = highest) |
| `Ensembl` | Ensembl gene ID |
| `Gene_symbol` | Gene name |
| `mcolor` | WGCNA module colour |
| `mlabel` | WGCNA module numeric label |
| `composite_score` | 0.70 × rank_norm(mean_TOM) + 0.30 × rank_norm(kME) |
| `mean_TOM` | Mean TOM connectivity to disease gene set (post q75 filter) |
| `sd_TOM` | Standard deviation of TOM connectivity |
| `kME_own_module` | Intramodular hub score (signed kME to own module eigengene) |
| `hub_gene` | TRUE if kME ≥ 0.5 |
| `stable_CV` | TRUE if sd_TOM / mean_TOM < 1 |
| `n_disease_genes` | Number of disease genes used in ranking |
| `q75_threshold` | 75th-percentile TOM threshold applied |

---

## Cytoscape Visualisation

Module 3 writes Cytoscape-compatible edge and node files for every colour module.

**Manual import:**
1. Open Cytoscape 3.9+
2. File → Import → Network from File → select `CytoscapeInput-edges-<color>.txt`
3. File → Import → Table from File → select `CytoscapeInput-nodes-<color>.txt`
4. Use the VizMapper to colour nodes by `group` (= module colour)

**Automated import (live Cytoscape):**
Re-run with `--skip_cytoscape false -profile local`. Docker cannot reach a host
Cytoscape instance without additional network bridging configuration.

---

## Running Individual Modules

Each script is self-contained and can be run independently:

```bash
# Module 1 — Preprocessing
Rscript scripts/01_preprocess.R \
    data/input/gene_expression_filtered.csv \
    results/01_preprocessed/

# Module 2 — Network construction (requires ~16 GB RAM)
Rscript scripts/02_build_network.R \
    results/01_preprocessed/gene_expression_filtered.csv \
    results/01_preprocessed/gene_annotation.rds \
    results/02_network/

# Module 3 — TOM + Cytoscape export (requires ~16 GB RAM)
Rscript scripts/03_export_network.R \
    results/02_network/network_objects.rds \
    results/03_network_export/ \
    TRUE

# Module 4 — Supplemental co-expression table
Rscript scripts/04_supplement_table.R \
    results/03_network_export/tom_objects.rds \
    results/02_network/network_objects.rds \
    results/03_network_export/cytoscape_files/ \
    results/04_supplement_table/

# Module 5 — GO enrichment (AnnotationDbi + phyper)
Rscript scripts/05_go_analysis.R \
    results/02_network/network_objects.rds \
    results/05_go_analysis/

# Module 6 — Gene ranking (5 disease comparisons)
Rscript scripts/06_gene_ranking.R \
    results/02_network/network_objects.rds \
    results/03_network_export/tom_objects.rds \
    results/06_gene_ranking/
```

---

## HPC / SLURM

Edit `nextflow/nextflow.config` to set your queue, partition, and resource
limits, then:

```bash
nextflow run nextflow/main.nf \
    -params-file conf/params.yaml \
    -profile slurm
```

The pipeline is configured for Singularity under SLURM. Convert the Docker
image to Singularity format with:

```bash
singularity pull 22q13-wgcna.sif docker://YOUR_DOCKERHUB_USERNAME/22q13-wgcna:latest
```

---

## Reproducibility Notes

Results produced by this pipeline may differ from the published paper. This
is expected. See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for a full
explanation of what differs, why, and what remains consistent with the paper.

Key differences at a glance:

| Analysis | Paper | This pipeline | Expected to match? |
|----------|-------|---------------|--------------------|
| WGCNA module assignments | Original run | Re-run (same data, same method) | Approximately — minor differences normal |
| GO enrichment tool | DAVID (whole-genome background) | `org.Hs.eg.db` + `phyper` (network background) | No — different method by design |
| Gene ranking | ASD only | All 5 disease phenotypes | Expanded scope |
| Bug fixes | Several bugs present | All fixed | Results intentionally differ where bugs were fixed |

---

## Citation

See [`CITATION.md`](CITATION.md) for complete citation details.

If you use this pipeline, please cite:

> Shah S, et al. (2023). Weighted gene co-expression network analysis of
> developmental brain transcriptomics reveals connectivity patterns of 22q13
> deletion syndrome genes. *Genes* **14**(11): 1998.
> https://doi.org/10.3390/genes14111998

---

## License

This project is licensed under the MIT License — see [`LICENSE`](LICENSE) for details.
