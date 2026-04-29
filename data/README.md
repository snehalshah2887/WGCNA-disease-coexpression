# Input Data

The processed input file (`gene_expression_filtered.csv`) is included in this
repository (11 MB). It is a filtered, annotated expression matrix derived from
the publicly available BrainSpan developmental transcriptomics atlas. The gene
list used for disease annotations was published in the supplemental material of
the associated paper (https://doi.org/10.3390/genes14111998).

This section describes the file format and how to reproduce it from raw
BrainSpan data if needed.

---

## Required input file

### `data/input/gene_expression_filtered.csv`

A combined gene annotation + expression matrix with the following structure:

| Columns | Content |
|---------|---------|
| 1–10 | Gene annotation (see column names below) |
| 11–534 | BrainSpan RNA-seq expression values (524 samples) |

**Annotation columns (1–10):**

| Column | Description |
|--------|-------------|
| `Ensembl` | Ensembl gene ID (e.g. `ENSG00000100393`) |
| `Gene_id` | Gene name / symbol |
| `Gene` | Alternative gene name |
| `Entrez` | Entrez gene ID |
| `ASD` | 1 if gene is associated with Autism Spectrum Disorder, else 0 |
| `ID` | 1 if gene is associated with Intellectual Disability, else 0 |
| `Seizures` | 1 if gene is associated with Seizures, else 0 |
| `Hypotonia` | 1 if gene is associated with Hypotonia, else 0 |
| `LangImp` | 1 if gene is associated with Language Impairment, else 0 |
| `X22q13` | 1 if gene is located in the 22q13 region (PMS gene), else 0 |

**Expression values:** RNA-seq RPKM from the BrainSpan developmental
transcriptomics atlas. Pre-filtered to genes with mean expression ≥ 0.3
across all 524 samples.

---

## Data sources

### BrainSpan expression data

Raw RNA-seq data are available from the BrainSpan Atlas of the Developing
Human Brain:

> https://www.brainspan.org/static/download.html

Download: **RNA-Seq Gencode v10 summarized to genes** (RPKM values).

Reference: Miller JA, et al. (2014). Transcriptional landscape of the
prenatal human brain. *Nature* **508**: 199–206.
https://doi.org/10.1038/nature13185

### Disease gene annotations

Disease gene lists (ASD, ID, Seizures, Hypotonia, Language Impairment) were
compiled from:

- **SFARI Gene database** (ASD): https://gene.sfari.org
- **SysID database** (ID): https://sysid.cmbi.umcn.nl
- **OMIM** (Seizures, Hypotonia, Language Impairment): https://omim.org
- **ClinVar**: https://www.ncbi.nlm.nih.gov/clinvar/

### 22q13 / PMS gene list

Genes in the 22q13.3 chromosomal region were identified using Ensembl
(GRCh38 coordinates: chr22:45,700,000–51,304,566).

---

## Data preparation

To reproduce the input file from raw BrainSpan downloads:

1. Download RNA-seq RPKM matrix from BrainSpan (link above)
2. Map probes/rows to Ensembl IDs using the provided gene metadata file
3. Filter to genes with mean RPKM ≥ 0.3 across all samples
4. Add disease annotation columns (ASD, ID, Seizures, Hypotonia, LangImp,
   X22q13) by joining to curated gene lists
5. Save the combined matrix as `data/input/gene_expression_filtered.csv`

The final file should have **2,116 rows** (genes) and **534 columns**
(10 annotation + 524 expression).
