# =============================================================================
# Module 1: Data Preprocessing
# 22q13 / Phelan-McDermid Syndrome WGCNA Pipeline
#
# Description: Reads the combined gene expression file (gene annotations in
#              columns 1-10, BrainSpan expression in columns 11-534).
#              Extracts the gene annotation block and passes the file through
#              to the results directory for downstream modules.
#
#              NOTE: The input file is already merged and filtered (mean >= 0.3
#              threshold applied prior to this pipeline). No merge or filter
#              step is performed here.
#
# Usage:
#   Rscript 01_preprocess.R <gene_expression_filtered.csv> <output_dir>
#
# Arguments:
#   gene_expression_filtered.csv : Combined file — columns 1-10 are gene
#                                  annotation, columns 11-534 are BrainSpan
#                                  expression values (524 samples).
#   output_dir                   : Directory to write outputs
#
# Outputs (written to output_dir):
#   gene_expression_filtered.csv  : Copy of input (for downstream modules)
#   gene_annotation.rds           : Gene annotation object (columns 1-10)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
})

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_preprocess.R <gene_expression_filtered.csv> <output_dir>",
       call. = FALSE)
}
expr_file  <- args[1]
output_dir <- args[2]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("=== Module 1: Preprocessing ===")
message("Input file     : ", expr_file)
message("Output directory: ", output_dir)

# ---------------------------------------------------------------------------
# Read the combined gene expression file
# ---------------------------------------------------------------------------
gene_final <- read.csv(expr_file)

# Extract gene annotation block (first 10 columns)
gene <- gene_final[, 1:10]
names(gene) <- c("Ensembl","Gene_id","Gene","Entrez","ASD","ID","Seizures","Hypotonia","LangImp","X22q13")

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
write.csv(gene_final,
          file.path(output_dir, "gene_expression_filtered.csv"),
          row.names = FALSE)

saveRDS(gene, file.path(output_dir, "gene_annotation.rds"))

message("Written: gene_expression_filtered.csv")
message("Written: gene_annotation.rds")
message("Genes in file: ", nrow(gene_final))
message("=== Module 1 complete ===")
