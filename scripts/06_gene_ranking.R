# =============================================================================
# Module 6: 22q13 Gene Ranking by Co-expression Connectivity
# 22q13 / Phelan-McDermid Syndrome WGCNA Pipeline
#
# Description: Ranks 22q13 (PMS) genes by their average weighted topological
#              overlap (TOM) connectivity to genes associated with each of five
#              disease phenotypes: ASD, Intellectual Disability (ID), Seizures,
#              Hypotonia, and Language Impairment.
#
#              For each disease the FULL network TOM is subset to rows
#              corresponding to that disease's annotated genes. The top-25%
#              (q75) strongest TOM values are retained; column means then give
#              each gene's connectivity score to that disease gene set. PMS
#              genes are ranked by this score and written to individual CSVs.
#              A combined wide-format summary across all diseases is also
#              written.
#
#              Module-eigengene connectivity (kME) is computed once and joined
#              to every disease ranking as a supplementary intramodular metric.
#
#              A composite score is computed per gene per disease:
#                composite = 0.70 * rank_norm(mean_TOM) + 0.30 * rank_norm(kME)
#              where rank_norm maps values to [0,1] via rank / n_genes.
#              Genes with kME_own_module < 0.5 or CV (sd/mean) >= 1 are flagged
#              but retained so the researcher can apply filters as needed.
#
# Usage:
#   Rscript 06_gene_ranking.R <network_objects.rds> \
#                              <tom_objects.rds>     \
#                              <output_dir>
#
# Arguments:
#   network_objects.rds : Output of Module 2 (datExpr, MEs, final1, ...)
#   tom_objects.rds     : Output of Module 3 (full TOM, final disease genes)
#   output_dir          : Directory to write all outputs
#
# Outputs (written to output_dir/):
#   ranked_ASD_PMS.csv           : PMS genes ranked by composite score vs ASD
#   ranked_ID_PMS.csv            : PMS genes ranked by composite score vs ID
#   ranked_Seizures_PMS.csv      : PMS genes ranked by composite score vs Seizures
#   ranked_Hypotonia_PMS.csv     : PMS genes ranked by composite score vs Hypotonia
#   ranked_LangImp_PMS.csv       : PMS genes ranked by composite score vs Lang Imp
#   ranked_all_diseases_PMS.csv  : Wide-format summary — all five disease scores
#                                   per PMS gene in one table
#   top5_PMS_disease_genes.csv   : Top-5 PMS genes per disease by composite score
#                                   (hub-filtered: kME >= 0.5, CV < 1)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
  library(matrixStats)
})

options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 06_gene_ranking.R <network_objects.rds> <tom_objects.rds> <output_dir>",
       call. = FALSE)
}
network_objects_file <- args[1]
tom_objects_file     <- args[2]
output_dir           <- args[3]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("=== Module 6: Gene Ranking ===")

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
net_obj  <- readRDS(network_objects_file)
tom_obj  <- readRDS(tom_objects_file)

datExpr  <- net_obj$datExpr   # samples x genes; colnames = Ensembl IDs
MEs      <- net_obj$MEs       # module eigengenes
final1   <- net_obj$final1    # ALL datExpr genes with module assignments
                               # columns: Ensembl, Gene_symbol, ASD, ID,
                               #          Seizures, Hypotonia, Lang_Imp,
                               #          PMS, mlabel, mcolor

# Full N x N TOM — Ensembl IDs as both row and column names.
# This is distinct from modTom (which is only the last colour module's
# sub-matrix and must NOT be used here).
full_TOM <- tom_obj$TOM

# TOMsimilarity() does not copy dimnames from the adjacency matrix, so the
# saved TOM has no row/column names. Restore them from datExpr, which holds
# genes as columns in exactly the same order as the TOM rows/columns.
if (is.null(rownames(full_TOM))) {
  gene_ids           <- colnames(datExpr)
  rownames(full_TOM) <- gene_ids
  colnames(full_TOM) <- gene_ids
  message("TOM dimnames restored from datExpr (", length(gene_ids), " genes)")
}

# Disease-annotated gene subset (filtered to genes present in datExpr).
# Column names use I_D (not ID) for intellectual disability.
# Columns: Ensembl, Gene_symbol, ASD, I_D, Seizures, Hypotonia, Lang_Imp,
#           PMS, mlabel, mcolor, total
final    <- tom_obj$final

# ---------------------------------------------------------------------------
# Disease map: output label -> column name in 'final'
# ---------------------------------------------------------------------------
diseases <- c(
  ASD       = "ASD",
  ID        = "I_D",
  Seizures  = "Seizures",
  Hypotonia = "Hypotonia",
  LangImp   = "Lang_Imp"
)

# Genes present in the TOM
tom_genes <- rownames(full_TOM)

# PMS gene Ensembl IDs (must be present in the TOM)
pms_ensembl <- intersect(final$Ensembl[final$PMS == 1], tom_genes)
message("PMS genes in TOM : ", length(pms_ensembl))

if (length(pms_ensembl) == 0) {
  stop("No PMS genes found in TOM. Verify that final$PMS and TOM row names use the same Ensembl ID format.")
}

# ---------------------------------------------------------------------------
# Module-eigengene connectivity (kME) — computed once for all genes.
# signedKME() returns a data frame: rows = Ensembl IDs, cols = ME<color>.
# For each PMS gene we extract the kME to its OWN module eigengene.
# ---------------------------------------------------------------------------
message("Computing signed kME for all genes ...")
kme_all <- signedKME(datExpr, MEs, outputColumnName = "ME", corFnc = "cor")

# Gene-level info for PMS genes (module colour, gene symbol, kME)
pms_info <- final1[final1$Ensembl %in% pms_ensembl,
                   c("Ensembl", "Gene_symbol", "mcolor", "mlabel")]

pms_info$kME_own_module <- mapply(
  function(ens, label) {
    col_name <- paste0("ME", label)   # MEs are named ME1, ME2, ... (numeric labels)
    if (ens %in% rownames(kme_all) && col_name %in% colnames(kme_all)) {
      kme_all[ens, col_name]
    } else {
      NA_real_
    }
  },
  pms_info$Ensembl, pms_info$mlabel
)

# ---------------------------------------------------------------------------
# Helper: rank PMS genes by TOM connectivity to a given disease gene set
#
#   full_TOM    : full N x N TOM matrix (Ensembl row/col names)
#   final       : disease-annotated gene data frame (from tom_obj)
#   pms_ensembl : character vector of PMS gene Ensembl IDs in TOM
#   disease_col : column name in 'final' for this disease (e.g. "ASD")
#   label       : short label used in messages and the 'disease' column
#
# Returns a data frame or NULL if no disease genes or no PMS genes found.
# ---------------------------------------------------------------------------
rank_pms_vs_disease <- function(full_TOM, final, pms_ensembl,
                                disease_col, label) {

  dis_ensembl <- intersect(final$Ensembl[final[[disease_col]] == 1],
                           rownames(full_TOM))
  message("  ", label, ": ", length(dis_ensembl), " disease genes in TOM")

  if (length(dis_ensembl) == 0) {
    message("  No ", label, " genes found in TOM — skipping.")
    return(NULL)
  }

  # Subset TOM: rows = disease genes, columns = all genes in network
  TOM_dis <- full_TOM[dis_ensembl, , drop = FALSE]

  # Retain only the top 25% strongest TOM connections (q75 threshold).
  # Threshold is computed from this disease-specific sub-matrix so the
  # cutoff reflects the distribution of connections to THESE disease genes.
  q75 <- quantile(TOM_dis, na.rm = TRUE)[4][[1]]
  TOM_thresh <- TOM_dis
  TOM_thresh[TOM_thresh < q75] <- NA

  # Column means and SDs: each gene's average connectivity to this disease set
  aves <- colMeans(TOM_thresh, na.rm = TRUE)
  sds  <- colSds(TOM_thresh,   na.rm = TRUE)
  names(sds) <- names(aves)   # colSds drops names; restore them

  pms_in <- intersect(pms_ensembl, names(aves))
  if (length(pms_in) == 0) {
    message("  PMS Ensembl IDs not found in TOM columns — skipping.")
    return(NULL)
  }

  data.frame(
    Ensembl         = pms_in,
    mean_TOM        = as.numeric(aves[pms_in]),
    sd_TOM          = as.numeric(sds[pms_in]),
    n_disease_genes = length(dis_ensembl),
    q75_threshold   = as.numeric(q75),
    disease         = label,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Run ranking for every disease
# ---------------------------------------------------------------------------
all_ranked <- list()

for (label in names(diseases)) {
  col <- diseases[[label]]
  message("Ranking PMS genes vs ", label, " ...")

  res <- rank_pms_vs_disease(full_TOM, final, pms_ensembl, col, label)
  if (is.null(res)) next

  # Join gene symbol, module colour, kME
  res <- merge(res, pms_info, by = "Ensembl", all.x = TRUE)

  # Composite score: 0.70 * rank_norm(mean_TOM) + 0.30 * rank_norm(kME)
  # rank_norm maps each metric to [0,1] so the two scales are comparable.
  # NAs in kME are treated as rank 0 (conservative — avoids promoting genes
  # whose module membership could not be confirmed).
  n <- nrow(res)
  tom_rank_norm <- rank(res$mean_TOM,        ties.method = "average",  na.last = "keep") / n
  kme_rank_norm <- rank(res$kME_own_module,  ties.method = "average",  na.last = "keep") / n
  kme_rank_norm[is.na(kme_rank_norm)] <- 0

  res$composite_score <- round(0.70 * tom_rank_norm + 0.30 * kme_rank_norm, 6)

  # Quality flags (retained in output; not used to drop genes)
  res$hub_gene  <- !is.na(res$kME_own_module) & res$kME_own_module >= 0.5
  res$stable_CV <- !is.na(res$sd_TOM) & !is.na(res$mean_TOM) &
                   (res$sd_TOM / res$mean_TOM) < 1

  # Sort by composite score (descending) and assign rank
  res <- res[order(res$composite_score, decreasing = TRUE), ]
  res$rank <- seq_len(nrow(res))

  res <- res[, c("rank", "Ensembl", "Gene_symbol", "mcolor", "mlabel",
                 "composite_score", "mean_TOM", "sd_TOM", "kME_own_module",
                 "hub_gene", "stable_CV",
                 "n_disease_genes", "q75_threshold", "disease")]

  out_file <- file.path(output_dir, paste0("ranked_", label, "_PMS.csv"))
  write.csv(res, out_file, row.names = FALSE)
  message("  Written: ", basename(out_file), " (", nrow(res), " PMS genes)")

  all_ranked[[label]] <- res
}

# ---------------------------------------------------------------------------
# Wide-format combined summary: one row per PMS gene, one score per disease
# ---------------------------------------------------------------------------
if (length(all_ranked) > 0) {
  # Base columns from the first available result (gene identity + kME)
  first <- all_ranked[[1]]
  wide  <- first[, c("Ensembl", "Gene_symbol", "mcolor", "kME_own_module",
                     "hub_gene", "stable_CV")]

  for (label in names(all_ranked)) {
    sub <- all_ranked[[label]][, c("Ensembl", "composite_score", "mean_TOM", "rank")]
    colnames(sub)[2:4] <- c(paste0("composite_",  label),
                             paste0("mean_TOM_",   label),
                             paste0("rank_",        label))
    wide <- merge(wide, sub, by = "Ensembl", all.x = TRUE)
  }

  # Sort by the first available composite score (descending)
  first_score_col <- paste0("composite_", names(all_ranked)[1])
  wide <- wide[order(wide[[first_score_col]], decreasing = TRUE,
                     na.last = TRUE), ]

  write.csv(wide,
            file.path(output_dir, "ranked_all_diseases_PMS.csv"),
            row.names = FALSE)
  message("Written: ranked_all_diseases_PMS.csv (",
          nrow(wide), " PMS genes x ", length(all_ranked), " diseases)")

  # -------------------------------------------------------------------------
  # Top-5 per disease: hub-filtered (kME >= 0.5, CV < 1), ranked by composite
  # -------------------------------------------------------------------------
  top5_list <- lapply(names(all_ranked), function(label) {
    df <- all_ranked[[label]]

    # Apply publication-grade filters; warn but don't crash if no genes pass
    filtered <- df[df$hub_gene & df$stable_CV, ]
    if (nrow(filtered) == 0) {
      message("  WARNING: No genes pass hub+CV filter for ", label,
              " — using unfiltered top 5")
      filtered <- df
    }

    # Top 5 by composite score (already sorted descending)
    top5 <- head(filtered, 5)
    top5$top5_rank <- seq_len(nrow(top5))
    top5[, c("top5_rank", "Ensembl", "Gene_symbol", "mcolor", "mlabel",
             "composite_score", "mean_TOM", "sd_TOM", "kME_own_module",
             "hub_gene", "stable_CV", "n_disease_genes", "q75_threshold",
             "disease")]
  })

  top5_all <- do.call(rbind, top5_list)
  rownames(top5_all) <- NULL

  write.csv(top5_all,
            file.path(output_dir, "top5_PMS_disease_genes.csv"),
            row.names = FALSE)
  message("Written: top5_PMS_disease_genes.csv (",
          nrow(top5_all), " rows — up to 5 per disease)")
} else {
  message("No disease rankings produced — check input data and disease flags.")
}

message("=== Module 6 complete ===")
