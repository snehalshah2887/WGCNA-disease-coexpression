# =============================================================================
# Module 5: GO Term Enrichment Analysis & Visualisation
# 22q13 / Phelan-McDermid Syndrome WGCNA Pipeline
#
# Description: Runs GO Biological Process enrichment analysis using
#              hypergeometric testing (phyper) with annotations from
#              org.Hs.eg.db and GO.db. Background universe is restricted to
#              the 2,116 genes present in the co-expression network, providing
#              a more stringent and biologically appropriate enrichment test
#              than whole-genome backgrounds.
#
# Usage:
#   Rscript 05_go_analysis.R <network_objects.rds> <output_dir>
#
# Outputs (written to output_dir/):
#   filtered_go/        : enrichment results per colour module (CSV)
#   Goterms_22q13.csv   : GO terms containing at least one PMS/22q13 gene
#   top5_go_all.csv     : Top 5 GO terms per module (all modules combined)
#   plots/              : Bar plots per module + combined bar plot
# =============================================================================

local({
  pkgs <- c("org.Hs.eg.db", "GO.db", "AnnotationDbi")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
})

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 05_go_analysis.R <network_objects.rds> <output_dir>",
       call. = FALSE)
}
network_objects_file <- args[1]
output_dir           <- args[2]

filtered_go_dir <- file.path(output_dir, "filtered_go")
plots_dir       <- file.path(output_dir, "plots")
for (d in c(output_dir, filtered_go_dir, plots_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("=== Module 5: GO Term Enrichment Analysis ===")

# ---------------------------------------------------------------------------
# Load network objects
# ---------------------------------------------------------------------------
obj     <- readRDS(network_objects_file)
final1  <- obj$final1
datExpr <- obj$datExpr

universe_genes <- colnames(datExpr)   # Ensembl IDs of all network genes
color_list     <- sort(unique(final1$mcolor[final1$mcolor != "grey"]))
message("Modules to analyse: ", paste(color_list, collapse = ", "))

# ---------------------------------------------------------------------------
# Build gene -> GO:BP lookup table (universe genes only)
# ---------------------------------------------------------------------------
message("Fetching GO:BP annotations from org.Hs.eg.db …")
go_ann <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = universe_genes,
  columns = c("ENSEMBL", "GO", "ONTOLOGY"),
  keytype = "ENSEMBL"
) %>%
  filter(ONTOLOGY == "BP", !is.na(GO)) %>%
  dplyr::select(ENSEMBL, GO) %>%
  distinct()

# GO term descriptions
go_desc <- AnnotationDbi::select(
  GO.db,
  keys    = unique(go_ann$GO),
  columns = c("GOID", "TERM"),
  keytype = "GOID"
) %>%
  dplyr::rename(GO = GOID, Description = TERM)

# For each GO term: which universe genes are annotated?
go_to_genes <- split(go_ann$ENSEMBL, go_ann$GO)   # list: GO_id -> [Ensembl...]

N <- length(universe_genes)   # total background size

# ---------------------------------------------------------------------------
# Helper: hypergeometric enrichment for one module
# ---------------------------------------------------------------------------
enrich_module <- function(module_genes, go_to_genes, universe_genes, go_desc,
                           p_cut = 0.05, q_cut = 0.20,
                           min_gs = 10, max_gs = 500) {
  n <- length(module_genes)
  N <- length(universe_genes)

  # Filter GO terms by size in universe
  go_sizes <- sapply(go_to_genes, length)
  go_to_genes <- go_to_genes[go_sizes >= min_gs & go_sizes <= max_gs]

  if (length(go_to_genes) == 0) return(NULL)

  rows <- lapply(names(go_to_genes), function(go_id) {
    go_genes <- go_to_genes[[go_id]]
    K  <- length(go_genes)
    k  <- sum(module_genes %in% go_genes)
    if (k == 0) return(NULL)
    pval <- phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
    hit_genes <- intersect(module_genes, go_genes)
    data.frame(
      GOID      = go_id,
      GeneRatio = paste0(k, "/", n),
      BgRatio   = paste0(K, "/", N),
      Count     = k,
      pvalue    = pval,
      geneID    = paste(hit_genes, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })

  res <- bind_rows(rows)
  if (is.null(res) || nrow(res) == 0) return(NULL)

  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res$qvalue   <- p.adjust(res$pvalue, method = "BH")   # same as p.adjust here

  res <- res %>%
    filter(p.adjust <= q_cut) %>%
    left_join(go_desc, by = c("GOID" = "GO")) %>%
    arrange(pvalue)

  if (nrow(res) == 0) return(NULL)
  res
}

# ---------------------------------------------------------------------------
# Run enrichment for every colour module
# ---------------------------------------------------------------------------
go_results <- list()
top5_list  <- list()

for (color in color_list) {
  module_genes <- final1$Ensembl[final1$mcolor == color]
  message("  ", color, " module — ", length(module_genes), " genes")

  res <- tryCatch(
    enrich_module(module_genes, go_to_genes, universe_genes, go_desc),
    error = function(e) {
      message("    Skipped (", conditionMessage(e), ")")
      NULL
    }
  )

  go_results[[color]] <- res

  if (is.null(res) || nrow(res) == 0) {
    message("    No significant GO terms for ", color)
    next
  }

  write.csv(res,
            file.path(filtered_go_dir, paste0(color, "_filtered_module.csv")),
            row.names = FALSE)

  res$module          <- color
  res$log_Adj_Pvalue  <- -log10(res$p.adjust)
  top5_list[[color]]  <- head(res, 5)

  # Bar plot — top 10 terms
  plot_df <- head(res, 10)
  plot_df$Description <- factor(
    plot_df$Description,
    levels = rev(plot_df$Description)
  )

  p <- ggplot(plot_df, aes(x = Description, y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = color, color = color) +
    coord_flip() +
    xlab("") + ylab("-log10(adjusted p-value)") +
    ggtitle(paste(color, "module — GO Biological Process")) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line        = element_line(color = "black"))

  ggsave(file.path(plots_dir, paste0("go_barplot_", color, ".png")),
         plot = p, width = 10, height = 6)
  message("    Saved bar plot for ", color)
}

# ---------------------------------------------------------------------------
# Goterms_22q13.csv — GO terms with at least one PMS/22q13 gene
# ---------------------------------------------------------------------------
pms_ensembl <- final1$Ensembl[final1$PMS == 1]

goterms_22q13 <- bind_rows(lapply(names(go_results), function(color) {
  res <- go_results[[color]]
  if (is.null(res) || nrow(res) == 0) return(NULL)
  has_pms <- sapply(res$geneID, function(ids) {
    any(strsplit(ids, "/")[[1]] %in% pms_ensembl)
  })
  if (!any(has_pms)) return(NULL)
  sub <- res[has_pms, c("Description", "geneID", "p.adjust")]
  colnames(sub) <- c("Term", "Genes", "FDR")
  sub$module <- color
  sub
}))

if (is.null(goterms_22q13) || nrow(goterms_22q13) == 0) {
  goterms_22q13 <- data.frame(Term = character(), Genes = character(),
                               FDR = numeric(), module = character())
}
goterms_22q13 <- distinct(goterms_22q13[order(goterms_22q13$Term), ])
write.csv(goterms_22q13, file.path(output_dir, "Goterms_22q13.csv"), row.names = FALSE)
message("Written: Goterms_22q13.csv (", nrow(goterms_22q13), " terms)")

# ---------------------------------------------------------------------------
# top5_go_all.csv + combined bar plot
# ---------------------------------------------------------------------------
if (length(top5_list) > 0) {
  go_df <- bind_rows(top5_list)
  write.csv(go_df, file.path(output_dir, "top5_go_all.csv"), row.names = FALSE)
  message("Written: top5_go_all.csv")

  go_df_plot <- go_df[!is.na(go_df$log_Adj_Pvalue), ]
  if (nrow(go_df_plot) > 0) {
    go_df_plot$Description <- make.unique(as.character(go_df_plot$Description))
    go_df_plot$Description <- fct_reorder(go_df_plot$Description, go_df_plot$log_Adj_Pvalue)

    p_all <- ggplot(go_df_plot,
                    aes(x = Description, y = log_Adj_Pvalue, fill = module)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      xlab("") + ylab("-log10(adjusted p-value)") +
      ggtitle("Top 5 GO BP terms per module") +
      geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line        = element_line(color = "black"),
            axis.text.y      = element_text(size = 7))

    ggsave(file.path(plots_dir, "go_barplot_all_modules.png"),
           plot = p_all, width = 14, height = max(8, nrow(go_df_plot) * 0.3))
    message("Written: go_barplot_all_modules.png")
  }
} else {
  message("No significant GO terms in any module")
  write.csv(data.frame(), file.path(output_dir, "top5_go_all.csv"), row.names = FALSE)
}

message("=== Module 5 complete ===")
