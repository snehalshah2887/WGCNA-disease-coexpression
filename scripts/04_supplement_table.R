# =============================================================================
# Module 4: Supplemental Co-expression Table (22q13 Genes)
# 22q13 / Phelan-McDermid Syndrome WGCNA Pipeline
#
# Description: Reads the consolidated Cytoscape edge files (CytoscapeInput-
#              edges-10_29_22.txt), identifies all co-expressed partners for
#              each 22q13/PMS gene per module, and builds a wide-format
#              supplemental table.
#
# Usage:
#   Rscript 04_supplement_table.R <tom_objects.rds> \
#                                  <network_objects.rds> \
#                                  <cytoscape_dir> \
#                                  <output_dir>
#
# Arguments:
#   tom_objects.rds     : Output of Module 3
#   network_objects.rds : Output of Module 2
#   cytoscape_dir       : Directory containing CytoscapeInput-edges-10_29_22.txt
#   output_dir          : Directory to write outputs
#
# Outputs (written to output_dir/):
#   Table22q13.csv  : Co-expression table for 22q13 genes across modules
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 04_supplement_table.R <tom_objects.rds> <network_objects.rds> <cytoscape_dir> <output_dir>",
       call. = FALSE)
}
tom_file             <- args[1]
network_objects_file <- args[2]
cytoscape_dir        <- args[3]
output_dir           <- args[4]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("=== Module 4: Supplemental Co-expression Table ===")

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
tom_obj    <- readRDS(tom_file)
net_obj    <- readRDS(network_objects_file)
final      <- tom_obj$final
color_list <- tom_obj$color_list

myvars  = c('Ensembl', 'Gene_id', 'Gene', 'Entrez', 'ASD', 'ID', 'Seizures', 'Hypotonia', 'LangImp', 'X22q13')
gene_final      <- net_obj$genedata
gene_annotation <- gene_final[myvars]

# Function to create supplement table
create_supplement_table <- function(col_mod)
{
  suppTableFlag = 0

  edges = read.delim2(file = file.path(cytoscape_dir, "CytoscapeInput-edges-10_29_22.txt"))
  edges = edges[,c(5,6,3,4,1,2)]
  edges = rename(edges, fromNode = fromAltName, toNode = toAltName, fromAltName = fromNode, toAltName = toNode)
  pms = gene_annotation[gene_annotation$X22q13 == 1,]

  pms_from = edges[edges$fromNode %in% pms$Gene ,]
  if(nrow(pms_from)>0)
  {
    pms_from = pms_from[,1:2]
    names(pms_from) = c("fromGene", "toGene")
    pms_from$value <- 1

    pms_from <- unique(pms_from)

    from_suppTable = pms_from %>%
      pivot_wider(names_from = toGene,
                  values_from = value)

    for(i in 2:ncol(from_suppTable))
    {
      from_suppTable[,i] <- ifelse(from_suppTable[,i] == 1, colnames(from_suppTable[,i]),from_suppTable[,i])
    }

    from_suppTable = unite(from_suppTable, 2:ncol(from_suppTable), col = "co_expression", sep = ", ", remove = TRUE, na.rm = TRUE)
    names(from_suppTable) = c("X22q13","co_expression")
  }

  pms_to = edges[edges$toNode %in% pms$Gene ,]
  if(nrow(pms_to) > 0)
  {
    pms_to = pms_to[,1:2]
    names(pms_to) = c("fromGene", "toGene")
    pms_to$value <- 1

    pms_to <- unique(pms_to)

    to_suppTable = pms_to %>%
      pivot_wider(names_from = fromGene,
                  values_from = value)

    for(i in 2:ncol(to_suppTable))
    {
      to_suppTable[,i] <- ifelse(to_suppTable[,i] == 1, colnames(to_suppTable[,i]),to_suppTable[,i])
    }

    to_suppTable = unite(to_suppTable, 2:ncol(to_suppTable), col = "co_expression", sep = ", ", remove = TRUE, na.rm = TRUE)
    names(to_suppTable) = c("X22q13","co_expression")
  }


  if(nrow(pms_from) > 0 & nrow(pms_to) > 0)
  {
    suppTable = as.data.frame(rbind(from_suppTable, to_suppTable))
    suppTable$value <- 1
    suppTableFlag = 1
  }
  if(nrow(pms_from) > 0 & nrow(pms_to) == 0)
  {
    suppTable = as.data.frame(from_suppTable)
    suppTable$value <- 1
    suppTableFlag = 1
  }
  if(nrow(pms_from) == 0 & nrow(pms_to) > 0)
  {
    suppTable = as.data.frame(to_suppTable)
    suppTable$value <- 1
    suppTableFlag = 1
  }

  if(suppTableFlag == 1)
  {
    suppTable = suppTable %>%
      pivot_wider(names_from = "co_expression",
                  values_from = value)

    for(i in 2:ncol(suppTable))
    {
      suppTable[,i] <- ifelse(suppTable[,i] == 1, colnames(suppTable[,i]),suppTable[,i])
    }

    suppTable = unite(suppTable, 2:ncol(suppTable), col = "co_expression", sep = ", ", remove = TRUE, na.rm = TRUE)
    names(suppTable) = c("X22q13","co_expression")
    suppTable$module = col_mod

    return(suppTable)
  }
  else
  {
    suppTable = data.frame(X22q13 = character(),
               co_expression = character(),
               module = character(),
               stringsAsFactors=FALSE)

    return(suppTable)
  }

}


table_22q13 <- data.frame(X22q13 = character(),
                          co_expression = character(),
                          module = character(),
                          stringsAsFactors=FALSE)

for(j in color_list)
{
  suppT = create_supplement_table(j)
  if(nrow(suppT) > 0)
  {
    table_22q13 <- rbind(table_22q13, suppT)
  }
}
table_22q13$degree = lengths(gregexpr(",", table_22q13$co_expression)) + 1

write.csv(table_22q13, file.path(output_dir, "Table22q13.csv"))
message("Written: Table22q13.csv (", nrow(table_22q13), " rows)")
message("=== Module 4 complete ===")
