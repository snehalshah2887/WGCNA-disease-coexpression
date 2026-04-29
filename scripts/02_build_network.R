# =============================================================================
# Module 2: WGCNA Network Construction, Module Annotation & Enrichment
# 22q13 / Phelan-McDermid Syndrome WGCNA Pipeline
#
# Description: Reads the filtered expression matrix, performs soft-threshold
#              selection, builds the co-expression network via blockwiseModules,
#              annotates modules with gene/disease metadata, runs Fisher exact
#              test enrichment analysis, and exports per-module gene tables.
#
# Usage:
#   Rscript 02_build_network.R <gene_expression_filtered.csv> \
#                               <gene_annotation.rds> \
#                               <output_dir>
#
# Arguments:
#   gene_expression_filtered.csv : Output of Module 1
#   gene_annotation.rds          : Output of Module 1
#   output_dir                   : Directory to write outputs
#
# Outputs (written to output_dir/):
#   color_modules/               : Per-color module gene files (V1 and V2)
#   colormodule.csv              : Full annotated module assignment table
#   summary_colors.csv           : Disease-count summary per module
#   plots/                       : Soft-threshold and dendrogram plots (PNG)
#   network_objects.rds          : Serialised R objects for downstream modules
#
# NOTE: blockwiseModules will also save TOM block files (datExpr-block.*.RData)
#       in output_dir/TOM/ — these are large files (~100s MB).
# =============================================================================

local({
  bioc_pkgs <- c("WGCNA", "biomaRt", "biomartr")
  missing   <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing missing Bioconductor packages: ", paste(missing, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
})

suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
  library(openxlsx)
  library(readtext)
  library(tidyr)
  library(ggplot2)
  library(gplots)
  library(RColorBrewer)
  library(VennDiagram)
  library(venn)
  library(stringr)
  library(biomartr)
  library(forcats)
  library(matrixStats)
})

options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 02_build_network.R <gene_expression_filtered.csv> <gene_annotation.rds> <output_dir>",
       call. = FALSE)
}
filtered_expr_file   <- args[1]
gene_annotation_file <- args[2]
output_dir           <- args[3]

# Create output sub-directories
plots_dir   <- file.path(output_dir, "plots")
modules_dir <- file.path(output_dir, "color_modules")
tom_dir     <- file.path(output_dir, "TOM")
for (d in c(output_dir, plots_dir, modules_dir, tom_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("=== Module 2: Network Construction + Annotation + Enrichment ===")

# ---------------------------------------------------------------------------
# Load inputs from Module 1
# ---------------------------------------------------------------------------
genedata = read.csv(filtered_expr_file)
gene     = readRDS(gene_annotation_file)

#Quick look of the data set
dim(genedata)
# Remove all annotation columns (1-10); keep expression values only (cols 11-534)
datExpr0 = genedata[,-c(1:10)]

################################# Pivoting the data ##########################################################

# Assigning the data to datExpr1(making a new variable for pivoting the data)
datExpr1 = datExpr0
# Adding a column row_num to datExpr1 and genedata
datExpr1$row_num = rownames(datExpr1)
genedata$row_num = rownames(genedata)
# Filtering the data to two column namely row_num and Ensembl to a new variable a
a <- genedata[,c("row_num","Ensembl")]
# Merging two dataframes a and datExpr1 for getting the ensembl id to the datExp variable for pivoting
b <- merge(a, datExpr1, by = "row_num")
# Dropping the row_num column
b <- b[,-1]
# Making the Ensembl column to naming the rows
rownames(b) <- b[,1]
# Dropping the Ensembl column for analysis
b <- b[,-1]
#Transposing the dataframe to desired format where gene instance are columns and samples are rows
datExpr = as.data.frame(t(b))
datExpr0 = datExpr
# Checking the matrix of the dataframe
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

# Checking the data for excessive missing values and identification of outlier samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# Removing offending genes and samples from the data
# if (!gsg$allOK)
# {
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0)
#     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
#   if (sum(!gsg$goodSamples)>0)
#     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
#   # Remove the offending genes and samples from the data:
#   datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
# }

# Checking for outliers in the data
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree
png(file.path(plots_dir, "sampleClustering.png"), width = 1600, height = 900)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#Plot a line to show the cut
abline(h = 40, col = "red")
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 40, minSize = 55)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

######################################### One step network and module analysis #############################################
######################################### Import gene names based on cluster developments ##################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
png(file.path(plots_dir, "soft_threshold.png"), width = 900, height = 1000)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = file.path(tom_dir, "datExpr"),
                       verbose = 3)
table(net$colors)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
png(file.path(plots_dir, "module_dendrogram.png"), width = 1200, height = 900)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
# Export module genes and eigengenes:
moduleGenes = data.frame(gene = names(datExpr), mlabel = moduleLabels, mcolor = moduleColors)
z = order(moduleGenes$mcolor)
x = moduleGenes[z,]
moduleGenes = x
moduleGenes = moduleGenes[,-c(1)]
colors = unique(moduleColors)

for (col in colors) {
  df = moduleGenes[moduleGenes$mcolor == col,]
  df$Ensembl = rownames(df)
  rownames(df) = seq(1:nrow(df))
  assign(paste0(col),df)
  write.csv(df,
            file.path(modules_dir, paste0(col, "_module_V1.csv")),
            row.names = FALSE)
}
moduleGenes$Ensembl = rownames(moduleGenes)
rownames(moduleGenes) = seq(1:nrow(moduleGenes))

# Extracting the columns like Ensembl, gene_symbol + asd, pms,id, seizure, hypotonia and lang_Imp
ensemble = gene
geneTree = net$dendrograms[[1]];

colnames(ensemble)
# Merging the dataset of ensembl and module gene for getting the modcolor and modlabel in one table
final = merge(ensemble,moduleGenes, by = "Ensembl")

final$Gene_id = NULL
final$Entrez = NULL
Sep = data.frame(Ensembl = c("ENSG00000100167"),
                 Gene = c("Sept-3"),
                 ASD = c(0),
                 ID = c(0),
                 Seizures = c(0),
                 Hypotonia = c(0),
                 LangImp = c(0),
                 X22q13 = c(1),
                 mlabel = c(1),
                 mcolor = c("turquoise"))
oth = data.frame(Ensembl = c("ENSG00000213683"),
                 Gene = c("AC002056.3"),
                 ASD = c(0),
                 ID = c(0),
                 Seizures = c(0),
                 Hypotonia = c(0),
                 LangImp = c(0),
                 X22q13 = c(1),
                 mlabel = c(9),
                 mcolor = c("magenta"))

final = rbind(oth,final)

# Renaming the final dataframe as per required column name
names(final) = c("Ensembl","Gene_symbol", "ASD", "ID", "Seizures", "Hypotonia", "Lang_Imp", "PMS","mlabel", "mcolor")
# Write the file in csv format
write.csv(final,
          file.path(output_dir, "colormodule.csv"),
          row.names = FALSE)

final = read.csv(file.path(output_dir, "colormodule.csv"))
test = final %>% dplyr::group_by(mcolor) %>% dplyr::summarise(ASD = sum(ASD), ID = sum(ID), Seizures = sum(Seizures), Hypotonia = sum(Hypotonia),
                                                              Lang_Imp = sum(Lang_Imp), PMS = sum(PMS))
write.csv(test,
          file.path(output_dir, "summary_colors.csv"))

# Gather function helps in combining the disease columns into one column. Here 3:8 means asd - pms columns
# and 2:7 means the 1st column is mcolor and followed by 6 disease which needs to be gathered.
c = gather(final[,c(10,3:8)], disease, count, 2:7)
# Summing the count of disease and color together
bar_plot_data = c %>% dplyr::group_by(mcolor, disease) %>% dplyr::summarise(count = sum(count))
bar_plot_data = filter(bar_plot_data, mcolor!= 'grey')
# Plotting the bar graph
p_bar <- ggplot(data=bar_plot_data, aes(x=mcolor, y=count, fill=disease)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c('blue','green','red','yellow','magenta','orange'))
ggsave(file.path(plots_dir, "disease_distribution_barplot.png"), plot = p_bar, width = 12, height = 6)

# Creating the heat map data for enrichment analysis
heat_map_data = as.data.frame(spread(bar_plot_data, disease, count))
# Making a table with mcolor and the count of genes in it
mod_table = as.data.frame(table(final$mcolor))
# Labelling the columns
names(mod_table) = c("mcolor",'mod_totals')
# Merging heat_map_tables and mcolor for getting the total count+disease wise distribution of genes
# in each color module
heat_map_data = merge(heat_map_data,mod_table, by = 'mcolor')
# Naming the color modules as rownames
rownames(heat_map_data) = heat_map_data[,1]
# removing the modcolor column and creating a table
heat_map_data = heat_map_data[,2:8]
# Assigning mat_data variable in matrix form the heat_map_data
mat_data = as.matrix(heat_map_data)

################ ENRICHMENT ANALYSIS USING FISHER EXACT TEST USING FDR FOR COLOR MODULES STARTS ####################

# The total number of genes
total_genes = nGenes
#Assigning the mod_sum variable the data of mod_data
mod_sum = mat_data
# Empty variable matrix formed for _p and _or for 6 disease columns
mat_p=matrix(data=NA, nrow=nrow(mod_sum), ncol=6)
mat_or=matrix(data=NA, nrow=nrow(mod_sum), ncol=6)

# for loop to calculate the calculation of the enrichment analysis basically the pvalue and OR results of the analysis
for (row in 1:nrow(mod_sum)){
  for (col in 1:6){
    mod_total<-as.numeric(mod_sum[row,7])
    mod_count<- as.numeric(mod_sum[row,col])
    mod_non<- mod_total-mod_count
    non_count<-sum(as.numeric(mod_sum[,col])) - mod_count
    non_non <- (total_genes-mod_total)-non_count

    contigency<-matrix(c(mod_count,mod_non,non_count,non_non),2,2)
    results<-fisher.test(contigency,)
    p_val<-results[[1]]
    OR<-results[[3]]

    mat_p[row,col]<-p_val
    mat_or[row,col]<-OR
  }
}

# making the colormodules and colnames of disease into the table
rownames(mat_p)<-rownames(mod_sum)
rownames(mat_or)<-rownames(mod_sum)
colnames(mat_p)<-colnames(mod_sum)[1:6]
colnames(mat_or)<-colnames(mod_sum)[1:6]
# Adjust the p value for the correctness
ad_p<-c()

### now correct P-values for each test using FDR method
for (i in 1:ncol(mat_p)){
  ad_p<-cbind(ad_p,p.adjust(mat_p[,i],"fdr"))
}
colnames(ad_p) <-colnames(mat_p)
# Log-10 transformation
mat_p<-ad_p
log_p<- log10(mat_p)*-1      ### -log10 transformation for figure
log_p <- round(log_p,2)
log_p<-log_p[!(row.names(log_p) %in% c("grey")),]

# creates a own color palette from red to white for enrichment analysis for table
my_palette <- colorRampPalette(c("white", "orange", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,1,length=100),  # for white
               seq(1.01,5.1,length=100),           # for orange
               seq(5.12,10.1,length=100)) # red

# making the heatmap for the analysis
png(file.path(plots_dir, "enrichment_heatmap.png"), width = 1500, height = 1500, res = 300, pointsize = 8)
heatmap.2(log_p,
          cellnote = log_p,  # same data set for cell labels
          main = "Gene Significance and module membership", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,6),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          #Colv="NA",
          key.xlab="-log10 p-value")            # turn off column clustering
# FLAG 7 (original): dev.off() was called without a matched png()/pdf() in the
# original script (those calls were commented out). Corrected here to match
# the png() opened above.
dev.off()

###################### ENRICHMENT ANALYSIS USING FISHER EXACT TEST ENDS ############################

###################### MERGED DATA OUTPUT GENE SYMBOL, COLOR MODULES , AND DISEASES #################
for (col in colors) {
  df = final[final$mcolor == col,]
  assign(paste0(col),df)
  write.csv(df,
            file.path(modules_dir, paste0(col, "_module_V2.csv")),
            row.names = FALSE)
}
###################### MERGED DATA OUTPUT GENE SYMBOL, COLOR MODULES , AND DISEASES ENDS ################################

# ---------------------------------------------------------------------------
# Construct final1: all datExpr genes with module assignments (needed by
# Module 6 gene ranking). Distinct from 'final' which holds disease genes only.
# ---------------------------------------------------------------------------
final1 = merge(genedata[, c("Ensembl","Gene","ASD","ID","Seizures","Hypotonia","LangImp","X22q13")],
               moduleGenes[, c("Ensembl","mlabel","mcolor")],
               by = "Ensembl")
names(final1) = c("Ensembl","Gene_symbol","ASD","ID","Seizures","Hypotonia","Lang_Imp","PMS","mlabel","mcolor")

# ---------------------------------------------------------------------------
# Save all objects needed by downstream modules
# ---------------------------------------------------------------------------
network_objects <- list(
  datExpr       = datExpr,
  datExpr1      = datExpr1,
  genedata      = genedata,
  moduleColors  = moduleColors,
  moduleLabels  = moduleLabels,
  mergedColors  = mergedColors,
  MEs           = MEs,
  net           = net,
  colors        = colors,
  nGenes        = nGenes,
  nSamples      = nSamples,
  final         = final,
  final1        = final1,
  gene          = gene,
  mat_data      = mat_data,
  log_p         = log_p
)
saveRDS(network_objects, file.path(output_dir, "network_objects.rds"))
message("Written: network_objects.rds")
message("=== Module 2 complete ===")
