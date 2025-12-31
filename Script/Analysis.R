# -------------------------------
# 1. Load required libraries
# -------------------------------
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(pheatmap)

# -------------------------------
# 2. Define file paths for WIG files
# -------------------------------
wig_files <- list(
  H3K4me3_TSA = "H3K4me3_TSA.wig.gz",
  H3K4me3_Vehicle = "H3K4me3_Vehicle.wig.gz",
  AcH3K9_14_TSA = "AcH3K9_14_TSA.wig.gz",
  AcH3K9_14_Vehicle = "AcH3K9_14_Vehicle.wig.gz",
  AcH4K12_TSA = "AcH4K12_TSA.wig.gz",
  AcH4K12_Vehicle = "AcH4K12_Vehicle.wig.gz",
  panAcH2B_TSA = "panAcH2B_TSA.wig.gz",
  panAcH2B_Vehicle = "panAcH2B_Vehicle.wig.gz"
)

# -------------------------------
# 3. Import all WIG files
# -------------------------------
wig_list <- lapply(wig_files, function(f) import(f, format="wig"))

# -------------------------------
# 4. Import mouse gene annotation and define promoters
# -------------------------------
gtf_file <- "Mus_musculus.GRCm39.109.gtf"  # update path
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")

# Define promoters (TSS +/- 1kb)
promoters_regions <- promoters(genes(txdb), upstream=1000, downstream=1000)
gene_names <- names(promoters_regions)

# -------------------------------
# 5. Compute mean signal per promoter per sample
# -------------------------------
signal_matrix <- matrix(0, nrow=length(promoters_regions), ncol=length(wig_list))
rownames(signal_matrix) <- gene_names
colnames(signal_matrix) <- names(wig_list)

for(i in seq_along(promoters_regions)){
  region <- promoters_regions[i]
  for(j in seq_along(wig_list)){
    wig <- wig_list[[j]]
    overlap <- findOverlaps(region, wig)
    hits <- subjectHits(overlap)
    if(length(hits) > 0){
      signal_matrix[i,j] <- mean(wig$score[hits])
    } else {
      signal_matrix[i,j] <- 0
    }
  }
}

# -------------------------------
# 6. Calculate Fold Change (TSA vs Vehicle)
# -------------------------------
fc_matrix <- matrix(0, nrow=length(promoters_regions), ncol=4)
rownames(fc_matrix) <- gene_names
colnames(fc_matrix) <- c("H3K4me3_FC", "AcH3K9_14_FC", "AcH4K12_FC", "panAcH2B_FC")

fc_matrix[,1] <- signal_matrix[, "H3K4me3_TSA"] / (signal_matrix[, "H3K4me3_Vehicle"] + 1e-6)
fc_matrix[,2] <- signal_matrix[, "AcH3K9_14_TSA"] / (signal_matrix[, "AcH3K9_14_Vehicle"] + 1e-6)
fc_matrix[,3] <- signal_matrix[, "AcH4K12_TSA"] / (signal_matrix[, "AcH4K12_Vehicle"] + 1e-6)
fc_matrix[,4] <- signal_matrix[, "panAcH2B_TSA"] / (signal_matrix[, "panAcH2B_Vehicle"] + 1e-6)

# -------------------------------
# 7. Generate individual heatmaps for each histone mark
# -------------------------------
for(i in 1:4){
  mark <- colnames(fc_matrix)[i]
  log_fc <- log2(fc_matrix[,i] + 1)
  
  # Select top 50 hyperacetylated and hypoacetylated genes
  top_hyper <- names(sort(log_fc, decreasing=TRUE))[1:50]
  top_hypo <- names(sort(log_fc, decreasing=FALSE))[1:50]
  top_genes <- c(top_hyper, top_hypo)
  
  heatmap_matrix <- matrix(log_fc[top_genes], ncol=1)
  rownames(heatmap_matrix) <- top_genes
  colnames(heatmap_matrix) <- mark
  
  pheatmap(heatmap_matrix,
           cluster_rows=TRUE,
           cluster_cols=FALSE,
           main=paste0("Log2 Fold Change TSA vs Vehicle - ", mark),
           color = colorRampPalette(c("blue","white","red"))(100),
           show_rownames=TRUE,
           show_colnames=TRUE)
}

# -------------------------------
# 8.Save results
# -------------------------------
write.csv(signal_matrix, file="ChIPSeq_promoter_signals.csv")
write.csv(fc_matrix, file="ChIPSeq_TSA_vs_Vehicle_FC.csv")
