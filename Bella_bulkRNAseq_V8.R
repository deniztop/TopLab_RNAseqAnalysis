# Function to install and load packages automatically
install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}
# Core analysis packages
install_load("DESeq2")
install_load("tidyverse") # Includes ggplot2, dplyr, readr, etc.
install_load("pheatmap")  # For heatmaps
install_load("Rsubread")  # Added for FASTQ alignment and counting
install_load("clusterProfiler") # For GSEA analysis
install_load("org.Dm.eg.db")    # Drosophila Genome Annotation

# Load packages if already installed
library("DESeq2")
library("tidyverse") # Includes ggplot2, dplyr, readr, etc.
library("pheatmap")  # For heatmaps
library("Rsubread")  # Added for FASTQ alignment and counting
library("clusterProfiler") # For GSEA analysis
library("org.Dm.eg.db")    # Drosophila Genome Annotation


#-------------------------------------------------------------------------------
# 1. Generate meta data sheet and counts matrix to use for analysis
#-------------------------------------------------------------------------------

# --- 1. Meta Data Sheet ---

# Create meta data file: Rows are samples, columns are conditions/factors.
# Create the Metadata Table immediately based on your design
# We create a 2x2 factorial design: Genotype (WT//S15A) and Timepoint (ZT2/ZT14)
col_data <- data.frame(
  row.names = sample_names,
  Genotype = factor(rep(c("wt", "wt", "S15A", "S15A",), each=4), levels = c("wt", "S15A",)),
  Timepoint = factor(rep(c("ZT2", "ZT14", "ZT2", "ZT14"), each=4), levels = c("ZT2", "ZT14"))
)

# Combine into a single Group factor for easier plotting/contrasts later - this combines the two factors (genotype and timepoint) into one group
col_data$Group <- factor(paste(col_data$Genotype, col_data$Timepoint, sep="_"))

# --- 2. Define Files and Groups ---
# NOTE: Ensure all 32 bam files are in your working directory

bam_files <- c("Top-1_DT_012726.bam", 
                 "Top-2_DT_012726.bam", 
                 "Top-3_DT_012726.bam", 
                 "Top-4_DT_012726.bam", 
                 "Top-5_DT_012726.bam", 
                 "Top-6_DT_012726.bam", 
                 "Top-7_DT_012726.bam", 
                 "Top-8_DT_012726.bam", 
                 "Top-9_DT_012726.bam", 
                 "Top-10_DT_012726.bam", 
                 "Top-11_DT_012726.bam", 
                 "Top-12_DT_012726.bam", 
                 "Top-13_DT_012726.bam", 
                 "Top-14_DT_012726.bam", 
                 "Top-15_DT_012726.bam", 
                 "Top-16_DT_012726.bam",)

# Define clean sample names for the columns
sample_names <- c("wt_ZT2_1", "wt_ZT2_2", "wt_ZT2_3", "wt_ZT2_4", "wt_ZT14_1", "wt_ZT14_2", "wt_ZT14_3", "wt_ZT14_4", "S15A_ZT2_1", "S15A_ZT2_2", "S15A_ZT2_3", "S15A_ZT2_4", "S15A_ZT14_1", "S15A_ZT14_2", "S15A_ZT14_3", "S15A_ZT14_4")

# --- 3. Counts Matrix ---

# Genome annotation file used for alignment must be in your wd (i.e. in the same file as your BAM files)
# Count Features (Generate Matrix)
# PairedEnd=TRUE because the reads were paired 
gtf_file <- "Drosophila_melanogaster.BDGP6.22.97.gtf"
fc <- featureCounts(files=bam_files, annot.ext=gtf_file, 
                    isGTFAnnotationFile=TRUE, isPairedEnd=TRUE) 
 
# Extract count matrix and rename the columns with your sample names
counts_data <- fc$counts
colnames(counts_data) <- sample_names

#--- 4. Save Counts and Metadata sheets ---

# Save generated files for future use. The next steps of the script are completed with these two files.
write.csv(counts_data, "generated_counts_matrix.csv")
write.csv(col_data, "generated_metadata.csv")


# ------------------------------------------------------------------------------
# 2. Pre-processing and Sanity Checks
# ------------------------------------------------------------------------------

counts_data <- read.csv("generated_counts_matrix.csv", row.names = 1)
col_data <- read.csv("generated_metadata.csv", row.names = 1)

# Ensure column names in counts match row names in metadata
if (!all(colnames(counts_data) == rownames(col_data))) {
  stop("Error: Column names of count matrix do not match row names of metadata.")
}

# Create DESeq2 dataset object
# Use a design that accounts for both Genotype and Timepoint
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ Timepoint + Genotype)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Additional DEseq2 dataset object to look at S15A vs. wt at each timepoint
# Create a new DESeq object with the design using Group factor
dds_group <- DESeqDataSetFromMatrix(countData = counts_data,
                                    colData = col_data,
                                    design = ~ Group)
# Filter low counts and run dds pipeline
dds_group <- dds_group[rowSums(counts(dds_group)) >= 10, ]


# ------------------------------------------------------------------------------
# 3. Differential Expression Analysis
# ------------------------------------------------------------------------------
# --- Results 1: Effect of Genotype (S15A vs wt) ---
# Run the main DESeq pipeline
dds <- DESeq(dds)

# This averages across both Timepoints
res_geno <- results(dds, contrast=c("Genotype", "S15A", "wt"))
summary(res_geno)

# Order from lowest to highest padj value 
resOrdered <- res_geno[order(res_geno$padj),]
write.csv(as.data.frame(resOrdered), file="deseq2_results_genotype.csv")

# Apply FDR<0.05 filter on results
res_geno_05 <- results(dds, contrast=c("Genotype", "S15A", "wt"), alpha = 0.05)
# Count total DE genes (both up and down)
sum(res_geno_05$padj < 0.05, na.rm = TRUE)

summary(res_geno, alpha = 0.05)

resOrdered_05 <- res_geno_05[order(res_geno_05$padj),]
# Subset to keep ONLY genes with padj < 0.05 (and remove NA values)
resSig_05 <- subset(resOrdered_05, padj < 0.05)

# Export the filtered data frame
write.csv(as.data.frame(resSig_05), file="deseq2_results_genotype_FDR0.05.csv")

# --- Results 2: Group Effect (S15A vs. wt at each timepoint) ---
# Run group design dds piepline
dds_group <- DESeq(dds_group)

# Comparison at ZT2
res_ZT2 <- results(dds_group, contrast=c("Group", "S15A_ZT2", "wt_ZT2"))
summary(res_ZT2)
res_ZT2_ordered <- res_ZT2[order(res_ZT2$padj),]
write.csv(as.data.frame(res_ZT2_ordered), file = "deseq2_results_ZT2_group.csv")

# Comparison at ZT14
res_ZT14 <- results(dds_group, contrast=c("Group", "S15A_ZT14", "wt_ZT14"))
summary(res_ZT14)
res_ZT14_ordered <- res_ZT14[order(res_ZT14$padj),]
write.csv(as.data.frame(res_ZT14_ordered), file = "deseq2_results_ZT14_group.csv")


# ------------------------------------------------------------------------------
# 4. Heatmap
# ------------------------------------------------------------------------------
# --- Heatmap of Top Genes ---
# --- Top 50 DE Genes ---
# Create a matrix with the top 50 most significantly differentially expressed genes
top_genes50 <- head(order(res_geno$padj), 50)
mat50 <- assay(vst_data)[top_genes50, ]

mat50_scaled <- t(scale(t(mat50)))

# Add Gene names from FlyBase query to the matrix
print(rownames(mat50_avg_scaled))
rownames(mat50_scaled) <- c("IntS12", "Mis12", "Unc-115a", "Ir76a", "Ugt35E2", "CG31205", "Fam161", "Pex3", "CG8539", "Prp39", "Obp99d", "hng3", "ppk20", "CG32523", "tld", "CR14798", "CG18542", "CG32939", "ras", "Ir75a", "AstC-R1", "CG3091", "Ar2", "CG42541", "ppk30", "CG13895", "CG14131", "Obp99a", "spartin", "Mrp4", "Lsp1β", "lncRNA:CR32553", "Ctr9", "CR43217", "CG32212", "fit", "Vps13B", "CG10962", "Cnx14D", "CG11052", "CkIIα-i3",
                                "GNBP1",
                                "Herc4",
                                "Loxl1",
                                "Mvb12",
                                "CG2604",
                                "CR45496",
                                "Mdr50",
                                "mRpL39",
                                "CR43215")

# Generate Heatmap
pheatmap(mat50_scaled, 
         annotation_col = df, 
         cluster_cols = FALSE,
         main = "Top 50 DE Genes (Genotype Effect) Z-score Normalized",
         fontsize_row = 10)

#--- Averaging Across 4 Biological Replicates ---
# 50 DE genes
mat50_avg <- t(apply(mat50, 1, function(x) tapply(x, vst_data$Group, mean)))
mat50_avg_scaled <- t(scale(t(mat50_avg)))

# Replace FBgns with gene names in the matrix
rownames(mat50_avg_scaled) <- c("IntS12", "Mis12", "Unc-115a", "Ir76a", "Ugt35E2", "CG31205", "Fam161", "Pex3", "CG8539", "Prp39", "Obp99d", "hng3", "ppk20", "CG32523", "tld", "CR14798", "CG18542", "CG32939", "ras", "Ir75a", "AstC-R1", "CG3091", "Ar2", "CG42541", "ppk30", "CG13895", "CG14131", "Obp99a", "spartin", "Mrp4", "Lsp1β", "lncRNA:CR32553", "Ctr9", "CR43217", "CG32212", "fit", "Vps13B", "CG10962", "Cnx14D", "CG11052", "CkIIα-i3",
                                "GNBP1",
                                "Herc4",
                                "Loxl1",
                                "Mvb12",
                                "CG2604",
                                "CR45496",
                                "Mdr50",
                                "mRpL39",
                                "CR43215")

pheatmap(mat50_avg_scaled[, rownames(df_unique)], annotation_col=df_unique, main="Updated Top 50 DE Genes (Averaged) Z-score Normalized")

write.csv(mat50_avg_scaled, file="Top50_DEgenes_Z-scoreNormalized_average_genenames.csv")


# ------------------------------------------------------------------------------
# 5. Extracting Normalized Counts for genes of interest
# ------------------------------------------------------------------------------
# --- 1. Circadian Genes ---
genes_of_interest <- c("FBgn0023076", "FBgn0023094", "FBgn0003068", "FBgn0014396", "FBgn0016076", "FBgn0016694", "FBgn0025680", "FBgn0259938")

# Extract all normalized counts from DESeq2 object
normalized_counts <- counts(dds, normalized = TRUE)
# Subset the matrix for only circadian
cg_normalized_counts <- normalized_counts[genes_of_interest, ]
# Save as csv
write.csv(as.data.frame(cg_normalized_counts), file = "normalized_counts_circadian_genes.csv")

# --- 2. Supplemental Figure Genes ---
supp_genes <- c("FBgn0032840", "FBgn0027109", "FBgn0039298", "FBgn0023178", "FBgn0035023")

# Normalized counts already is already set as data frame, therefore extract normalized counts for other genes of interest
supp_normalized_counts <- normalized_counts[supp_genes, ]
# Save as csv
write.csv(as.data.frame(supp_normalized_counts), "sup_normalized_counts.csv")


# ------------------------------------------------------------------------------
# 6. Gene Set Enrichment Analysis (GSEA)
# ------------------------------------------------------------------------------
# --- 1. Prepare Ranked List of Genes ---
# Use the Wald statistic (stat) to rank genes. 
res_gsea <- res_geno
res_gsea <- res_gsea[!is.na(res_gsea$stat), ] # Remove NAs

gene_list <- res_gsea$stat
names(gene_list) <- rownames(res_gsea)
gene_list <- sort(gene_list, decreasing = TRUE)

set.seed(12)
# --- 2. Run GSEA (GO: Biological Process) ---
# Run GSEA against Drosophila Database
gse <- gseGO(geneList     = gene_list,
               OrgDb        = org.Dm.eg.db,
               ont          = "BP",        # Biological Process
               keyType      = "ENSEMBL",   # Match your ID type
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE,
               seed         = TRUE)
  
# --- 3. Visualization ---
# Visualizing the Top 10 UP and Top 10 DOWN pathways
top_pathways <- rbind(head(gse@result, 10), tail(gse@result, 10))

# Sort the results by NES (highest to lowest) to ensure head() is UP and tail() is DOWN
gse_res_sorted <- gse@result[order(gse@result$NES, decreasing = TRUE), ]
# Now safely extract the top and bottom 10
top_pathways_sorted <- rbind(head(gse_res_sorted, 10), tail(gse_res_sorted, 10))

#---- Visualize top 20 up and top 20 down pathways
# Filter for statistically significant results first
sig_res <- gse@result[gse@result$p.adjust < 0.05, ]

# Sort by NES
sig_res_sorted <- sig_res[order(sig_res$NES, decreasing = TRUE), ]

# Extract top and bottom 
top_up <- head(sig_res_sorted, 20)
top_down <- tail(sig_res_sorted, 20)
top_pathways_40 <- rbind(top_up, top_down)
# Calculate the core gene count from the core_enrichment string
top_pathways_40$Count <- sapply(strsplit(top_pathways_40$core_enrichment, "/"), length)
top_pathways_40$Description <- factor(top_pathways_40$Description, 
                                          levels = unique(top_pathways_40$Description[order(top_pathways_40$NES)]))

# Generate dot plot
ggplot(top_pathways_40, aes(x = NES, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(
    title = "GSEA Dot Plot - Top 20 UP/DOWN by NES",
    subtitle = "Top GO Biological Processes",
    x = "Normalized Enrichment Score (NES)",
    y = "Gene Ontology Term",
    size = "Gene Count",
    color = "Adjusted p-value"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14)
  )

# Save GSEA Results
write.csv(as.data.frame(gse), "gsea_results.csv")
write.csv(top_pathways_40, file = "top_20_up_down_pathways.csv", row.names = FALSE)
