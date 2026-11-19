# ==============================================================================
# Differential Expression Analysis with DESeq2 (WT vs. DKO, Control vs. Case)
# Includes both group contrasts (~ Group) and genotype x infection interaction
# (~ genotype + status + genotype:status)
# ==============================================================================

# Install packages if not already installed:
# install.packages(c("ggplot2", "dplyr", "tibble", "pheatmap", "RColorBrewer"))
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)

# ---- 1. Paths (Use forward slashes "/" to avoid Windows escape errors) ----
counts_file <- "C:/Users/gmont/Desktop/RNASeqData/featureCounts/gene_counts.txt"
out_dir     <- "C:/Users/gmont/Desktop/RNASeqData/DESeq2Contrast"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load featureCounts Table ----
raw_counts <- read.delim(counts_file, comment.char = "#", check.names = FALSE)

count_mat <- as.matrix(raw_counts[, 7:ncol(raw_counts)])
rownames(count_mat) <- raw_counts$Geneid

# Clean column names so they match Sample IDs (e.g., SRR7821949)
colnames(count_mat) <- sub("\\.sorted\\.bam$", "", basename(colnames(count_mat)))
storage.mode(count_mat) <- "integer"

# ---- 3. Sample-to-Condition Metadata ----
# Case = Infected, Control = Uninfected
# WT = Wildtype, DKO = Double-Knockout
sample_info <- data.frame(
  Sample = c(
    "SRR7821949", "SRR7821950", "SRR7821951", "SRR7821952", "SRR7821953",
    "SRR7821968", "SRR7821969", "SRR7821970",
    "SRR7821954", "SRR7821955", "SRR7821956", "SRR7821957",
    "SRR7821971", "SRR7821972", "SRR7821973"
  ),
  Group = c(
    rep("Blood_WT_Case", 5),
    rep("Blood_WT_Control", 3),
    rep("Blood_DKO_Case", 4),
    rep("Blood_DKO_Control", 3)
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    genotype = factor(
      ifelse(grepl("_WT_", Group), "WT", "DKO"),
      levels = c("WT", "DKO")
    ),
    status = factor(
      ifelse(grepl("_Case$", Group), "Case", "Control"),
      levels = c("Control", "Case")
    ),
    Group = factor(
      Group,
      levels = c("Blood_WT_Control", "Blood_WT_Case", "Blood_DKO_Control", "Blood_DKO_Case")
    )
  )

rownames(sample_info) <- sample_info$Sample

# Align count matrix columns to match sample_info row order
missing_samples <- setdiff(rownames(sample_info), colnames(count_mat))
if (length(missing_samples) > 0) {
  stop("Missing count columns for samples: ", paste(missing_samples, collapse = ", "))
}
count_mat <- count_mat[, rownames(sample_info), drop = FALSE]
stopifnot(identical(colnames(count_mat), rownames(sample_info)))

# ---- 4. Group Model (~ Group): Pre-Filter and Run DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = ~ Group
)

# Keep genes with at least 10 counts in at least 3 samples (smallest group size is 3)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)

save_contrast <- function(dds_obj, numerator_grp, denominator_grp, prefix) {
  res_obj <- results(dds_obj, contrast = c("Group", numerator_grp, denominator_grp), alpha = 0.05)
  
  res_tbl <- as.data.frame(res_obj) %>%
    rownames_to_column(var = "gene_id") %>%
    arrange(padj, pvalue)
  
  sig_tbl <- res_tbl %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  
  write.csv(res_tbl, file.path(out_dir, paste0("deseq2_all_", prefix, ".csv")), row.names = FALSE)
  write.csv(sig_tbl, file.path(out_dir, paste0("deseq2_sig_", prefix, ".csv")), row.names = FALSE)
  
  cat(sprintf("[%s] Total genes: %d | Sig DEGs (padj < 0.05 & |log2FC| > 1): %d\n",
              prefix, nrow(res_tbl), nrow(sig_tbl)))
  
  invisible(list(all = res_tbl, sig = sig_tbl))
}

wt_inf    <- save_contrast(dds, "Blood_WT_Case", "Blood_WT_Control", "WT_Case_vs_WT_Control")
dko_inf   <- save_contrast(dds, "Blood_DKO_Case", "Blood_DKO_Control", "DKO_Case_vs_DKO_Control")
ctrl_geno <- save_contrast(dds, "Blood_DKO_Control", "Blood_WT_Control", "DKO_Control_vs_WT_Control")
case_geno <- save_contrast(dds, "Blood_DKO_Case", "Blood_WT_Case", "DKO_Case_vs_WT_Case")

# ---- 4b. Interaction Model (~ genotype + status + genotype:status) ----
dds_int <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = ~ genotype + status + genotype:status
)
dds_int <- dds_int[keep, ]
dds_int <- DESeq(dds_int)

print(resultsNames(dds_int))

# Interaction term: difference in infection response (Case vs. Control) between DKO and WT
res_int <- results(dds_int, name = "genotypeDKO.statusCase", alpha = 0.05)

int_all <- as.data.frame(res_int) %>%
  rownames_to_column(var = "gene_id") %>%
  arrange(padj, pvalue)

int_sig <- int_all %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

write.csv(int_all, file.path(out_dir, "deseq2_all_interaction_genotypeDKO_x_statusCase.csv"), row.names = FALSE)
write.csv(int_sig, file.path(out_dir, "deseq2_sig_interaction_genotypeDKO_x_statusCase.csv"), row.names = FALSE)

cat(sprintf("[Interaction genotypeDKO:statusCase] Total genes: %d | Sig DEGs (padj < 0.05 & |log2FC| > 1): %d\n",
            nrow(int_all), nrow(int_sig)))

# Optional: Likelihood Ratio Test (LRT) comparing full interaction model vs. additive model
dds_lrt <- DESeq(dds_int, test = "LRT", reduced = ~ genotype + status)
res_lrt <- results(dds_lrt, alpha = 0.05)
lrt_tbl <- as.data.frame(res_lrt) %>%
  rownames_to_column(var = "gene_id") %>%
  arrange(padj, pvalue)
write.csv(lrt_tbl, file.path(out_dir, "deseq2_LRT_interaction.csv"), row.names = FALSE)

# ---- 5. PCA and Sample Distance Heatmap Colored by Group ----
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_data$Group <- factor(pca_data$Group, levels = levels(sample_info$Group))
pca_data$genotype <- factor(
  ifelse(grepl("_WT_", pca_data$Group), "WT", "DKO"),
  levels = c("WT", "DKO")
)
pca_data$status <- factor(
  ifelse(grepl("_Case$", pca_data$Group), "Case", "Control"),
  levels = c("Control", "Case")
)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, shape = genotype, label = name)) +
  geom_point(size = 3.5) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of Variance-Stabilized Counts by Group") +
  theme_bw()

ggsave(file.path(out_dir, "pca_by_group.png"), p_pca, width = 8.5, height = 6.5, dpi = 300)

sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- paste(vsd$Sample, vsd$Group, sep = " | ")
colnames(sample_dist_mat) <- NULL

colors_hm <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
png(file.path(out_dir, "sample_distance_heatmap.png"), width = 2200, height = 1900, res = 300)
pheatmap(
  sample_dist_mat,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col = colors_hm,
  main = "Sample-to-Sample Euclidean Distances (VST)"
)
dev.off()

# ---- 6. Volcano Plot for Interaction Contrast ----
volcano_int_df <- int_all %>%
  mutate(
    neg_log10_padj = -log10(ifelse(is.na(padj) | padj == 0, 1e-300, padj)),
    regulation = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 1  ~ "Stronger Case Induction in DKO",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Weaker Case Induction in DKO",
      TRUE                                            ~ "Not Significant"
    )
  )

p_volcano_int <- ggplot(volcano_int_df, aes(x = log2FoldChange, y = neg_log10_padj, color = regulation)) +
  geom_point(alpha = 0.65, size = 1.5) +
  scale_color_manual(values = c(
    "Stronger Case Induction in DKO" = "#D95F02",
    "Weaker Case Induction in DKO"   = "#2B5C8F",
    "Not Significant"                = "#BDBDBD"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    title = "Interaction Volcano: Infection Response in DKO vs. WT",
    x = expression(log[2] ~ "Interaction Fold Change (DKO Case/Control vs. WT Case/Control)"),
    y = expression(-log[10] ~ "(adjusted p-value)"),
    color = "Interaction"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "volcano_interaction_genotype_x_status.png"), p_volcano_int, width = 9, height = 6.5, dpi = 300)

# ---- 7. Top 30 Interaction Genes Heatmap Across All Groups ----
top_int_genes <- head(int_sig$gene_id, 30)
if (length(top_int_genes) >= 2) {
  mat_int <- assay(vsd)[top_int_genes, , drop = FALSE]
  mat_int_scaled <- t(scale(t(mat_int)))
  annot_col <- data.frame(
    Group = sample_info$Group,
    row.names = rownames(sample_info)
  )
  
  png(file.path(out_dir, "top30_interaction_heatmap.png"), width = 2300, height = 2600, res = 300)
  pheatmap(
    mat_int_scaled,
    annotation_col = annot_col,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = "Top 30 Genotype x Infection Interaction Genes (VST Z-Score)"
  )
  dev.off()
}
