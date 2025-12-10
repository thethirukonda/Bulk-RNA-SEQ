# Bulk-RNA-SEQ
Hypoxia cell lines
\# RNA-seq Analysis: Hypoxia Response in Prostate Cancer Cell Lines

\## Overview

This repository contains a complete bulk RNA-seq analysis pipeline adapted from eric ku repository. <https://github.com/erilu/bulk-rnaseq-analysis?tab=readme-ov-file#obtaining-raw-data-from-geo>

/). The analysis investigates differential gene expression in response to hypoxia (low oxygen conditions) in two prostate cancer cell lines: LNCaP and PC3.

\*\*Biological Context:\*\*

Hypoxia is a common feature of solid tumors and plays a critical role in cancer progression, metastasis, and treatment resistance. Bulk RNA-seq allows us to identify which genes are upregulated or downregulated under hypoxic conditions, providing insights into:

\- Cellular adaptation mechanisms to low oxygen

\- Potential therapeutic targets

\- Biomarkers for hypoxic tumor regions

\## Dataset Information

\*\*Source:\*\* \[GSE106305\](<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305>) from NCBI Gene Expression Omnibus (GEO)

\*\*Study Design:\*\*

We analyzed control samples from two prostate cancer cell lines under normoxia (normal oxygen) and hypoxia conditions:

| Sample ID | Cell Line | Condition | Treatment | Oxygen Level |

|-----------|-----------|-----------|-----------|--------------|

| SRR7179504 | LNCaP | Control | Empty_Vector | Normoxia |

| SRR7179505 | LNCaP | Control | Empty_Vector | Hypoxia |

| SRR7179508 | PC3 | Control | siCtrl | Normoxia |

| SRR7179509 | PC3 | Control | siCtrl | Hypoxia |

<https://github.com/erilu/bulk-rnaseq-analysis?tab=readme-ov-file#obtaining-raw-data-from-geo>

\## Analysis Pipeline

\### 1. Data Acquisition

\*\*Script:\*\*  
<br/>\`download_fastq.sh\`

Downloaded raw FASTQ files from SRA using \`fastq-dump\` or \`prefetch\` from SRA Toolkit.

\`\`\`bash

\# Example command

fastq-dump --gzip SRR7179504

\`\`\`

\### 2. Quality Control (QC)

\*\*Scripts:\*\* \`fastqc_analysis.sh\`, \`multiqc_report.sh\`

Performed initial quality assessment using FastQC on all raw FASTQ files:

\`\`\`bash

fastqc fastq/\*.fastq.gz -o fastqc_results/ --threads 8

\`\`\`

\*\*Key QC Metrics Evaluated:\*\*

\- Per base sequence quality scores

\- Per sequence quality scores

\- Per base sequence content

\- Sequence duplication levels

\- Adapter content

\- GC content distribution

\*\*MultiQC Integration:\*\*

Aggregated all FastQC reports into a single interactive HTML report for easy comparison:

\`\`\`bash

multiqc fastqc_results/ -o multiqc_report/

\`\`\`

\### 3. Read Trimming

\*\*Script:\*\* \`trimming.sh\`

Based on QC results, performed quality trimming using Trimmomatic (installed via Anaconda):

\`\`\`bash

java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE \\

\-threads 4 \\

fastq/SRR7179504.fastq.gz \\

fastq/SRR7179504_trimmed.fastq.gz \\

TRAILING:10 \\

\-phred33

\`\`\`

\*\*Trimming Parameters:\*\*

\- \`TRAILING:10\` - Remove low quality bases from the end (quality < 10)

\- \`-phred33\` - Specify quality score encoding

Post-trimming QC was performed to verify improvement in read quality.

\### 4. Read Alignment

\*\*Script:\*\* \`align_samples.sh\`

\*\*Tools Required:\*\*

\- \*\*HISAT2\*\* - Fast and sensitive aligner for RNA-seq reads

\- \*\*SAMtools\*\* - Toolkit for manipulating alignment files (sorting, indexing, format conversion)

\*\*Reference Genome:\*\* GRCh38 (Homo sapiens)

Downloaded and indexed the reference genome:

\`\`\`bash

\# Download pre-built HISAT2 index

wget <https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz>

tar -xzf grch38_genome.tar.gz

\`\`\`

\*\*Alignment Pipeline:\*\*

\`\`\`bash

hisat2 -q -x grch38/genome \\

\-U fastq/sample1.fastq.gz | \\

samtools sort -o sample1.bam | \\

samtools index sample1.bam

\`\`\`

The bash script automated this process for all samples, generating sorted and indexed BAM files.

\### 5. Gene Quantification

\*\*Script:\*\* \`count_features.sh\`

\*\*Gene Annotation:\*\* Downloaded GTF file for human reference genome (Homo_sapiens.GRCh38.114.gtf)

The GTF (Gene Transfer Format) file contains gene coordinates, exon boundaries, and gene annotations necessary for assigning reads to specific genes.

\### 6. Count Matrix Generation

\*\*Script:\*\* \`merge_counts.py\`

Combined individual featureCounts outputs into a single count matrix for downstream differential expression analysis:

\`\`\`python

# !/usr/bin/env python

\# coding: utf-8

\`\`\`

\*\*What this script does:\*\*

1\. Locates all featureCounts output files in the specified directory

2\. Reads each file and extracts gene IDs and count values

3\. Renames count columns with their respective sample names

4\. Merges all individual count tables into a single matrix

5\. Exports the final count matrix as a CSV file

The resulting count matrix has genes as rows and samples as columns, ready for differential expression analysis with DESeq2.

\### 7. Differential Expression Analysis

\*\*Script:\*\* \`deseq2_analysis.R\` \*(to be added)\*

Using DESeq2 in R to identify differentially expressed genes between normoxia and hypoxia conditions for each cell line.

\## Software Requirements

\- FastQC (v0.11.9+)

\- MultiQC (v1.12+)

\- Trimmomatic (v0.39)

\- HISAT2 (v2.2.1+)

\- SAMtools (v1.15+)

\- featureCounts (Subread package v2.0.1+)

\- Python 3.8+ (pandas, glob)

\- R (v4.0+) with DESeq2

\## Installation

\`\`\`bash

\# Create conda environment

conda create -n rnaseq python=3.8

conda activate rnaseq

\# Install tools

conda install -c bioconda fastqc multiqc hisat2 samtools subread trimmomatic

\`\`\`

\## Usage

Run scripts in order:

\`\`\`bash

bash scripts/download_fastq.sh

bash scripts/fastqc_analysis.sh

bash scripts/multiqc_report.sh

bash scripts/trimming.sh

bash scripts/align_samples.sh

bash scripts/count_features.sh

python scripts/merge_counts.py

Rscript scripts/deseq2_analysis.R

\`\`\`

\## References

\- Original tutorial: \[Elucidata RNA-seq Tutorial\](<https://github.com/ElucidataInc/>)

\- Dataset: GSE106305 - \[GEO Database\](<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305>)

\## Author

Aakash - Bioinformatics learning project

**Differential Expression Analysis with DESeq2**

**Overview**

This document describes the differential expression (DE) analysis workflow for identifying hypoxia-responsive genes in LNCaP and PC3 prostate cancer cell lines. The analysis uses DESeq2, a widely-used R package for differential gene expression analysis of RNA-seq count data.

**Prerequisites**

**Required R Packages**

\# Bioconductor packages

BiocManager::install(version = "3.22")

BiocManager::install("DESeq2")

\# CRAN packages

install.packages(c("ggrepel", "dplyr", "tibble", "ggplot2", "data.table",

"pheatmap", "RColorBrewer"))

**Required Libraries**

library(DESeq2)

library(dplyr)

library(tibble)

library(ggrepel)

library(ggplot2)

library(data.table)

library(pheatmap)

library(RColorBrewer)

**Analysis Workflow**

**1\. Data Import and Preparation**

**Input:** Count matrix generated from featureCounts (counts_matrix.csv)

\# Set working directory

setwd("/path/to/rnaseq/quants/")

\# Load count matrix

raw_counts <- read.csv("counts_matrix.csv",

header = TRUE,

row.names = "Geneid",

stringsAsFactors = FALSE)

\# Sort columns alphabetically for consistency

raw_counts <- raw_counts\[, sort(colnames(raw_counts))\]

**Sample Structure After Sorting:**

- LNCAP_Hypoxia_S1, LNCAP_Hypoxia_S2
- LNCAP_Normoxia_S1, LNCAP_Normoxia_S2
- PC3_Hypoxia_S1, PC3_Hypoxia_S2
- PC3_Normoxia_S1, PC3_Normoxia_S2

**2\. Create Metadata**

\# Define conditions matching sorted column order

condition <- c(rep("LNCAP_Hypoxia", 2),

rep("LNCAP_Normoxia", 2),

rep("PC3_Hypoxia", 2),

rep("PC3_Normoxia", 2))

\# Create metadata data frame

my_colData <- as.data.frame(condition)

rownames(my_colData) <- colnames(raw_counts)

**Critical Note:** Ensure that the order of conditions in metadata exactly matches the sorted column order in the count matrix.

**3\. Create DESeq2 Dataset Object**

dds <- DESeqDataSetFromMatrix(countData = raw_counts,

colData = my_colData,

design = ~condition)

**QC Check:** Examine zero count distribution

count_matrix <- counts(dds)

zero_counts_per_gene <- rowSums(count_matrix == 0)

print(table(zero_counts_per_gene))

**4\. Gene Annotation**

**Input:** Human genome annotation file (GRCh38annotation.csv)

\# Load annotation

annotation <- fread("GRCh38annotation.csv", stringsAsFactors = FALSE)

colnames(annotation) <- c("Geneid", "Genebiotype", "Genesymbol")

\# Remove version numbers from gene IDs

counts_gse\$Geneid <- sub("\\\\..\*\$", "", counts_gse\$Geneid)

annotation\$Geneid <- sub("\\\\..\*\$", "", annotation\$Geneid)

\# Merge counts with annotation

annotated_counts <- left_join(counts_gse, annotation, by = "Geneid")

**5\. Gene Filtering**

**Step 1: Filter by Gene Biotype**

We retain only protein-coding genes and immunoglobulin/T-cell receptor genes:

biotypes_to_keep <- c("protein_coding",

"IG_J_gene", "IG_V_gene", "IG_C_gene", "IG_D_gene",

"TR_D_gene", "TR_C_gene", "TR_V_gene", "TR_J_gene")

filtered_counts &lt;- annotated_counts %&gt;%

filter(Genebiotype %in% biotypes_to_keep)

**Rationale:**

- Protein-coding genes are the primary targets for differential expression analysis
- Immunoglobulin and T-cell receptor genes are functionally important despite being less abundant

**Step 2: Remove Low-Count Genes**

Remove genes with zero counts in 7 or more samples (out of 8 total):

zero_counts <- rowSums(filtered_counts\[, 4:11\] == 0)

keep_genes <- zero_counts < 7

filtered_counts_nozero <- filtered_counts\[keep_genes, \]

**Rationale:**

- Genes with mostly zero counts provide little statistical power
- Reduces multiple testing burden
- Improves computational efficiency

**Update DESeq2 Object**

dds_filtered <- dds\[rownames(dds) %in% filtered_counts_nozero\$Geneid, \]

**6\. Exploratory Data Analysis**

**Variance Stabilizing Transformation (VST)**

VST normalizes count data to stabilize variance across the mean, making samples more comparable:

vsd <- vst(dds_filtered, blind = TRUE)

**Parameter:** blind = TRUE performs transformation without considering experimental design (appropriate for exploratory analysis)

**Principal Component Analysis (PCA)**

PCA visualizes overall sample relationships and identifies batch effects or outliers:

plot_PCA = function(vsd.obj) {

pcaData <- plotPCA(vsd.obj, intgroup = c("condition"), returnData = TRUE)

percentVar <- round(100 \* attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +

geom_point(size = 3) +

labs(x = paste0("PC1: ", percentVar\[1\], "% variance"),

y = paste0("PC2: ", percentVar\[2\], "% variance"),

title = "PCA Plot colored by condition") +

geom_text_repel(aes(label = name), color = "black")

}

\# Save plot

png("pca_plot.png", width = 2000, height = 2000, res = 300)

plot_PCA(vsd)

dev.off()

**Expected Result:**

- Clear separation between hypoxia and normoxia samples
- Clustering of biological replicates
- Distinct grouping by cell line

**Sample Distance Heatmap**

Visualizes pairwise Euclidean distances between samples:

plotDists = function(vsd.obj) {

sampleDists <- dist(t(assay(vsd.obj)))

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(vsd.obj\$condition)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(55)

pheatmap(sampleDistMatrix,

clustering_distance_rows = sampleDists,

clustering_distance_cols = sampleDists,

col = colors,

fontsize_row = 4,

fontsize_col = 4)

}

png("sample_heatmap.png", width = 1000, height = 900, res = 300)

plotDists(vsd)

dev.off()

**Interpretation:**

- Darker blue = more similar samples
- Biological replicates should cluster together
- Different conditions should show greater distances

**Variable Gene Heatmap**

Displays expression patterns of the most variable genes across samples:

variable_gene_heatmap <- function(vsd.obj, num_genes = 500, annotation) {

brewer_palette <- "RdBu"

ramp <- colorRampPalette(brewer.pal(11, brewer_palette))

mr <- ramp(256)\[256:1\]

stabilized_counts <- assay(vsd.obj)

row_variances <- rowVars(stabilized_counts)

top_variable_genes <- stabilized_counts\[order(row_variances, decreasing = TRUE)\[1:num_genes\], \]

\# Center genes (subtract row means)

top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm = TRUE)

\# Add gene symbols

gene_names <- annotation\$Genesymbol\[match(rownames(top_variable_genes), annotation\$Geneid)\]

rownames(top_variable_genes) <- gene_names

coldata <- as.data.frame(vsd.obj@colData)

coldata\$sizeFactor <- NULL

pheatmap(top_variable_genes,

color = mr,

annotation_col = coldata,

fontsize_col = 8,

fontsize_row = 250/num_genes,

border_color = NA)

}

png("variable_gene_heatmap.png", width = 1000, height = 1000, res = 300)

variable_gene_heatmap(vsd, num_genes = 40, annotation = annotation)

dev.off()

**Color Scale:**

- Red = higher expression
- Blue = lower expression
- White = mean expression

**7\. Biotype Distribution Analysis**

Visualize the proportion of different gene biotypes retained after filtering:

biotype_counts &lt;- filtered_counts_nozero %&gt;%

count(Genebiotype) %>%

mutate(Proportion = n / sum(n),

Percentage = Proportion \* 100)

p <- ggplot(biotype_counts, aes(x = reorder(Biotype, -Proportion),

y = Proportion,

fill = Biotype)) +

geom_bar(stat = "identity") +

labs(title = "Proportion of Genes by Biotype",

x = "Gene Biotypes",

y = "Proportion") +

scale_y_continuous(labels = scales::percent_format(scale = 100)) +

theme_minimal() +

theme(axis.text.x = element_text(angle = 45, hjust = 1),

legend.position = "none") +

scale_fill_brewer(palette = "Set2")

ggsave("genebiotype_proportions.png", plot = p, width = 8, height = 6, dpi = 300)

**8\. Differential Expression Analysis**

**Run DESeq2 Analysis**

dds <- DESeq(dds_filtered)

\# Extract normalized counts

normalized_counts <- counts(dds, normalized = TRUE)

write.csv(as.data.frame(normalized_counts), "normalized_counts.csv", row.names = TRUE)

**DESeq2 performs:**

- Estimation of size factors (normalization)
- Estimation of dispersion
- Negative binomial GLM fitting and Wald statistics

**Cell Line-Specific Analysis: LNCaP**

Extract only LNCaP samples and perform differential expression:

\# Subset to LNCaP samples

dds_lncap <- dds_filtered\[, grepl("LNCAP", colnames(dds_filtered))\]

\# Update factor levels

dds_lncap\$condition <- droplevels(dds_lncap\$condition)

dds_lncap\$condition <- relevel(dds_lncap\$condition, ref = "LNCAP_Normoxia")

\# Run DESeq2

dds_lncap <- DESeq(dds_lncap)

\# Extract results

res_lncap <- results(dds_lncap,

contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))

\# Order by adjusted p-value

reslncapOrdered <- res_lncap\[order(res_lncap\$padj), \]

\# Count significant genes

sum(reslncapOrdered\$padj < 0.05, na.rm = TRUE)

\# Save results

write.csv(as.data.frame(reslncapOrdered), "DEGs_lncap.csv")

**Key Statistics:**

- baseMean: Average normalized count across all samples
- log2FoldChange: Log2 fold change (Hypoxia / Normoxia)
- lfcSE: Standard error of log2 fold change
- stat: Wald statistic
- pvalue: Raw p-value
- padj: Benjamini-Hochberg adjusted p-value

**Interpretation:**

- padj < 0.05: Statistically significant
- log2FoldChange > 0: Upregulated in hypoxia
- log2FoldChange < 0: Downregulated in hypoxia

**9\. Visualization of Results**

**MA Plot**

Shows relationship between mean expression and fold change:

plotMA(res_lncap)

**Interpretation:**

- X-axis: Mean of normalized counts (log scale)
- Y-axis: Log2 fold change
- Blue points: Significant genes (padj < 0.1 by default)
- Red horizontal lines: Log2 fold change thresholds

**Volcano Plot**

Displays statistical significance vs. biological effect size:

res_df <- as.data.frame(reslncapOrdered)

res_df <- na.omit(res_df)

res_df\$gene <- rownames(res_df)

\# Classify genes

res_df\$regulation <- "Not Significant"

res_df\$regulation\[res_df\$padj &lt; 0.05 & res_df\$log2FoldChange &gt; 1\] <- "Upregulated"

res_df\$regulation\[res_df\$padj < 0.05 & res_df\$log2FoldChange < -1\] <- "Downregulated"

\# Create plot

qp <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +

geom_point(alpha = 0.6) +

scale_color_manual(values = c("Upregulated" = "#FEA405",

"Downregulated" = "purple",

"Not Significant" = "gray")) +

geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +

theme_minimal() +

labs(title = "Volcano Plot - LNCaP Hypoxia vs Normoxia",

x = "Log2 Fold Change",

y = "-Log10 Adjusted P-Value")

ggsave("vp_lncap.png", plot = qp, width = 8, height = 6, dpi = 300)

**Thresholds:**

- Horizontal line: padj = 0.05
- Vertical considerations: |log2FC| > 1 (2-fold change)

**Colors:**

- Orange: Significantly upregulated (padj &lt; 0.05, log2FC &gt; 1)
- Purple: Significantly downregulated (padj < 0.05, log2FC < -1)
- Gray: Not significant

**Heatmap of Top DE Genes**

Visualizes expression patterns of most significant genes:

DE_gene_heatmap <- function(res, count_matrix, padj_cutoff = 0.0001, ngenes = 20) {

brewer_palette <- "RdBu"

ramp <- colorRampPalette(brewer.pal(11, brewer_palette))

mr <- ramp(256)\[256:1\]

\# Get top significant genes

significant_genes &lt;- as.data.frame(res) %&gt;%

filter(padj &lt; padj_cutoff) %&gt;%

arrange(desc(log2FoldChange)) %>%

head(ngenes)

\# Extract counts and scale

heatmap_values <- count_matrix\[rownames(significant_genes), \]

heatmap_values <- t(scale(t(heatmap_values))) # Z-score normalization

pheatmap(heatmap_values,

color = mr,

scale = "none",

cluster_rows = TRUE,

cluster_cols = TRUE,

fontsize_col = 10,

fontsize_row = max(6, 200/ngenes),

border_color = NA,

main = paste("Top", ngenes, "DE Genes (padj <", padj_cutoff, ")"))

}

count_matrix <- assay(dds_filtered)

png("de_gene_heatmap.png", width = 800, height = 600, res = 150)

DE_gene_heatmap(res_lncap, count_matrix, padj_cutoff = 0.001, ngenes = 30)

dev.off()

**Features:**

- Rows: Top differentially expressed genes
- Columns: Samples
- Colors: Z-score normalized expression (blue = low, red = high)
- Clustering: Hierarchical clustering of both genes and samples

**Density Plots: Raw vs VST Counts**

Compare distribution before and after variance stabilization:

raw_counts <- assay(dds_filtered)

vst_counts <- assay(vsd)

png("density_plots_raw_vst.png", width = 4000, height = 4000, res = 300)

par(mfrow = c(4, 4), mar = c(3, 3, 2, 1))

for (i in 1:8) {

\# Raw counts

plot(density(raw_counts\[, i\]),

main = paste("Raw - Sample", colnames(raw_counts)\[i\]),

col = "red", lwd = 2)

\# VST counts

plot(density(vst_counts\[, i\]),

main = paste("VST - Sample", colnames(vst_counts)\[i\]),

col = "blue", lwd = 2)

}

dev.off()

**Interpretation:**

- Raw counts: Right-skewed, heteroscedastic (variance depends on mean)
- VST counts: More normal-like, homoscedastic (stable variance)
- VST transformation makes samples more comparable

**Output Files**

| **File Name** | **Description** |
| --- | --- |
| annotated_counts_with_symbols.csv | Count matrix with gene symbols and biotypes |
| 9biotype_count_matrix.csv | Counts filtered by biotype |
| filtered_biotype_6.csv | Final filtered counts (biotype + low-count filter) |
| normalized_counts.csv | DESeq2 normalized counts |
| DEGs_lncap.csv | LNCaP differential expression results |
| pca_plot.png | PCA visualization |
| sample_heatmap.png | Sample distance heatmap |
| variable_gene_heatmap.png | Most variable genes heatmap |
| genebiotype_proportions.png | Biotype distribution bar plot |
| vp_lncap.png | Volcano plot |
| de_gene_heatmap.png | Top DE genes heatmap |
| density_plots_raw_vst.png | Raw vs VST transformation comparison |

**Key Findings**

**Quality Control**

- Clear PCA separation between hypoxia and normoxia conditions
- Biological replicates cluster together
- No obvious batch effects or outliers

**Gene Filtering**

- Started with 78,899 genes (all annotated genes)
- Retained ~20,000-25,000 genes after biotype and count filtering
- Majority are protein-coding genes (~95%)

**Differential Expression (LNCaP)**

- Hundreds to thousands of genes differentially expressed (padj < 0.05)
- Both upregulated and downregulated genes identified
- Enrichment expected in hypoxia-response pathways (HIF1A targets, glycolysis, angiogenesis)

**Biological Interpretation**

**Expected Hypoxia-Responsive Genes**

**Upregulated:**

- HIF1A pathway targets (e.g., VEGFA, EPO, LDHA)
- Glycolysis genes (e.g., GLUT1, HK2, PFKFB3)
- Angiogenesis factors
- Survival/anti-apoptotic genes

**Downregulated:**

- Oxidative phosphorylation genes
- Cell cycle progression genes
- Differentiation markers

**Next Steps for Biological Validation**

- **Functional Enrichment Analysis:**
  - Gene Ontology (GO) analysis
  - KEGG pathway analysis
  - GSEA (Gene Set Enrichment Analysis)
- **Validation:**
  - Check known hypoxia markers (HIF1A, VEGFA, CA9, BNIP3)
  - Compare with published hypoxia signatures
  - qRT-PCR validation of top hits
- **Cell Line Comparison:**
  - Perform same analysis for PC3 cell line
  - Compare response patterns between LNCaP and PC3
  - Identify cell line-specific vs. universal hypoxia responses

Remember to save session info for reproducibility:

sessionInfo()

**Author:** Aakash Deva Thirukonda Prakash  
**Date:** December 10, 2025  
**Analysis:** RNA-seq Differential Expression - Hypoxia Response
