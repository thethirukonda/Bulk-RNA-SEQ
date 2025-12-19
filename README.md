# README

# Bulk RNA-Seq Data Analysis – Hypoxia Response in Prostate Cancer Cell Lines

## 🎯 Project Overview

I'm just trying out a bulk rna seq pipeline- starting from raw SRA files and culminating in differential gene expression analysis with biological pathway interpretation. This is an adaptation of [Eric Lu’s Tutorial on bulk rna seq](https://github.com/erilu/bulk-rnaseq-analysis?tab=readme-ov-file#obtaining-raw-data-from-geo), which is in turn based on  ONECUT2 is a driver of neuroendocrine prostate cancer by  [Guo et al 2019 work-](https://www.ncbi.nlm.nih.gov/pubmed/30655535) 

The analysis investigates how **hypoxia** (low oxygen conditions) alters gene expression in two prostate cancer cell lines:

**1. LNCaP (Lymph Node Carcinoma of the Prostate)**

- Origin: Derived from a lymph node metastasis of human prostate adenocarcinoma
- Characteristics: Androgen-sensitive, represents hormone-responsive prostate cancer, androgen-dependent 

**2. PC3 (Prostate Cancer 3)**

- Origin: Derived from a bone metastasis of grade IV prostate adenocarcinoma

- Characteristics: Androgen-independent, highly aggressive phenotype

**Research Question:** How does hypoxia affect gene expression in these prostate cancer cell lines, and what biological pathways are activated in response to low oxygen?

---

## 📌 Objectives

- Master the complete workflow of Bulk RNA-Seq analysis from raw data to biological insights
- Perform quality control, alignment, quantification, and differential expression analysis
- Visualize results using industry-standard plots (PCA, volcano plots, heatmaps)
- Conduct pathway enrichment analysis to interpret biological significance
- Document the entire process for reproducibility and educational purposes

---

## 📂 Dataset Information

**Dataset Source:** [GSE106305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305) (NCBI GEO)

**Publication:** Guo et al., Nature Communications, 2019

### Experimental Design

The original study examined multiple conditions, but for this analysis, we focused on **control samples** to investigate the basic hypoxia response:

| Cell Line | Condition | Control Treatment | Oxygen Level | Biological Replicates |
| --- | --- | --- | --- | --- |
| LNCaP | Control | Empty_Vector | Normoxia | 2 |
| LNCaP | Control | Empty_Vector | Hypoxia | 2 |
| PC3 | Control | siCtrl | Normoxia | 2 |
| PC3 | Control | siCtrl | Hypoxia | 2 |

**Total samples analyzed:** 8 (4 conditions × 2 replicates)

### Data Organization in SRA Database

SRA stores the actual raw sequence reads (FASTQ files) from various sequencing experiments.

Raw sequencing data is distributed across multiple technical runs (lanes). Each biological sample was sequenced in multiple runs and stored as separate SRR files:

**Downloaded SRR Files:**

```
SRR7179504, SRR7179505, SRR7179506, SRR7179507  (LNCaP Normoxia Rep1)
SRR7179508, SRR7179509, SRR7179510, SRR7179511  (LNCaP Normoxia Rep2)
SRR7179520, SRR7179521, SRR7179522, SRR7179523  (LNCaP Hypoxia Rep1)
SRR7179524, SRR7179525, SRR7179526, SRR7179527  (LNCaP Hypoxia Rep2)
SRR7179536 (PC3 Normoxia Rep1)
SRR7179537 (PC3 Normoxia Rep2)
SRR7179540 (PC3 Hypoxia Rep1)
SRR7179541 (PC3 Hypoxia Rep2)
```

**Total raw files:** 20 SRR runs → merged into 8 biological samples

 These represent 8 biological samples (2 replicates for each condition).

---

## ⚙️ Analysis Pipeline Overview

Ideally, create a conda environment to isolate, version-control, and ensure reproducibility of all dependencies, tools, and software packages required for the RNA-seq analysis pipeline across different systems.

I didn't create an environment this time, this is something I realised in hindsight, but it is a standard good practice to make it more replicable and ensure cross-platform functionality. 

```
Raw SRA Files
    ↓ [1. Data Download & Conversion]
FASTQ Files (Sequencing Reads)
    ↓ [2. Quality Control]
QC Reports (FastQC/MultiQC)
    ↓ [3. Read Trimming (Optional)]
Cleaned FASTQ Files
    ↓ [4. Sample Merging]
8 Biological Samples
    ↓ [5. Reference Genome Preparation]
Genome Index + GTF Annotation
    ↓ [6. Alignment]
BAM Files (Aligned Reads)
    ↓ [7. Post-Alignment QC]
Alignment Quality Reports
    ↓ [8. Gene Quantification]
Count Matrix (Genes × Samples)
    ↓ [9. Differential Expression Analysis]
DESeq2 Results + Visualizations
    ↓ [10. Pathway Enrichment]
Biological Insights
```

---

## 🔬 Detailed Methodology

### Step 1: Data Download & Conversion

**Tools:** SRA Toolkit (`prefetch`, `fastq-dump`)

**Installation:**

```bash
#download from github wiki page of sratoolkit 

cd rnaseq 

# Download version 3.3.0
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.3.0/sratoolkit.3.3.0-mac-x86_64.tar.gz

# Extract the file you actually downloaded (matching filename)
tar -xzf sratoolkit.3.3.0-mac-x86_64.tar.gz

# Navigate to the correct extracted folder (version 3.3.0, not 3.0.0)
cd sratoolkit.3.3.0-mac-x86_64/bin

# Make it usable
export PATH=$PATH:$(pwd)

# Test it works
./fastq-dump --version
```

This is if you are creating in a conda enviroment, i didnt make in a conda enviroment but for future refrence.

```bash

# 1. Create conda environment
conda create -n rnaseq python=3.9 -y

# 2. Activate the environment
conda activate rnaseq

# 3. Navigate to your rnaseq directory (or create it)
mkdir -p ~/rnaseq
cd ~/rnaseq

# 4. Download SRA Toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.3.0/sratoolkit.3.3.0-mac-x86_64.tar.gz

# 5. Extract it
tar -xzf sratoolkit.3.3.0-mac-x86_64.tar.gz

# 6. Add to conda environment's PATH permanently
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo "export PATH=\"$(pwd)/sratoolkit.3.3.0-mac-x86_64/bin:\$PATH\"" > $CONDA_PREFIX/etc/conda/activate.d/sra-tools.sh

# 7. Reactivate environment to apply changes
conda deactivate
conda activate rnaseq

# 8. Test it works (no need for ./ or full path anymore!)
fastq-dump --version
which fastq-dump

# 9. Install other RNA-seq tools
conda install -c bioconda fastqc star samtools multiqc -y
```

If you have installed you would get an output like this 

```
Read 2 spots for SRR390728
Written 2 spots for SRR390728
@SRR390728.1 1 length=72
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
+SRR390728.1 1 length=72
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
@SRR390728.2 2 length=72
AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
+SRR390728.2 2 length=72
;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
```

This is the **official test** from NCBI. If you see this exact output, the toolkit is installed correctly! 

**Download Example (Single Sample):**

```bash
# Download SRA fileprefetch SRR7179504
# Convert to FASTQ format

prefetch SRR7179504
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra

```

## Automated Download and Conversion

For processing all 20 samples, I used two Python scripts due to a path error in the initial attempt.

### What Happened

**First Script: `download_fastq.py`**

- Successfully **downloaded** all 20 `.sra` files using `prefetch`
- **Failed** to convert them to FASTQ due to incorrect file path

```python

# WRONG PATH (what I initially wrote):
sra_path = os.path.expanduser(f"~/ncbi/public/sra/{sra_id}.sra")
# This looks for: ~/ncbi/public/sra/SRR7179504.sra ❌

```

But `prefetch` actually saves files to:

```python
# CORRECT PATH (what it should be):
sra_path = os.path.expanduser(f"~/ncbi/public/sra/{sra_id}/{sra_id}.sra")
# This finds: ~/ncbi/public/sra/SRR7179504/SRR7179504.sra ✅
```

**Visual comparison:**
```
WRONG:   ~/ncbi/public/sra/SRR7179504.sra  ❌ File doesn't exist here
CORRECT: ~/ncbi/public/sra/SRR7179504/SRR7179504.sra  ✅ Actual location
                          ^^^^^^^^^
                          Missing this folder level!
```

So when we use just `download_fastq.py` we get this 

![Screenshot 2025-12-15 at 19.30.55.png](README/Screenshot_2025-12-15_at_19.30.55.png)

### Results

**After `download_fastq.py`:**

- ✅ Created folders with `.sra` files: `~/ncbi/public/sra/SRR7179504/`, `SRR7179505/`, etc.
- ❌ No FASTQ files created (conversion failed)

where as when we execute this `download_fastq_fixed.py`

- ✅ Successfully converted all `.sra` files
- ✅ Output: 40 compressed FASTQ files in `fastq/` directory (paired-end: `_1.fastq.gz` and `_2.fastq.gz` for each sample)

![Screenshot 2025-12-15 at 17.23.22.png](README/Screenshot_2025-12-15_at_17.23.22.png)

### Key Lesson

When using `prefetch` + `fastq-dump`, remember that `prefetch` creates a **subfolder** for each accession. Always include this extra directory level in your path: `{accession}/{accession}.sra`

---

### Step 2: Quality Control

**Tools:** FastQC, MultiQC

**Installation:**

```bash
conda install -c bioconda fastqc multiqc
```

**Execution:**

```bash
# Create output directory
mkdir -p fastqc_results
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8

#For Multi QC
multiqc fastqc_results/ -o multiqc_report/   
```

**Metrics Evaluated:**
- Per-base sequence quality-Shows the quality scores across all base positions in reads, helping identify if quality degrades toward read ends.
- Per-sequence quality scores- Displays the distribution of average quality scores across all reads to identify low-quality read subsets.
- GC content distribution-  Compares the observed GC content distribution against the expected theoretical distribution to detect contamination or bias.
- Sequence duplication levels- Measures the degree of read duplication, which may indicate PCR over-amplification or low library complexity.
- Adapter contamination-Detects the presence of sequencing adapter sequences that should be trimmed before downstream analysis.
- Overrepresented sequences-Identifies sequences appearing unusually frequently, which may indicate contamination, adapters, or highly abundant transcripts.

**Key Learning:** Some FastQC “warnings” (e.g., per-base sequence content, duplication levels) are **normal for RNA-seq data** due to the nature of gene expression (some genes are highly expressed).

**Output:** HTML reports showing comprehensive quality metrics

---

### Step 3: Read Trimming (Optional)

**Tool:** Trimmomatic

**Installation:** Ensure Java is installed first 

```bash
# Install Java 

echo 'export PATH="/opt/homebrew/opt/openjdk/bin:$PATH"' >> ~/.bash_profile
source ~/.bash_profile

# Test it
java -version

```

```bash
#Move to directory where fastq is 
cd /Users/aakash/Downloads/rnaseq
java -jar ~/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 fastq/SRR7179504_pass.fastq.gz fastq/SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33
```

**Single Sample Example:**

Parameters used:

- LEADING:3 - Remove leading bases with quality < 3- The first few bases of reads often have lower quality due to sequencing chemistry issues during the initial cycles. Quality score 3 means ~50% error rate - these bases are essentially random noise. Removing them prevents incorrect alignments and false variant calls at read starts.
- TRAILING:10 - Remove trailing bases with quality < 10.Sequencing quality typically degrades toward the end of reads as the sequencing chemistry becomes exhausted. Quality score 10 means ~10% error rate. End-of-read errors can cause misalignments and reduce mapping rates. Trimming them improves alignment accuracy.
- SLIDINGWINDOW:4:15 - Cut when average quality in 4-base window drops below 15Cut when average quality in 4-base window drops below 15.his catches quality drops in the **middle** of reads (not just ends). It scans a 4-base window and cuts when average quality falls below 15 (~3% error rate). Preserves good-quality portions while removing problematic regions
- MINLEN:36 - Drop reads shorter than 36 bases after trimming ery short reads (after trimming) are problematic: Increase false positives in analysis and  Can align to multiple genome locations (non-unique mapping). 36 bases is typically long enough for unique alignment to the transcriptome while allowing some trimming flexibility. For reference genomes with low complexity or repetitive regions, you might need even longer (50-75bp).

**Decision:** After testing, trimming was **not used** for final analysis because:
- Initial quality scores were already high
- HISAT2 aligner handles low-quality bases effectively
- Minimal improvement observed in post-trimming QC

---

### Step 4: Sample Merging

**Purpose:** Combine technical replicates (multiple sequencing runs of the same biological sample)

What’s  happening here ? 

1. Concatenating Multiple Sequencing Runs (LNCaP samples-The LNCaP samples were sequenced across 4 separate **runs** that need to be merged into single files.
2. Renaming- The PC3 samples were single runs that just need meaningful names.

```bash
# LNCaP Normoxia Replicate 1cat SRR7179504.fastq.gz SRR7179505.fastq.gz \    SRR7179506.fastq.gz SRR7179507.fastq.gz > LNCAP_Normoxia_S1.fastq.gz
# LNCaP Normoxia Replicate 2cat SRR7179508.fastq.gz SRR7179509.fastq.gz \    SRR7179510.fastq.gz SRR7179511.fastq.gz > LNCAP_Normoxia_S2.fastq.gz
# LNCaP Hypoxia Replicate 1cat SRR7179520.fastq.gz SRR7179521.fastq.gz \    SRR7179522.fastq.gz SRR7179523.fastq.gz > LNCAP_Hypoxia_S1.fastq.gz
# LNCaP Hypoxia Replicate 2cat SRR7179524.fastq.gz SRR7179525.fastq.gz \    SRR7179526.fastq.gz SRR7179527.fastq.gz > LNCAP_Hypoxia_S2.fastq.gz
# PC3 samples (single runs, just renamed)mv SRR7179536.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541.fastq.gz PC3_Hypoxia_S2.fastq.gz
```

**Final Sample Structure:**

| Sample Name | Cell Line | Condition | Replicate | Technical Runs Merged |
| --- | --- | --- | --- | --- |
| LNCAP_Normoxia_S1 | LNCaP | Normoxia | 1 | 4 |
| LNCAP_Normoxia_S2 | LNCaP | Normoxia | 2 | 4 |
| LNCAP_Hypoxia_S1 | LNCaP | Hypoxia | 1 | 4 |
| LNCAP_Hypoxia_S2 | LNCaP | Hypoxia | 2 | 4 |
| PC3_Normoxia_S1 | PC3 | Normoxia | 1 | 1 |
| PC3_Normoxia_S2 | PC3 | Normoxia | 2 | 1 |
| PC3_Hypoxia_S1 | PC3 | Hypoxia | 1 | 1 |
| PC3_Hypoxia_S2 | PC3 | Hypoxia | 2 | 1 |

---

### Step 5: Reference Genome & Annotation

Hisat2 is chosen over  STAR Simply because of practicality issues- I'm working on my personal computer and it would be easier for me to run this simply because hisat2 is memory efficient . [HISAT2 typically requires less RAM for genome indexing and alignment compared to STAR](https://biologyinsights.com/hisat2-vs-star-which-aligner-to-use-for-rna-seq/).

The human reference genome GRCh38 (Homo sapiens) was used. Prebuilt HISAT2 genome indices were downloaded to avoid local genome indexing and reduce memory usage. Gene annotation was obtained from Ensembl release 115 in GTF format.

**Reference:** Human genome GRCh38 (Homo sapiens)

**Download HISAT2 Prebuilt Index:**

```bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
```

**Download Gene Annotation (Ensembl GTF):**

```bash
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz
```

**Install Alignment Tools:**

```bash
# HISAT2 for RNA-seq alignmentconda
 install -c bioconda hisat2
# SAMtools for BAM manipulationconda 
install -c bioconda samtools
```

---

### Step 6: Genome Alignment

**Tool:** HISAT2 (splice-aware aligner for RNA-seq)

**Single Sample Example:**

```bash
hisat2 -q -x grch38/genome -U fastq/LNCAP_Hypoxia_S1.fastq.gz | \  samtools sort -o bam_files/LNCAP_Hypoxia_S1.bam
samtools index bam_files/LNCAP_Hypoxia_S1.bam
```

**Automated Alignment:**
Used bash script: `align_samples.sh`

**Key Alignment Metrics:**
- Overall alignment rate: ~92%
- Uniquely aligned reads: ~84%
- Multiple alignments: ~8%
- Unaligned reads: ~8%

**Output:** 8 sorted and indexed BAM files 

![Screenshot 2025-12-16 at 09.31.23.png](README/Screenshot_2025-12-16_at_09.31.23.png)

---

### Step 7: Post-Alignment Quality Control

**Tool:** Qualimap RNA-seq

**Installation:**

```bash
conda install bioconda::qualimap
```

**Execution:**

```bash
qualimap rnaseq \  -bam bam_files/LNCAP_Hypoxia_S1.bam \  -gtf Homo_sapiens.GRCh38.115.gtf \  -outdir Qualimap_results/LNCAP_Hypoxia_S1 \  --java-mem-size=8G
```

**Quality Metrics Assessed:**
- Mapping quality distribution
- Gene body coverage (5’ to 3’)
- Junction analysis
- Library strandedness
- Genomic origin of reads (exonic/intronic/intergenic)

**Results:**
- LNCaP samples: 67-68% reads assigned to genes ✓
- PC3 samples: 22-27% reads assigned to genes (concerning but usable)

When you load a BAM along with BAI file in IGV, you can visualize:

1. **Individual reads** aligned to the genome
2. **Coverage/depth** (how many reads cover each position)
3. **Splice junctions** ( shows exon-exon connections)
4. **Gene annotations** (if you load a GTF/GFF file alongside)
5. **Variants/mismatches** (differences between reads and reference)

---

### Step 8: Gene Quantification

**Tool:** featureCounts (Subread package)

**Installation:**

```bash
conda install -c bioconda subread
```

**Single Sample Example:**

```bash
featureCounts -S 2 -a Homo_sapiens.GRCh38.115.gtf \  -o quants/LNCAP_Hypoxia_S1_counts.txt \  bam_files/LNCAP_Hypoxia_S1.bam
```

**Parameters:**
- `-S 2`: Reverse-stranded library
- `-a`: GTF annotation file
- `-o`: Output file

**Automated Quantification:**
Used bash script: `featurecounts.sh`

**Output:** Individual count files for each sample

**Count Matrix Generation:**
Used Python script (`merge_counts.py`) to combine all individual featureCounts outputs into a single matrix:
- Rows: ~60,000 genes (Ensembl IDs)
- Columns: 8 samples
- Values: Raw read counts

**Final Output:** `counts_matrix.csv` (ready for DESeq2)

---

### 

## Analysis Pipeline

### 1. Data Preprocessing

- **Count Matrix Loading**: Raw gene counts were imported from `counts_matrix.csv`
- **Metadata Creation**: Sample conditions were assigned and organized
- **Gene Annotation**: Ensembl gene IDs were mapped to gene symbols and biotypes using GRCh38 annotation

### 2. Quality Control and Filtering

### Gene Biotype Filtering

Genes were retained based on functional relevance:

- `protein_coding`
- Immunoglobulin genes: `IG_J_gene`, `IG_V_gene`, `IG_C_gene`, `IG_D_gene`
- T-cell receptor genes: `TR_D_gene`, `TR_C_gene`, `TR_V_gene`, `TR_J_gene`

### Low Expression Filtering

- **Criterion**: Genes with zero counts in ≥7 out of 8 samples were removed
- **Rationale**: Reduces noise and improves statistical power by focusing on reliably detected genes

**Results:**

- Starting genes: 78,899
- After biotype filtering: ~20,000 genes
- After zero-count filtering: ~18,000 genes

![genebiotype_proportions1.png](README/genebiotype_proportions1.png)

Bar chart showing the proportion of different gene biotypes retained after quality control filtering. Protein-coding genes dominate the dataset (>99%), while immunoglobulin (IG) and T-cell receptor (TR) gene variants represent <1% of the total. The final filtered dataset contains approximately 18,000 genes suitable for differential expression analysis.

### 3. Data Normalization and Transformation

- **Method**: Variance Stabilizing Transformation (VST) from DESeq2
- **Purpose**: Stabilizes variance across the dynamic range of gene expression for downstream analyses

### 4. Exploratory Data Analysis

### Principal Component Analysis (PCA)

PCA was performed to assess:

- Overall sample similarity and clustering patterns
- Separation between experimental conditions
- Potential batch effects or outliers

![pcab.png](README/pcab.png)

**Principal component analysis of RNA-seq samples colored by experimental condition.** PCA plot based on variance-stabilized gene expression data shows clear separation of samples along two principal components. PC1 (98% variance) primarily separates samples by cell line identity, with LNCaP cells (red and green) clustering on the left and PC3 cells (cyan and purple) on the right. PC2 (1% variance) separates samples by oxygen condition, with hypoxia samples (red and cyan) at the bottom and normoxia samples (green and purple) at the top. Biological replicates cluster tightly together, demonstrating high experimental reproducibility. The dominant PC1 variance reflects inherent transcriptional differences between androgen-sensitive (LNCaP) and androgen-independent (PC3) prostate cancer cell lines, while PC2 captures the hypoxia-induced transcriptional response.

### 5. Differential Expression Analysis

### Method: DESeq2

DESeq2 performs:

1. Size factor normalization
2. Dispersion estimation (gene-wise and fitted)
3. Negative binomial generalized linear model (GLM) fitting
4. Wald test for significance

### Comparisons Performed

**Primary Comparison: LNCaP Hypoxia vs. Normoxia**

- **Threshold**: |log2FC| > 1, padj < 0.05
- **Significant genes**: Identified differentially expressed genes (DEGs) responding to hypoxia in LNCaP cells

**Secondary Comparison: PC3 Hypoxia vs. Normoxia**

- Similar thresholds applied to identify PC3-specific hypoxia response

![vp_lncap.png](README/vp_lncap.png)

**Volcano plot of differentially expressed genes in response to hypoxia.** Each point represents a single gene plotted by log2 fold change (x-axis) and statistical significance as -log10 adjusted p-value (y-axis). Significantly upregulated genes (orange, right side) show increased expression under hypoxia, while downregulated genes (purple, left side) show decreased expression, both meeting thresholds of |log2FC| > 1 and padj < 0.05 (horizontal dashed line). Non-significant genes are shown in gray. The distribution reveals a robust transcriptional response to hypoxia, with numerous genes showing strong upregulation (several exceeding 8-fold change) and moderate downregulation. The asymmetric distribution suggests hypoxia predominantly activates rather than suppresses gene expression in this cell line.

### 6. Functional Enrichment Analysis

### Gene Set Enrichment Analysis (GSEA)

- **Database**: Reactome pathways
- **Method**: Pre-ranked GSEA using log2 fold changes
- **Purpose**: Identify coordinated changes in biological pathways rather than individual genes

**Key Metrics:**

- **NES (Normalized Enrichment Score)**: Magnitude and direction of pathway enrichment
    - Positive NES: Pathway upregulated in hypoxia
    - Negative NES: Pathway downregulated in hypoxia
- **FDR (False Discovery Rate)**: Statistical significance (padj < 0.05 considered significant)

### Hallmark Gene Set Analysis (fgsea)

- **Database**: MSigDB Hallmark gene sets (50 well-defined biological processes)
- **Method**: Fast GSEA implementation
- **Visualization**: Waterfall plot showing all hallmark pathways ranked by NES

See `results/hallmark_waterfall_plot.png` for comprehensive pathway overview.

![Screenshot 2025-12-18 at 18.04.23.png](README/Screenshot_2025-12-18_at_18.04.23.png)

**Global pathway analysis of hypoxia-induced transcriptional changes.** GSEA waterfall plot ranking all Hallmark pathways by normalized enrichment score. Positively enriched pathways (right, cyan bars) include key hypoxia-adaptive responses: direct hypoxia response, angiogenesis, glycolysis, and epithelial-mesenchymal transition. Negatively enriched pathways (left, salmon bars) predominantly involve oxygen-dependent processes including oxidative phosphorylation, fatty acid metabolism, and MYC targets, indicating suppression of energy-intensive biosynthetic programs. Statistically significant pathways (padj < 0.05, darker colors) demonstrate that hypoxia triggers a coordinated cellular response involving metabolic reprogramming, altered signaling cascades, and adaptive stress responses characteristic of tumor cell adaptation to low oxygen environments.

### 7. Pathway-Specific Visualization

Individual enrichment plots were generated for key pathways:

- **HALLMARK_GLYCOLYSIS**: Expected to be upregulated under hypoxia (metabolic adaptation)
- **HALLMARK_OXIDATIVE_PHOSPHORYLATION**: Expected to be downregulated (oxygen-dependent metabolism)

These plots show:

- Running enrichment score (ES) across ranked genes
- Gene positions within the pathway (barcode-like marks)
- Statistical significance and enrichment magnitude

![Screenshot 2025-12-18 at 18.08.17.png](README/Screenshot_2025-12-18_at_18.08.17.png)

### 8. Gene-Level Expression Analysis

Custom function `plot_counts()` visualizes normalized expression of individual genes across conditions:

- Box plots showing distribution across samples
- Individual data points for each replicate
- Useful for validating DEGs and exploring candidate genes

---

## 📊 Key Results Summary

### Alignment & Quantification Metrics

| Metric | Value |
| --- | --- |
| Total samples processed | 8 |
| Total raw reads | ~108 million |
| Average alignment rate | 92% |
| Genes quantified | ~60,000 |
| Genes after filtering | ~20,000 |

### Differential Expression (LNCaP Hypoxia vs Normoxia)

| Category | Count |
| --- | --- |
| Total DEGs (padj < 0.05) | ~2,500 |
| Upregulated (log2FC > 1) | ~1,200 |
| Downregulated (log2FC < -1) | ~1,300 |
| Highly significant (padj < 0.001) | ~1,800 |

### Biological Interpretation

**Hypoxia activates:**
1. **HIF-1α pathway** - Master regulator of hypoxia response
2. **Glycolysis** - Metabolic switch from oxidative to glycolytic metabolism
3. **Angiogenesis** - VEGF signaling to promote blood vessel formation
4. **Survival pathways** - Anti-apoptotic gene expression
5. **pH regulation** - Carbonic anhydrase genes (CA9, CA12)

**Hypoxia suppresses:**
1. **Oxidative phosphorylation** - Mitochondrial respiration genes
2. **Cell proliferation** - Cell cycle genes
3. **Differentiation markers** - Tissue-specific genes

**Biological Significance:**
- Results align perfectly with established hypoxia biology
- Confirms LNCaP cells mount a robust transcriptional response to hypoxia
- Identified potential therapeutic targets (e.g., VEGFA, CA9, LDHA)
- Provides biomarkers for hypoxic tumor regions

---

## 📁 Repository Structure

### 

---

## 🛠️ Software & Tools

### Command Line Tools

- **SRA Toolkit** (3.2.1) - Data download from NCBI
- **FastQC** (0.12.1) - Quality control of raw reads
- **MultiQC** (1.14) - Aggregate QC reports
- **Trimmomatic** (0.39) - Read trimming (tested but not used)
- **HISAT2** (2.2.1) - Splice-aware genome alignment
- **SAMtools** (1.17) - BAM file manipulation
- **featureCounts** (Subread 2.0.6) - Gene quantification
- **Qualimap** (2.3) - Post-alignment quality control

### R/Bioconductor Packages

- **DESeq2** (1.42.0) - Differential expression analysis
- **clusterProfiler** (4.10.0) - ID conversion and enrichment
- **ReactomePA** (1.46.0) - Reactome pathway analysis
- **fgsea** (1.28.0) - Fast GSEA implementation
- **org.Hs.eg.db** (3.18.0) - Human gene annotation
- **ggplot2** (3.4.4) - Data visualization
- **pheatmap** (1.0.12) - Heatmap generation
- **ggrepel** (0.9.4) - Label repulsion for plots

### Python Libraries

- **Pandas** - Data manipulation for count matrix merging
- **BeautifulSoup** - FastQC HTML report parsing

---

### Biological Insights Gained

1. **Hypoxia Biology**
    - Hypoxia triggers profound transcriptional reprogramming
    - HIF-1α is a master regulator of the hypoxia response
    - Metabolic shift from oxidative to glycolytic metabolism
    - Activation of pro-survival and pro-angiogenic pathways
2. **Cancer Biology**
    - Hypoxia promotes tumor aggression and metastasis
    - Identified potential therapeutic targets (VEGFA, CA9, LDHA)
    - Understanding of prostate cancer cell line differences
3. **Systems Biology**
    - Pathway-level analysis provides biological context
    - Individual genes work as coordinated networks
    - Importance of functional enrichment for interpretation

### 

---

## 🚀 Future Directions

### Potential Extensions

1. Complete the rest of dataset
2. **Advanced Pathway Analysis**
    - Gene Ontology (GO) enrichment
    - KEGG pathway analysis
    - Network analysis of gene interactions
    - Transcription factor enrichment
3. **Additional Visualizations**
    - Interactive plots with plotly
    - Gene set enrichment plots (enrichment scores)
    - Compare with other hypoxia datasets
4. **Integration with Other Data**
    - ChIP-seq for HIF-1α binding sites
    - Proteomics data for post-transcriptional regulation
    

---

## 📚 References & Resources

### Dataset & Publication

- **Original Study:** Guo et al. (2019). “Inhibiting histone deacetylases suppresses glucose metabolism and hepatocellular carcinoma growth by restoring FBP1 expression.” *Nature Communications* 10:1622
- **GEO Accession:** GSE106305
- **SRA Project:** SRP133473

### Tutorial & Learning Resources

- **Eric Liu RNA-seq Tutorial:**   [https://github.com/erilu/bulk-rnaseq-analysis?tab=readme-ov-file#obtaining-raw-data-from-geo](https://github.com/erilu/bulk-rnaseq-analysis?tab=readme-ov-file#obtaining-raw-data-from-geo)
- **DESeq2 Vignette:** http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- **HISAT2 Manual:** http://daehwankimlab.github.io/hisat2/manual/
- **featureCounts Documentation:** [https://subread.sourceforge.net/SubreadUsersGuide.pdf](https://subread.sourceforge.net/SubreadUsersGuide.pdf)

---

## 📬 Contact & Collaboration

**Author:** Aakash Deva Thirukonda Prakash

**Purpose:** Educational project demonstrating proficiency in bioinformatics and computational biology

**Skills Demonstrated:**
- RNA-seq data analysis
- Statistical programming in R
- Bash scripting and automation
- Biological data interpretation
- Scientific visualization
- Reproducible research practices

For questions, suggestions, or potential collaborations, please feel free to reach out through GitHub.

---

## 📄 License

This project is created for educational purposes as part of learning bioinformatics analysis workflows. The data used is publicly available from NCBI GEO (GSE106305).

---

**Future Work:**
- ⏳ PC3 differential expression analysis
- ⏳ Comparative analysis between cell lines
- ⏳ Additional pathway databases
- ⏳ Interactive visualizations

---
