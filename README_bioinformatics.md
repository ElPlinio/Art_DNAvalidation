# Bioinformatics Analyses — Detailed Description

This document describes the full bioinformatics pipeline used to process 16S rRNA amplicon sequencing data for:

> Molina-Viramontes JP, Gómez-Acata ES, Hereira-Pacheco S, Estrada-Torres A, Navarro-Noya YE. (2025). Comparison of Methods to Characterize the Microbiota of Myxomycete Plasmodia. *Journal of Eukaryotic Microbiology*, 73, e70058. https://doi.org/10.1111/jeu.70058

Raw 16S sequencing data: [NCBI BioProject PRJNA1281717](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1281717)

---

## 1. Sequencing and raw data

Paired-end 16S rRNA amplicon sequencing targeting the V3–V4 hypervariable region was performed on an Illumina MiSeq (300 bp paired-end).

Samples were divided into two groups:
- **Commercial kit** (Q1–Q6): DNA extracted with DNeasy PowerSoil (Qiagen)
- **Non-commercial** (Fe1–Fe6): DNA extracted using phenol-chloroform protocol
- **Controls**: FeN and Negkit (negative controls) — excluded from downstream analysis

---

## 2. QIIME 2 processing pipeline

All preprocessing was done in **QIIME 2 v2023.4**.

### 2.1 Denoising — DADA2
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired_end_demux.qza \
  --p-trim-left-f [F] --p-trim-left-r [R] \
  --p-trunc-len-f [LF] --p-trunc-len-r [LR] \
  --o-table table.qza \
  --o-representative-sequences rep-seq-255195.qza \
  --o-denoising-stats dada2-stats.qza
```
Output: `table-255195-clean_predoc.qza`

### 2.2 Taxonomic classification — GreenGenes 2
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier gg2-classifier.qza \
  --i-reads rep-seq-255195.qza \
  --o-classification taxgg255195gg.qza

qiime metadata tabulate \
  --m-input-file taxgg255195gg.qza \
  --o-visualization taxgg255195gg.qzv
```
Nomenclature follows GTDB/GreenGenes 2 (e.g., *Pseudomonadota*, *Bacillota*, *Actinomycetota*).

### 2.3 Removal of mitochondrial sequences
Mitochondrial 16S sequences were identified and removed from the feature table.

### 2.4 Decontamination
Negative control samples (blanks) were used to identify contaminant ASVs with `decontam`:

```bash
# Identify contaminant ASVs from blank samples
qiime feature-table filter-samples \
  --i-table table-255195-clean_predoc.qza \
  --m-metadata-file meta.txt \
  --p-where '[edad]="negativo"' \
  --o-filtered-table filtered-table_blank.qza

qiime feature-table summarize \
  --i-table filtered-table_blank.qza \
  --o-visualization summ_filtered_table_blank.qzv

# Remove identified contaminant features
echo "feature-id" > id_decontam.txt
echo "64bf3513a819a26ce7a3be59ac8725e9" >> id_decontam.txt
echo "9c911aa0503e607d716c2f11ad37dec8" >> id_decontam.txt
echo "6003d85eeb90ee5f382be29a1c01124f" >> id_decontam.txt

qiime feature-table filter-features \
  --p-exclude-ids TRUE \
  --i-table table-255195-clean_predoc.qza \
  --m-metadata-file id_decontam.txt \
  --o-filtered-table table-255195-clean_predoc_decontam.qza
```

### 2.5 Rarefaction
```bash
qiime feature-table rarefy \
  --i-table table-255195-clean_predoc_decontam.qza \
  --p-sampling-depth 40625 \
  --o-rarefied-table table-255195-clean_predoc_decontam_rare40625.qza
```

### 2.6 Additional ASV filtering (for heatmap)
For the heatmap analysis, ASVs were further filtered by prevalence (≥2 samples) and unassigned taxa were removed, producing `table_filt_no_unassigned_2_filterasv.qza`.

### 2.7 Phylogenetic tree construction

```bash
# Filter representative sequences to match the rarefied table
qiime feature-table filter-seqs \
  --i-data rep-seq-255195.qza \
  --i-table table-255195-clean_predoc_decontam_rare40625.qza \
  --o-filtered-data rep-seqs-clean_predoc_decontam_rare40625.qza

# Align + mask + build tree (MAFFT + FastTree, one-step)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-clean_predoc_decontam_rare40625.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# Export tree for R
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path tree/
```
Output: `tree/tree.nwk`

---

## 3. Downstream R analyses

All scripts are in `Code/`. Run them from the project root directory.

| Script | Description | Output |
|--------|-------------|--------|
| `taxonomic_diversity.R` | Taxonomic bar plots (Phylum, Family, Genus) and alpha diversity boxplots (Hill numbers q0, q1, q2) | `figures/barplot_phylum.png`, `figures/boxplot_alpha_div.png` |
| `sankey_plot.R` | Hierarchical Sankey diagram (K→P→C→G→S) of mean relative abundances | `Net13.html` |
| `heatmap_plot.R` | ASV heatmap ordered by phylogenetic tree, annotated by extraction method | `heatmap_circlize_all_mod2.pdf` |

---

## 4. Data files

| File | Used in | Description |
|------|---------|-------------|
| `data/table-255195-clean_predoc_rar40625.qza` | `taxonomic_diversity.R` | Rarefied, decontaminated feature table |
| `data/taxgg255195gg/taxonomy.tsv` | `taxonomic_diversity.R` | GreenGenes 2 taxonomy (exported TSV) |
| `data/metadata.txt` | `taxonomic_diversity.R` | Sample metadata |
| `data/ehi_phylum_colors.tsv` | `taxonomic_diversity.R` | Phylum colour palette |
| `data/table-255195-clean_predoc.qza` | `sankey_plot.R` | Pre-decontam feature table |
| `data/taxonomy/taxonomy.tsv` | `sankey_plot.R` | GreenGenes 2 taxonomy (TSV) |
| `data/table_filt_no_unassigned_2_filterasv.qza` | `heatmap_plot.R` | Prevalence-filtered feature table |
| `data/taxgg255195gg/taxonomy.tsv` | `heatmap_plot.R` | GreenGenes 2 taxonomy (TSV) |
| `data/tree.nwk` | `heatmap_plot.R` | Rooted phylogenetic tree |
| `data/new_names.RDS` | `heatmap_plot.R` | Curated species labels |
| `data/metadata.txt` | `heatmap_plot.R` | Sample metadata |

> **Note:** `ehi_phylum_colors.tsv` is a custom colour palette file included in `Data/`; `meta.txt` is not tracked in this repository — add it before running `taxonomic_diversity.R`.

---

## 5. Software versions

| Tool | Version |
|------|---------|
| QIIME 2 | 2023.4 |
| DADA2 | embedded in QIIME 2 |
| FastTree | embedded in QIIME 2 (`align-to-tree-mafft-fasttree`) |
| R | ≥ 4.2 |
| hilldiv / hilldiv2 | CRAN |
| hillR | CRAN |
| ComplexHeatmap | ≥ 2.14 (Bioconductor) |
| ANCOMBC | Bioconductor |
| phyloseq | Bioconductor |
| circlize, ape, vegan, phytools | CRAN |
| tidyverse, ggpubr, cowplot, ggh4x | CRAN |
| qiime2R | GitHub (jgrembi/qiime2R) |
| sankeyD3 | GitHub (fbreitwieser/sankeyD3) |

---

*For questions about the pipeline, contact the corresponding author (see article).*
