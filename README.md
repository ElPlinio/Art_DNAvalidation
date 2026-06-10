# Microbiome Analysis of Myxomycete Plasmodia

R scripts and QIIME 2 pipeline for the analyses and figures in:

> Molina-Viramontes JP, Gómez-Acata ES, Hereira-Pacheco S, Estrada-Torres A, Navarro-Noya YE. (2025). Comparison of Methods to Characterize the Microbiota of Myxomycete Plasmodia. *Journal of Eukaryotic Microbiology*, 73, e70058. https://doi.org/10.1111/jeu.70058

Raw 16S sequencing data: [NCBI BioProject PRJNA1281717](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1281717)

See [`README_bioinformatics.md`](README_bioinformatics.md) for the full QIIME 2 processing pipeline used to generate the input data tables.

---

## Repository structure

```
.
├── Code/
│   ├── taxonomic_diversity.R   # Bar plots (Phylum/Family/Genus) + alpha diversity boxplots
│   ├── sankey_plot.R           # Sankey taxonomic composition diagram
│   └── heatmap_plot.R          # ASV heatmap with phylogenetic ordering
├── Data/                       # Input data files (see data/README.md)
```

---

## Scripts

### `taxonomic_diversity.R`
Generates stacked bar plots of relative abundance at Phylum, Family, and Genus level, and paired boxplots of alpha diversity (Hill numbers q0, q1, q2) comparing commercial kit vs. non-commercial extraction methods.

**Key packages:** `qiime2R`, `hilldiv`, `hilldiv2`, `hillR`, `tidyverse`, `ggpubr`, `cowplot`, `ggh4x`

### `sankey_plot.R`
Hierarchical Sankey diagram (Kingdom → Phylum → Class → Genus → Species) showing average bacterial relative abundance across all samples.


**Key packages:** `qiime2R`, `tidyverse`, `sankeyD3`, `plyr`

### `heatmap_plot.R`
Heatmap of ASV relative abundances (12 discretised bins) with rows ordered by phylogenetic tree, annotated by Phylum and prevalence in each extraction method.


**Key packages:** `ComplexHeatmap`, `circlize`, `ape`, `qiime2R`, `tidyverse`, `viridis`, `RColorBrewer`

---

## Setup


### Install R packages
```r
install.packages(c("tidyverse", "plyr", "scales", "viridis", "RColorBrewer",
                   "reshape2", "ggh4x", "cowplot", "vegan", "phytools", "ape", "circlize"))

# Bioconductor
BiocManager::install(c("ComplexHeatmap", "ANCOMBC", "phyloseq"))

# GitHub
remotes::install_github("jgrembi/qiime2R")
remotes::install_github("fbreitwieser/sankeyD3")
remotes::install_github("anttonalberdi/hilldiv2")
install.packages(c("hilldiv", "hillR", "ggpubr"))
```

---

## Authors

- Juan Pablo Molina-Viramontes
- Elizabeth Selene Gómez-Acata
- Stephanie Hereira-Pacheco
- Arturo Estrada-Torres
- Yendi E. Navarro-Noya

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
