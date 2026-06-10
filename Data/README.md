# Data Files

This folder contains the small metadata and taxonomy files required to run the scripts. Large QIIME 2 artifacts (`.qza`) are not included due to size; they should be placed here following the structure below before running the scripts.

## Files included in this repo

| File | Description |
|------|-------------|
| `metadata.txt` | Sample metadata (individual ID, extraction method, sample group) |
| `taxonomy/taxonomy.tsv` | 16S taxonomic assignments — GreenGenes standard |
| `taxgg255195gg/taxonomy.tsv` | 16S taxonomic assignments — GreenGenes 2 (used in heatmap script) |

## Large files — place manually (not tracked by Git)

These files must be in `PARA_ARTICULO/` relative to the working directory:

| File | Used in script | Description |
|------|---------------|-------------|
| `PARA_ARTICULO/table-255195-clean_predoc.qza` | `sankey_gg_nofo.R` | Raw ASV feature table |
| `PARA_ARTICULO/table_filt_no_unassigned_2_filterasv.qza` | `complexheatmap-circlize_all_mod.R` | Filtered feature table (≥2 samples) |
| `PARA_ARTICULO/export_rooted_seqs_table_filt_no_unassigned_2_filterasv.qza/tree.nwk` | `complexheatmap-circlize_all_mod.R` | Rooted phylogenetic tree (Newick) |

And in the project root:

| File | Used in script | Description |
|------|---------------|-------------|
| `new_names.RDS` | `complexheatmap-circlize_all_mod.R` | Curated species-level label names (RDS object) |

## Notes

- `.qza` files are QIIME 2 artifacts. They can be opened with QIIME 2 or read directly in R with the `qiime2R` package.
- Raw sequencing data are deposited at NCBI SRA under accession **[SRA accession]** (see article).
- If you use Git LFS, `.qza` files can be tracked with: `git lfs track "*.qza"`
