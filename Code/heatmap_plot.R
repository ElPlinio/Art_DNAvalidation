# ============================================================
# Heatmap — ASV relative abundance by DNA extraction method
# Article: Juan Pablo et al. (see README)
# Script: complexheatmap-circlize_all_mod.R
# ============================================================
# Produces a ComplexHeatmap with:
#   - Rows: ASVs (species-level labels), ordered by phylogenetic tree
#   - Columns: samples (Fe = non-commercial; Q = commercial kit)
#   - Left annotations: Phylum, Commercial kit prevalence,
#                       Non-commercial prevalence
#   - Color scale: discretised relative abundance (viridis "B")
# Output: heatmap_circlize_all_mod2.pdf
# ============================================================

library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(scales)
library(qiime2R)
library(ape)
library(RColorBrewer)

# --- Helper functions ---
relabunda   <- function(x) { (as.data.frame(t(t(x) / colSums(x))) * 100) }
prevalence  <- function(x) { x %>% mutate_if(is.numeric, ~ifelse(. > 0, 1, 0)) }

# --- Load data ---
otu  <- read_qza("data/table_filt_no_unassigned_2_filterasv.qza")$data
tax  <- read.delim("data/taxgg255195gg/taxonomy.tsv")
map  <- read.delim("data/metadata.txt") %>% dplyr::select(-Individuo) %>% rename(Individo = X)
tree <- ape::read.tree("data/tree.nwk")

# --- Taxonomy from tree tip labels ---
names_ASV <- tree$tip.label %>% as.data.frame() %>%
  rename("Feature.ID" = ".") %>%
  inner_join(tax) %>%
  mutate(
    Specie = str_extract(Taxon, "[^_]+$"),
    Phylum = str_extract(Taxon, "(?<=p__)\\w+"),
    Class  = str_extract(Taxon, "(?<=c__)\\w+"),
    Phylum = case_when(
      Class  == "Alphaproteobacteria"  ~ "Pseudomonadota",
      Class  == "Gammaproteobacteria"  ~ "Pseudomonadota",
      Phylum == "Actinobacteriota"     ~ "Actinomycetota",
      Phylum == "Firmicutes_D"         ~ "Bacillota_D",
      .default = as.character(Phylum)
    )
  )

# --- Abundance and prevalence ---
abundance       <- otu %>% as.data.frame() %>% relabunda()
preval_kit      <- otu %>% as.data.frame() %>% prevalence() %>%
  dplyr::select(starts_with("Q")) %>%
  rowSums() %>% as.data.frame() %>% prevalence() %>%
  dplyr::rename("Kit" = ".")
preval_fenol    <- otu %>% as.data.frame() %>% prevalence() %>%
  dplyr::select(starts_with("Fe")) %>%
  rowSums() %>% as.data.frame() %>% prevalence() %>%
  rename("Fenol" = ".")

wt.master <- cbind(abundance, preval_kit, preval_fenol)

heat <- wt.master %>%
  rownames_to_column(var = "Feature.ID") %>%
  inner_join(names_ASV)

# Load curated species names (manually edited)
new_names    <- read_rds("data/new_names.RDS")
heat$Specie  <- new_names

# --- Build discretised heatmap matrix (Fe1:Q6 columns) ---
heatm <- heat %>%
  column_to_rownames(var = "Specie") %>%
  dplyr::select(Fe1:Q6) %>%
  as.data.frame() %>%
  mutate_all(., funs(data = case_when(
    .              <= 0.001  ~  0,
    . >  0.001 & . <= 0.005  ~  1,
    . >  0.005 & . <= 0.01   ~  2,
    . >  0.01  & . <= 0.10   ~  3,
    . >  0.10  & . <= 0.20   ~  4,
    . >  0.20  & . <= 1.00   ~  5,
    . >  1.00  & . <= 2.00   ~  6,
    . >  2.00  & . <= 5.00   ~  7,
    . >  5.00  & . <= 10.00  ~  8,
    . >  10.00 & . <= 25.00  ~  9,
    . >  25.00 & . <= 50.00  ~ 10,
    . >  50.00 & . <= 75.00  ~ 11,
    . >  75.00               ~ 12
  ))) %>%
  select(ends_with("_data")) %>%
  rename_all(~gsub("_data", "", .))

my_palette <- viridis::viridis(n = 12, option = "B", direction = -1)

# --- Phylogenetic tree → dendrogram for row ordering ---
mytree           <- ape::read.tree("data/tree.nwk")
taxo             <- mytree$tip.label %>% as.data.frame() %>% rename(Feature.ID = ".")
taxa             <- taxo %>% left_join(heat) %>%
  mutate_if(is.character, ~str_extract(., "[^_]+$"))

mytree$tip.label <- taxa$Specie
ultrametric_tree <- chronopl(mytree, lambda = 0)
hc               <- as.hclust(ultrametric_tree)
newTree          <- as.dendrogram(hc)

# Reorder heatmap rows to match tree
heatms <- heatm %>%
  rownames_to_column(var = "ID") %>%
  { .[match(mytree$tip.label, .$ID), ] } %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(var = "ID")

# --- Row annotations ---
annotation_columns <- heat %>%
  inner_join(wt.master) %>%
  dplyr::select(Specie, Phylum,
                "Commercial kit" = Kit,
                "Non-commercial" = Fenol) %>%
  column_to_rownames(var = "Specie")

annotation_columns2 <- annotation_columns %>% rownames_to_column(var = "ID")
annotation_columns  <- annotation_columns2[match(mytree$tip.label, annotation_columns2$ID), ] %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(var = "ID")

# --- Colour palettes ---
cols_phyl <- list('Phylum' = c(
  "Actinomycetota"  = "#D72542",
  "Bacillota_A"     = "#04A8BE",
  "Bacillota_D"     = "#F3C50D",
  "Bacteroidota"    = "#EF5085",
  "Patescibacteria" = "#D96C29",
  "Pseudomonadota"  = "#2D3FA6"
))

annphyl <- HeatmapAnnotation(
  "Phylum"              = annotation_columns$Phylum,
  which                 = "row",
  show_legend           = FALSE,
  annotation_name_gp    = gpar(fontsize = 11, fontface = "bold"),
  gp                    = gpar(col = "white"),
  show_annotation_name  = TRUE,
  col                   = cols_phyl
)

annkit <- HeatmapAnnotation(
  "Commercial kit"      = annotation_columns$`Commercial kit`,
  which                 = "row",
  show_legend           = FALSE,
  annotation_name_gp    = gpar(fontsize = 11, fontface = "bold"),
  gp                    = gpar(col = "white"),
  show_annotation_name  = TRUE,
  col                   = list("Commercial kit" = colorRamp2(c(0, 1), c("white", "#E69F00")))
)

annfenol <- HeatmapAnnotation(
  "Non-commercial"      = annotation_columns$`Non-commercial`,
  which                 = "row",
  show_legend           = FALSE,
  annotation_name_gp    = gpar(fontsize = 11, fontface = "bold"),
  gp                    = gpar(col = "white"),
  show_annotation_name  = TRUE,
  col                   = list("Non-commercial" = colorRamp2(c(0, 1), c("white", "#44803F")))
)

# --- Final heatmap ---
heats <- Heatmap(
  heatms,
  col              = my_palette,
  cluster_rows     = newTree,
  width            = unit(6, "cm"),
  heatmap_legend_param = list(
    direction      = "vertical",
    labels_gp      = gpar(fontsize = 7),
    legend_gp      = gpar(fontsize = 9),
    title          = "Relab(%)",
    title_position = "topcenter",
    break_dist     = 1
  ),
  rect_gp          = gpar(col = "gray", lwd = 1),
  row_names_gp     = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp  = gpar(fontsize = 12),
  cluster_columns  = FALSE,
  show_column_names = TRUE,
  show_heatmap_legend = TRUE,
  left_annotation  = c(annphyl, annkit, annfenol)
)

# --- Save output ---
pdf("heatmap_circlize_all_mod2.pdf", width = 8.5, height = 8)
print(heats)
dev.off()
