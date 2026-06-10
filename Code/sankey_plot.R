# ============================================================
# Sankey Plot — Taxonomic composition across ranks
# Article: Juan Pablo et al. (see README)
# Script: sankey_gg_nofo.R
# ============================================================
# Generates a hierarchical Sankey diagram (K → P → C → G → S)
# showing average relative abundance across all samples,
# saved as an interactive HTML file (Net13.html).
# ============================================================

library(qiime2R)
library(tidyverse)
library(sankeyD3)

# --- Load data ---
otu  <- read_qza("data/table-255195-clean_predoc.qza")$data %>% as.data.frame()
otu  <- otu %>% dplyr::select(!c(FeN, Negkit))
taxa <- read.delim("data/taxonomy/taxonomy.tsv")

# --- Relative abundance function ---
relabunda <- function(x) { (as.data.frame(t(t(x) / colSums(x))) * 100) }

otu_kit <- otu %>% relabunda()

otu_kit_parse <- otu_kit %>%
  rownames_to_column(var = "Feature.ID") %>%
  inner_join(taxa) %>%
  rownames_to_column(var = "ids") %>%
  separate(Taxon, c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate_if(is.character, ~str_extract(., "[^_]+$")) %>%
  dplyr::select(-Confidence)

# --- Summarise by taxonomic rank ---

# Kingdom
bacterias <- otu_kit_parse %>%
  group_by(k) %>% summarise_if(is.numeric, sum) %>%
  column_to_rownames(var = "k") %>% t() %>% as.data.frame() %>%
  { colMeans(.) } %>% as.data.frame() %>% dplyr::rename(abund = ".") %>%
  mutate(taxRank = "K") %>% filter(round(.$abund, digits = 3) > 0)

# Phylum
phylum <- otu_kit_parse %>%
  unite("phylum", c("k", "p")) %>% dplyr::select(-c:-s) %>%
  group_by(phylum) %>% summarise_if(is.numeric, sum) %>%
  column_to_rownames(var = "phylum") %>% t() %>% as.data.frame() %>%
  { colMeans(.) } %>% as.data.frame() %>% dplyr::rename(abund = ".") %>%
  mutate(taxRank = "P") %>%
  rownames_to_column(var = "ids") %>%
  filter(!endsWith(ids, "NA")) %>%
  column_to_rownames(var = "ids") %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    Taxon = gsub("Firmicutes",       "Bacillota-D",     Taxon),
    Taxon = gsub("Proteobacteria",   "Pseudomonadota",  Taxon),
    Taxon = gsub("Actinobacteriota", "Actinomycetota",  Taxon),
    Taxon = gsub("Cyanobacteria",    "Cyanobacteriota", Taxon)
  ) %>%
  column_to_rownames(var = "Taxon")

# Class
class <- otu_kit_parse %>%
  unite("class", c("k", "p", "c")) %>% dplyr::select(-o:-s) %>%
  group_by(class) %>% summarise_if(is.numeric, sum) %>%
  column_to_rownames(var = "class") %>% t() %>% as.data.frame() %>%
  { colMeans(.) } %>% as.data.frame() %>% dplyr::rename(abund = ".") %>%
  mutate(taxRank = "C") %>%
  rownames_to_column(var = "ids") %>%
  filter(!endsWith(ids, "NA")) %>%
  column_to_rownames(var = "ids") %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    Taxon = gsub("Firmicutes",       "Bacillota-D",     Taxon),
    Taxon = gsub("Proteobacteria",   "Pseudomonadota",  Taxon),
    Taxon = gsub("Actinobacteriota", "Actinomycetota",  Taxon),
    Taxon = gsub("Cyanobacteria",    "Cyanobacteriota", Taxon)
  ) %>%
  column_to_rownames(var = "Taxon")

# Genus
genus <- otu_kit_parse %>%
  unite("genus", c("k", "p", "c", "o", "f", "g")) %>% dplyr::select(-s) %>%
  group_by(genus) %>% summarise_if(is.numeric, sum) %>%
  column_to_rownames(var = "genus") %>% t() %>% as.data.frame() %>%
  { colMeans(.) } %>% as.data.frame() %>% dplyr::rename(abund = ".") %>%
  mutate(taxRank = "G") %>%
  rownames_to_column(var = "ids") %>%
  filter(!endsWith(ids, "NA")) %>%
  column_to_rownames(var = "ids") %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    Taxon = gsub("Firmicutes",       "Bacillota-D",     Taxon),
    Taxon = gsub("Proteobacteria",   "Pseudomonadota",  Taxon),
    Taxon = gsub("Actinobacteriota", "Actinomycetota",  Taxon),
    Taxon = gsub("Cyanobacteria",    "Cyanobacteriota", Taxon)
  ) %>%
  column_to_rownames(var = "Taxon")

# Species (filtered to key taxa only)
specie <- otu_kit_parse %>%
  unite("specie", c("k", "p", "c", "o", "f", "g", "s")) %>%
  group_by(specie) %>% summarise_if(is.numeric, sum) %>%
  column_to_rownames(var = "specie") %>% t() %>% as.data.frame() %>%
  { colMeans(.) } %>% as.data.frame() %>% dplyr::rename(abund = ".") %>%
  mutate(taxRank = "S") %>%
  rownames_to_column(var = "ids") %>%
  filter(!endsWith(ids, "NA")) %>%
  filter(ids %in% c(
    "Bacteria_Proteobacteria_Gammaproteobacteria_Burkholderiales_Burkholderiaceae_Achromobacter_Achromobacter denitrificans",
    "Bacteria_Proteobacteria_Gammaproteobacteria_Xanthomonadales_Xanthomonadaceae_Stenotrophomonas_Stenotrophomonas nitritireducens",
    "Bacteria_Proteobacteria_lphaproteobacteria_Caulobacterales_Caulobacteraceae_Brevundimonas_Brevundimonas diminuta",
    "Bacteria_Bacteroidota_Bacteroidia_Flavobacteriales_Flavobacteriaceae_Flavobacterium_Flavobacterium injenense",
    "Bacteria_Proteobacteria_lphaproteobacteria_Rhizobiales_Beijerinckiaceae_Chelatococcus_Chelatococcus saccharovorans",
    "Bacteria_Firmicutes_Bacilli_Bacillales_Bacillaceae_Bacillus_Bacillus spizizenii",
    "Bacteria_Proteobacteria_Gammaproteobacteria_595422_595422_Achromobacter_593100",
    "Bacteria_Proteobacteria_Alphaproteobacteria_Caulobacterales_Caulobacteraceae_Brevundimonas_Brevundimonas diminuta",
    "Bacteria_Bacteroidota_Bacteroidia_877923_Flavobacteriaceae_Flavobacterium_Flavobacterium injenense",
    "Bacteria_Actinobacteriota_Thermoleophilia_Solirubrobacterales_405341_Conexibacter_Conexibacter stalactiti"
  )) %>%
  column_to_rownames(var = "ids") %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    Taxon = gsub("Firmicutes",       "Bacillota-D",     Taxon),
    Taxon = gsub("Proteobacteria",   "Pseudomonadota",  Taxon),
    Taxon = gsub("Actinobacteriota", "Actinomycetota",  Taxon),
    Taxon = gsub("Cyanobacteria",    "Cyanobacteriota", Taxon)
  ) %>%
  column_to_rownames(var = "Taxon")

# --- Combine ranks and filter ---
my_report <- rbind(bacterias, phylum, class, genus, specie) %>%
  rownames_to_column(var = "ids") %>%
  mutate(names = str_extract(.$ids, "[^_]+$")) %>%
  dplyr::filter(!names == "NA") %>%
  filter(!ids %in% c(
    "Bacteria_Patescibacteria_Saccharimonadia",
    "Bacteria_Patescibacteria",
    "Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Beijerinckiaceae_Beijerinckiaceae",
    "Bacteria_Actinobacteriota_Thermoleophilia_Solirubrobacterales_Solirubrobacteraceae_Conexibacter_stalactiti",
    "Bacteria_Bacteroidota_Bacteroidia_Flavobacteriales_Flavobacteriaceae_uncultured_Flavobacteriia",
    "Bacteria_Bacteroidota_Bacteroidia_Sphingobacteriales_Sphingobacteriaceae_Sphingobacterium_mizutaii",
    "Bacteria_Bacteroidota_Bacteroidia_Chitinophagales_Chitinophagaceae_Sediminibacterium_Bacteroidetes",
    "Bacteria_Bacteroidota_Bacteroidia_Flavobacteriales_Weeksellaceae_Chryseobacterium_Bacteroidetes",
    "Bacteria_Bacteroidota_Bacteroidia_Flavobacteriales_Flavobacteriaceae_Myroides_injenensis",
    "Bacteria_Firmicutes_Bacilli_Bacillales_Planococcaceae_Sporosarcina_bacterium",
    "Bacteria_Patescibacteria_Saccharimonadia_Saccharimonadales_Saccharimonadaceae_Saccharimonas_bacterium",
    "Bacteria_Proteobacteria_Gammaproteobacteria_Burkholderiales_Neisseriaceae_uncultured_bacterium",
    "Bacteria_Patescibacteria_Saccharimonadia_Saccharimonadales",
    "Bacteria"
  ))

# --- Build Sankey network ---
taxRanks <- c("K", "P", "C", "G", "S")
maxn <- 25

my_report <- subset(my_report, taxRank %in% taxRanks)
my_report <- plyr::ddply(my_report, "taxRank",
                         function(x) x[utils::tail(order(x$abund), n = maxn), , drop = FALSE])
my_report <- my_report[!my_report$name %in% c('-_root'), ]

splits <- strsplit(my_report$ids, "\\_")

root_nodes <- sapply(splits[sapply(splits, length) == 2], function(x) x[2])
sel    <- sapply(splits, length) >= 3
splits <- splits[sel]

links <- data.frame(
  do.call(rbind, lapply(splits, function(x) utils::tail(x[x %in% my_report$name], n = 2))),
  stringsAsFactors = FALSE
)
colnames(links) <- c("source", "target")
links$value <- my_report[sel, "abund"]

my_taxRanks        <- taxRanks[taxRanks %in% my_report$taxRank]
taxRank_to_depth   <- stats::setNames(seq_along(my_taxRanks) - 1, my_taxRanks)

nodes <- data.frame(
  name  = my_report$name,
  depth = taxRank_to_depth[my_report$taxRank],
  value = my_report$abund,
  stringsAsFactors = FALSE
)

names_id     <- stats::setNames(seq_len(nrow(nodes)) - 1, nodes[, 1])
links$source <- names_id[links$source]
links$target <- names_id[links$target]
links        <- links[links$source != links$target, ]

nodes$name        <- sub("^._", "", nodes$name)
links$source_name <- nodes$name[links$source + 1]
links$type        <- sub(' .*', '', nodes[links$source + 1, 'name'])

# --- Final Sankey plot (saved as interactive HTML) ---
sankeyD3::sankeyNetwork(
  Links                    = links,
  Nodes                    = nodes,
  doubleclickTogglesChildren = TRUE,
  LinkGroup                = "type",
  fontFamily               = "Helvetica",
  Source                   = "source",
  Target                   = "target",
  Value                    = "value",
  NodeID                   = "name",
  NodeGroup                = "name",
  NodePosX                 = "depth",
  NodeValue                = "value",
  dragY                    = TRUE,
  xAxisDomain              = my_taxRanks,
  numberFormat             = "pavian",
  title                    = NULL,
  nodeWidth                = 30,
  linkGradient             = TRUE,
  nodeShadow               = TRUE,
  nodeCornerRadius         = 5,
  units                    = "abund",
  fontSize                 = 10,
  nodePadding              = 8,
  width                    = 1200,
  iterations               = 10 * 100,
  align                    = "none",
  highlightChildLinks      = TRUE,
  orderByPath              = TRUE,
  scaleNodeBreadthsByString = TRUE
) %>% saveNetwork(file = 'Net13.html')
