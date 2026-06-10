# ============================================================
# Taxonomic Profile & Alpha Diversity
# Article: Molina-Viramontes et al. (2025). J. Eukaryot. Microbiol., 73, e70058.
#          https://doi.org/10.1111/jeu.70058
# Script: taxonomic_diversity.R
# Authors: Juan Pablo Molina-Viramontes, Yendi E. Navarro-Noya
# ============================================================

# --- Libraries ---
library(qiime2R)
library(hilldiv)
library(hillR)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(hilldiv2)
library(vegan)
library(ape)
library(phytools)
library(reshape2)
library(ggh4x)
library(ANCOMBC)
library(phyloseq)

# --- Colour palettes ---
pal<- c("#D72542", "#04A8BE",  "#F3C50D","#EF5085", "#D96C29",
  "#2D3FA6", "#F2C744", "#699B30", "#F2DEA2", "#4DE3C9" ,
  "#745BB0", "#146152", "#526AB0", "#F21313", "#B09141",
  "#44803F", "#E3906F", "#58B1E3",  "#B4CF66", "#1A4E95")

phylum_col<- read_tsv("data/ehi_phylum_colors.tsv")%>%
  mutate(phylum = gsub("p__", "", phylum))

# --- Load data (paths adjusted: original used ../data/, script runs from project root) ---
tabla<- read_qza("data/table-255195-clean_predoc_rar40625.qza")$data

taxonomy <- read_tsv("data/taxgg255195gg/taxonomy.tsv") %>%
  rename(Feature.ID = `Feature ID`) %>%
  column_to_rownames("Feature.ID") %>%
  mutate(
    Taxon = gsub("Firmicutes",     "Bacillota",       Taxon),
    Taxon = gsub("Proteobacteria", "Pseudomonadota",  Taxon),
    Taxon = gsub("Cyanobacteria",  "Cyanobacteriota", Taxon)
  )

meta<- read_tsv("data/metadata.txt") %>% rename(SampleID="sample-id") %>%
  mutate(metodo = gsub("fenol", "Non-commercial", metodo),
         metodo = gsub("kit", "Commercial kit", metodo))

parse<- parse_taxonomy(taxonomy)

tax<- parse_taxonomy(taxonomy)%>%
  rownames_to_column("Feature.ID")

sum<- summarize_taxa(tabla, parse)

# ============================================================
# 2. TAXONOMIC PROFILE — Bar plots
# ============================================================

p<- sum$Phylum %>% rownames_to_column("phylum") %>%
  mutate(phylum = gsub("d__Bacteria; ", "", phylum)) %>%
  mutate_at(vars(-phylum),~./sum(.)) %>% # relative abundance
  pivot_longer(-phylum, names_to = "SampleID", values_to = "count") %>%
  left_join(., meta, by = join_by(SampleID == SampleID)) %>%
  filter(count > 0) %>%
  mutate(Individuo=factor(Individuo,levels=c("P1","P2","P3","P4","P5","P6"))) %>%
  ggplot(., aes(x=metodo, y=count, fill=phylum, group=phylum)) +
    geom_bar(stat="identity", colour="white", linewidth=0.1) +
    scale_fill_manual(values=pal) +
    facet_nested(. ~ Individuo,  scales="free") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black")) +
   labs(fill="Phylum",y = "Relative abundance",x="Samples")

p
ggsave('figures/barplot_phylum.png',width = 6, height = 4, dpi = 300, plot =p)

sum_fam<- sum$Family %>%
  rownames_to_column() %>%
  filter(!str_detect(rowname, "NA")) %>%
 filter(!str_detect(rowname, "uncultured")) %>%
  column_to_rownames(var="rowname")
taxa <- rownames(sum_fam)
new_taxa <- gsub(".*;", "", taxa)
rownames(sum_fam) <- new_taxa

b<-taxa_barplot(sum_fam, ntoplot = 16, metadata = meta,"metodo")+
  theme_bw()+
  scale_fill_manual(values = pal)+
  theme(axis.text.x = element_text(angle = 90))+
ylab("Relative abundance (%)")+
  xlab("Sample")
b
#ggsave('figures/barplot_class.png',width = 12, height = 6, dpi = 300, plot =b)

sum_gen<- sum$Genus %>%
  rownames_to_column() %>%
  filter(!str_detect(rowname, "NA")) %>%
 filter(!str_detect(rowname, "uncultured")) %>%
  column_to_rownames(var="rowname")
taxa <- rownames(sum_gen)
new_taxa <- gsub(".*;", "", taxa)
rownames(sum_gen) <- new_taxa

g<-taxa_barplot(sum_gen, ntoplot = 16, metadata = meta,"metodo")+
  theme_bw()+
  scale_fill_manual(values = pal)+
  theme(axis.text.x = element_text(angle = 90))+
ylab("Relative abundance (%)")+
  xlab("Sample")
g
#ggsave('figures/barplot_genus.png',width = 12, height = 6, dpi = 300, plot =g)

# ============================================================
# 3. ALPHA DIVERSITY — Hill numbers
# ============================================================

otu_table <- tabla
q0 <- hill_div(otu_table, qvalue = 0)
q1 <- hill_div(otu_table, qvalue = 1)
q2 <- hill_div(otu_table, qvalue = 2)
q012 <- cbind(q0, q1, q2) %>% as.data.frame() %>% rownames_to_column(var = "SampleID")

#write.table(q012, file="data/q012_hilldiv.txt", sep = "\t")

Micro_div <- q012 %>%
  inner_join(meta, by = c("SampleID"="SampleID"))

Micro_div_mean <- Micro_div %>%
  filter(edad!="negativo") %>%
    summarise(
    q0_mean = mean(q0),
    q0_sd = sd(q0, na.rm = TRUE),
    q1_mean = mean(q1),
    q1_sd = sd(q1, na.rm = TRUE),
    q2_mean = mean(q2),
    q2_sd = sd(q2, na.rm = TRUE)
  ) %>%
  mutate_if(is.numeric, ~ round(., digits = 1))

Micro_div_sum <- Micro_div %>%
  filter(edad!="negativo") %>%
  group_by(metodo) %>%
    summarise(
    q0_mean = mean(q0),
    q0_sd = sd(q0, na.rm = TRUE),
    q1_mean = mean(q1),
    q1_sd = sd(q1, na.rm = TRUE),
    q2_mean = mean(q2),
    q2_sd = sd(q2, na.rm = TRUE)
  ) %>%
  mutate_if(is.numeric, ~ round(., digits = 1))

q0.p<-Micro_div %>%
 ggpaired(x = "metodo", y="q0", fill = "metodo", line.color = "black", line.size = 0.3)+
  ylab("Effective number of total ASVs")+
  xlab(element_blank())+
  scale_fill_manual(values = c("#E69F00","#44803F"))+
      theme_classic()+
   theme(legend.position = "none")+
  stat_compare_means(paired = TRUE)

q1.p<-Micro_div %>%
 ggpaired(x = "metodo", y="q1", fill = "metodo", line.color = "black", line.size = 0.3)+
  ylab("Effective number of frequent ASVs")+
  xlab(element_blank())+
  scale_fill_manual(values = c("#E69F00","#44803F"))+
      theme_classic()+
   theme(legend.position = "none")+
  stat_compare_means(paired = TRUE)

q2.p<-Micro_div %>%
 ggpaired(x = "metodo", y="q2", fill = "metodo", line.color = "black", line.size = 0.3)+
  ylab("Effective number of dominant ASVs")+
  xlab(element_blank())+
  scale_fill_manual(values = c("#E69F00","#44803F"))+
      theme_classic()+
   theme(legend.position = "none")+
  stat_compare_means(paired = TRUE)

alpha.div<- plot_grid(q0.p, q1.p, q2.p,
                          nrow = 1)

ggsave('figures/boxplot_alpha_div.png', width = 7.4, height = 3, dpi = 300, plot =alpha.div)

alpha.div

# ============================================================
# 3.4 ALPHA DIVERSITY — Phylogenetic Hill numbers
# ============================================================

filename= "data/tree.nwk"
# para calcular diversidad alfa
#cols --- species     rows ---- sites

comm =t(tabla)
IDs= colnames(comm)   # para obtener los nombres de las sp
tree = ape::read.tree(filename)   # leer el árbol
tree= ape::rtree(n=ncol(comm), tip.label = paste0(IDs))  # tener el árbol enraizado, con los nombres de los sitios y los nombres de las sp.

div0=hill_phylo(comm, tree, q = 0)
div1=hill_phylo(comm, tree, q = 1)
div2=hill_phylo(comm, tree, q = 2)

phy012 <- cbind(div0, div1, div2) %>% as.data.frame() %>% rownames_to_column(var = "SampleID")

write.table(phy012, file="data/phy012_hilldiv.txt", sep = "\t")

Micro_phy <- phy012 %>%
  inner_join(meta, by = c("SampleID"="SampleID"))

Micro_phy_mean <- Micro_phy %>%
  filter(edad!="negativo") %>%
  summarise(
    phy_q0_mean = mean(div0),
    phy_q0_sd = sd(div0, na.rm = TRUE),
    phy_q1_mean = mean(div1),
    phy_q1_sd = sd(div1, na.rm = TRUE),
    phy_q2_mean = mean(div2),
    phy_q2_sd = sd(div2, na.rm = TRUE)
  ) %>%
  mutate_if(is.numeric, ~ round(., digits = 1))

Micro_phy_sum <- Micro_phy %>%
  filter(edad!="negativo") %>%
  group_by(metodo) %>%
    summarise(
    phy_q0_mean = mean(div0),
    phy_q0_sd = sd(div0, na.rm = TRUE),
    phy_q1_mean = mean(div1),
    phy_q1_sd = sd(div1, na.rm = TRUE),
    phy_q2_mean = mean(div2),
    phy_q2_sd = sd(div2, na.rm = TRUE)
  ) %>%
  mutate_if(is.numeric, ~ round(., digits = 1))

q0.p<-Micro_phy %>%
 ggpaired(x = "metodo", y="div0", fill = "metodo", line.color = "black", line.size = 0.3)+
  ylab("Effective total branch lenght")+
  xlab(element_blank())+
  scale_fill_manual(values = c("#E69F00","#44803F"))+
      theme_classic()+
   theme(legend.position = "none")+
  stat_compare_means(paired = TRUE)

q1.p<-Micro_phy %>%
 ggpaired(x = "metodo", y="div1", fill = "metodo", line.color = "black", line.size = 0.3)+
  ylab("Effective frequent branch lenght")+
  xlab(element_blank())+
  scale_fill_manual(values = c("#E69F00","#44803F"))+
      theme_classic()+
   theme(legend.position = "none")+
  stat_compare_means(paired = TRUE)

q2.p<-Micro_div %>%
 ggpaired(x = "metodo", y="q2", fill = "metodo", line.color = "black", line.size = 0.3)+
  ylab("Effective dominant branch lenght")+
  xlab(element_blank())+
  scale_fill_manual(values = c("#E69F00","#44803F"))+
      theme_classic()+
   theme(legend.position = "none")+
  stat_compare_means(paired = TRUE)

alpha.phy<- plot_grid(q0.p, q1.p, q2.p,
                           nrow = 1)

ggsave('figures/boxplot_alpha_phy.png', width = 7.4, height = 3, dpi = 300, plot =alpha.phy)

alpha.phy

# ============================================================
# 4. BETA DIVERSITY ANALYSIS
# ============================================================

otutable= as.data.frame(otu_table)
hill_pair_dis <- hillpair(data = otutable, q=c(0,1,2))

# Order q=0
hill_0_dis_nmds <- hill_pair_dis$q0S %>%
  metaMDS(.,trymax = 500, k = 2, verbosity = FALSE) %>%
  vegan:::scores.metaMDS() %>%
  as_tibble(., rownames = "sample")

hill_0_dis_nmds <- hill_0_dis_nmds %>%
  left_join(meta, by = join_by(sample == SampleID)) %>%
  group_by(metodo) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup()

# Order q=1
hill_1_dis_nmds <- hill_pair_dis$q1S %>%
  metaMDS(.,trymax = 500, k = 2, verbosity = FALSE) %>%
  vegan:::scores.metaMDS() %>%
  as_tibble(., rownames = "sample")

hill_1_dis_nmds <- hill_1_dis_nmds %>%
  left_join(meta, by = join_by(sample == SampleID)) %>%
  group_by(metodo) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup()

# Order q=2
hill_2_dis_nmds <- hill_pair_dis$q2S %>%
  metaMDS(.,trymax = 500, k = 2, verbosity = FALSE) %>%
  vegan:::scores.metaMDS() %>%
  as_tibble(., rownames = "sample")

hill_2_dis_nmds <- hill_2_dis_nmds %>%
  left_join(meta, by = join_by(sample == SampleID)) %>%
  group_by(metodo) %>%
  mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
  mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
  ungroup()

NMDS_Plot_q0 <- ggplot(hill_0_dis_nmds, aes(x = NMDS1, y = NMDS2,
                                               color = metodo)) +
  geom_point(size = 6) +
  geom_segment(aes(x = x_cen, y = y_cen, xend = NMDS1, yend = NMDS2), alpha = 0.2) +
  theme_classic() +
    labs(title = "Order q=0") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(-0.7, 1.1))+
  scale_colour_manual(values = c("#E69F00","#44803F"))+
  annotate("text", x = -1, y = 1, label = expression("perMANOVA: F=0.99; R"^2*"=0.09; p=1"),
          color = "black", size = 4, fontface = "bold", hjust = 0, vjust = 0)

NMDS_Plot_q1 <- ggplot(hill_1_dis_nmds, aes(x = NMDS1, y = NMDS2,
                                               color = metodo)) +
  geom_point(size = 6) +
  geom_segment(aes(x = x_cen, y = y_cen, xend = NMDS1, yend = NMDS2), alpha = 0.2) +
  theme_classic() +
    labs(title = "Order q=1") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))+
  scale_x_continuous(limits = c(-0.6, 0.6))+
  scale_y_continuous(limits = c(-0.5, 0.6))+
  scale_colour_manual(values = c("#E69F00","#44803F"))+
 annotate("text", x = -0.6, y = 0.5, label = expression("perMANOVA: F=2.94; R"^2*"=0.22; p=0.031"),
          color = "black", size = 4, fontface = "bold", hjust = 0, vjust = 0)

NMDS_Plot_q2 <- ggplot(hill_2_dis_nmds, aes(x = NMDS1, y = NMDS2,
                                               color = metodo)) +
  geom_point(size = 6) +
  geom_segment(aes(x = x_cen, y = y_cen, xend = NMDS1, yend = NMDS2), alpha = 0.2) +
  theme_classic() +
    labs(title = "Order q=2") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
  scale_x_continuous(limits = c(-0.6, 0.6))+
  scale_y_continuous(limits = c(-0.5, 0.6))+
    scale_colour_manual(values = c("#E69F00","#44803F"))+
  annotate("text", x = -0.6, y = 0.5, label = expression("perMANOVA: F=3.68; R"^2*"=0.27; p=0.031"),
          color = "black", size = 4, fontface = "bold", hjust = 0, vjust = 0)

beta.div<- plot_grid(NMDS_Plot_q0, NMDS_Plot_q1, NMDS_Plot_q2,
                          ncol = 1,
                      rel_heights = c(0.9,0.9,1.1))
ggsave("figures/NMDS_Plots.jpeg", width = 4, height = 9.0, dpi = 300, plot = beta.div)

beta.div

set.seed(123)
metadata<-Micro_div %>% column_to_rownames(var = "SampleID") %>% rename(Método="metodo")
perm_q0 <- adonis2((hill_pair_dis$q0S %>% as.dist) ~ Método,
                   strata = metadata$Individuo,
                   data = metadata, permutations = 999) %>%
  round(., digits = 3)

perm_q1 <- adonis2((hill_pair_dis$q1S %>% as.dist) ~ Método,
                   strata = metadata$Individuo,
                   data = metadata, permutations = 999) %>%
  round(., digits = 3)

perm_q2 <- adonis2((hill_pair_dis$q2S %>% as.dist) ~ Método,
                   strata = metadata$Individuo,
                   data = metadata, permutations = 999) %>%
  round(., digits = 3)

per0<-data.frame(perm_q0, check.names = F)%>%
  replace(is.na(.), "-") %>%
  rownames_to_column(var="Factor") %>%
  ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(4), row.side = "bottom",
                linewidth = 3, linetype = 1) %>%
  tab_add_title(text = "Diversidad beta al orden q = 0",
                face = "plain", size = 10)

per1<-data.frame(perm_q1, check.names = F)%>%replace(is.na(.), "-") %>%
  rownames_to_column(var="Factor") %>%
  ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(4), row.side = "bottom",
                linewidth = 3, linetype = 1) %>%
  tab_add_title(text = "Diversidad beta al orden q = 1",
                face = "plain", size = 10)

per2<-data.frame(perm_q2, check.names = F)%>%
  replace(is.na(.), "-") %>%
  rownames_to_column(var="Factor") %>%
  ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(4), row.side = "bottom",
                linewidth = 3, linetype = 1) %>%
  tab_add_title(text = "Diversidad beta al orden q = 2",
                face = "plain", size = 10)

beta.divs<- plot_grid(NMDS_Plot_q0, NMDS_Plot_q1, NMDS_Plot_q2,
                      per0, per1, per2,
                          labels = "AUTO", nrow = 2, ncol = 3,
                      rel_widths = c(0.8,0.8,1.2),
                     rel_heights = c(1,0.4))
#ggsave("figures/NMDS_Plots_w_perm.jpeg", width = 13, height = 5.5, dpi = 300, plot = beta.divs)

beta.divs

disp.species.0 = betadisper((hill_pair_dis$q0S %>% as.dist),
                            metadata$Método)
dis.per.0<-permutest(disp.species.0)
disp.species.0$centroids

disp.species.1 = betadisper((hill_pair_dis$q1S %>% as.dist),
                            metadata$Método)
dis.per.1<-permutest(disp.species.1)

disp.species.2 = betadisper((hill_pair_dis$q2S %>% as.dist),
                            metadata$Método)
dis.per.2<-permutest(disp.species.2)

par(mfrow = c(1, 3))
boxplot(disp.species.0, ylab= "Distancia del centroide", xlab = "")
boxplot(disp.species.1, ylab= "Distancia del centroide", xlab = "")
boxplot(disp.species.2, ylab= "Distancia del centroide", xlab = "")

disp.0<-data.frame(dis.per.0$tab, check.names = F) %>%
  replace(is.na(.), "-") %>%
    mutate(F = as.numeric(F)) %>%
    mutate_if(is.numeric, ~ifelse(is.na(.), NA,
                                  round(., digits = 3))) %>%
  rownames_to_column(var="Factor") %>%
  ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(3), row.side = "bottom",
                linewidth = 3, linetype = 1) %>%
  tab_add_title(text = "Diversidad beta al orden q = 0",
                face = "plain", size = 10)

 disp.1<-data.frame(dis.per.1$tab, check.names = F) %>%
   replace(is.na(.), "-") %>%
   mutate(F = as.numeric(F)) %>%
   mutate_if(is.numeric, ~ifelse(is.na(.), NA,
                                 round(., digits = 3))) %>%
   rownames_to_column(var="Factor") %>%
   ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
   tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
   tab_add_hline(at.row = c(3), row.side = "bottom",
                 linewidth = 3, linetype = 1) %>%
   tab_add_title(text = "Diversidad beta al orden q = 1",
                 face = "plain", size = 10)

 disp.2<-data.frame(dis.per.2$tab, check.names = F) %>%
   replace(is.na(.), "-") %>%
   mutate(F = as.numeric(F)) %>%
   mutate_if(is.numeric, ~ifelse(is.na(.), NA,
                                 round(., digits = 3))) %>%
   rownames_to_column(var="Factor") %>%
   ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
   tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
   tab_add_hline(at.row = c(3), row.side = "bottom",
                 linewidth = 3, linetype = 1) %>%
   tab_add_title(text = "Diversidad beta al orden q = 2",
                 face = "plain", size = 10)

dispers<- plot_grid(disp.0, disp.1, disp.2, nrow = 1, ncol = 3)
#ggsave("figures/dispersion_anova.jpeg", width = 14, height = 2, dpi = 300, plot = dispers)

dispers

# ============================================================
# 4.4 PHYLOGENETIC BETA DIVERSITY
# NOTE: dist_matrix.0/1/2 are used below but not defined in the original script.
# ============================================================

NMDS_Plot_q0 <- ggplot(hill_0_dis_nmds, aes(x = NMDS1, y = NMDS2,
                                               color = metodo)) +
  geom_point(size = 6) +
  geom_segment(aes(x = x_cen, y = y_cen, xend = NMDS1, yend = NMDS2), alpha = 0.2) +
  theme_classic() +
    labs(title = "Order q=0") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))+
  scale_x_continuous(limits = c(-0.2, 0.18))+
  scale_y_continuous(limits = c(-0.15, 0.2))+
  scale_colour_manual(values = c("#E69F00","#44803F"))+
  annotate("text", x = -0.2, y = 0.15, label = expression("perMANOVA: F=1.90; R"^2*"=0.16; p=0.031"),
          color = "black", size = 4, fontface = "bold", hjust = 0, vjust = 0)

NMDS_Plot_q1 <- ggplot(hill_1_dis_nmds, aes(x = NMDS1, y = NMDS2,
                                               color = metodo)) +
  geom_point(size = 6) +
  geom_segment(aes(x = x_cen, y = y_cen, xend = NMDS1, yend = NMDS2), alpha = 0.2) +
  theme_classic() +
    labs(title = "Order q=1") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))+
  scale_x_continuous(limits = c(-0.2, 0.3))+
  scale_y_continuous(limits = c(-0.2, 0.2))+
  scale_colour_manual(values = c("#E69F00","#44803F"))+
  annotate("text", x = -0.2, y = 0.15, label = expression("perMANOVA: F=3.61; R"^2*"=0.26; p=0.031"),
          color = "black", size = 4, fontface = "bold", hjust = 0, vjust = 0)

NMDS_Plot_q2 <- ggplot(hill_2_dis_nmds, aes(x = NMDS1, y = NMDS2,
                                               color = metodo)) +
  geom_point(size = 6) +
  geom_segment(aes(x = x_cen, y = y_cen, xend = NMDS1, yend = NMDS2), alpha = 0.2) +
  theme_classic() +
    labs(title = "Order q=2") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
  scale_x_continuous(limits = c(-0.2, 0.3))+
  scale_y_continuous(limits = c(-0.15, 0.18))+
    scale_colour_manual(values = c("#E69F00","#44803F"))+
  annotate("text", x = -0.2, y = 0.15, label = expression("perMANOVA: F=2.72; R"^2*"=0.21; p=0.062"),
          color = "black", size = 4, fontface = "bold", hjust = 0, vjust = 0)

beta.div<- plot_grid(NMDS_Plot_q0, NMDS_Plot_q1, NMDS_Plot_q2,
                       ncol = 1,
                      rel_heights = c(0.9,0.9,1.1))
ggsave("figures/NMDS_Plots_phylo.jpeg", width = 4, height = 9.0, dpi = 300, plot = beta.div)

beta.div

set.seed(123)
metadata<-Micro_phy %>% column_to_rownames(var = "SampleID") %>% rename(Método="metodo")

perm_q0 <- adonis2(((1-dist_matrix.0) %>% as.dist) ~ Método,
                   strata = metadata$Individuo,
                   data = metadata, permutations = 999) %>%
  round(., digits = 3)

perm_q1 <- adonis2(((1-dist_matrix.1) %>% as.dist) ~ Método,
                   strata = metadata$Individuo,
                   data = metadata, permutations = 999) %>%
  round(., digits = 3)

perm_q2 <- adonis2(((1-dist_matrix.2) %>% as.dist) ~ Método,
                   strata = metadata$Individuo,
                   data = metadata, permutations = 999) %>%
  round(., digits = 3)

per0<-data.frame(perm_q0, check.names = F)%>%
  replace(is.na(.), "-") %>%
  rownames_to_column(var="Factor") %>%
  ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(4), row.side = "bottom",
                linewidth = 3, linetype = 1) %>%
  tab_add_title(text = "Phylogenetic beta diversity at order q = 0",
                face = "plain", size = 10)

per1<-data.frame(perm_q1, check.names = F)%>%replace(is.na(.), "-") %>%
  rownames_to_column(var="Factor") %>%
  ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(4), row.side = "bottom",
                linewidth = 3, linetype = 1) %>%
  tab_add_title(text = "Phylogenetic beta diversity at order q = 1",
                face = "plain", size = 10)

per2<-data.frame(perm_q2, check.names = F)%>%
  replace(is.na(.), "-") %>%
  rownames_to_column(var="Factor") %>%
  ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(4), row.side = "bottom",
                linewidth = 3, linetype = 1) %>%
  tab_add_title(text = "Phylogenetic beta diversity at order q = 2",
                face = "plain", size = 10)

beta.divs<- plot_grid(NMDS_Plot_q0, NMDS_Plot_q1, NMDS_Plot_q2,
                      per0, per1, per2,
                          labels = "AUTO", nrow = 2, ncol = 3,
                      rel_widths = c(0.9,0.9,1.1),
                     rel_heights = c(1,0.4))
#ggsave("figures/NMDS_Plots_w_perm_phylogenetic.jpeg", width = 13, height = 5.5, dpi = 300, plot = beta.divs)

beta.divs

disp.species.0 = betadisper(((1-dist_matrix.0) %>% as.dist),
                            metadata$Método)
dis.per.0<-permutest(disp.species.0)
disp.species.0$centroids

disp.species.1 = betadisper(((1-dist_matrix.1) %>% as.dist),
                            metadata$Método)
dis.per.1<-permutest(disp.species.1)

disp.species.2 = betadisper(((1-dist_matrix.2) %>% as.dist),
                            metadata$Método)
dis.per.2<-permutest(disp.species.2)

par(mfrow = c(1, 3))
boxplot(disp.species.0, ylab= "Distance from the centroide", xlab = "")
boxplot(disp.species.1, ylab= "Distance from the centroide", xlab = "")
boxplot(disp.species.2, ylab= "Distance from the centroide", xlab = "")

disp.0<-data.frame(dis.per.0$tab, check.names = F) %>%
  replace(is.na(.), "-") %>%
    mutate(F = as.numeric(F)) %>%
    mutate_if(is.numeric, ~ifelse(is.na(.), NA,
                                  round(., digits = 3))) %>%
  rownames_to_column(var="Factor") %>%
  ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = c(3), row.side = "bottom",
                linewidth = 3, linetype = 1) %>%
  tab_add_title(text = "Diversidad beta al orden q = 0",
                face = "plain", size = 10)

 disp.1<-data.frame(dis.per.1$tab, check.names = F) %>%
   replace(is.na(.), "-") %>%
   mutate(F = as.numeric(F)) %>%
   mutate_if(is.numeric, ~ifelse(is.na(.), NA,
                                 round(., digits = 3))) %>%
   rownames_to_column(var="Factor") %>%
   ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
   tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
   tab_add_hline(at.row = c(3), row.side = "bottom",
                 linewidth = 3, linetype = 1) %>%
   tab_add_title(text = "Diversidad beta al orden q = 1",
                 face = "plain", size = 10)

 disp.2<-data.frame(dis.per.2$tab, check.names = F) %>%
   replace(is.na(.), "-") %>%
   mutate(F = as.numeric(F)) %>%
   mutate_if(is.numeric, ~ifelse(is.na(.), NA,
                                 round(., digits = 3))) %>%
   rownames_to_column(var="Factor") %>%
   ggtexttable(., rows = NULL, theme = ttheme("blank")) %>%
   tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
   tab_add_hline(at.row = c(3), row.side = "bottom",
                 linewidth = 3, linetype = 1) %>%
   tab_add_title(text = "Diversidad beta al orden q = 2",
                 face = "plain", size = 10)

dispers<- plot_grid(disp.0, disp.1, disp.2, nrow = 1, ncol = 3)
#ggsave("figures/dispersion_anova_phylogenetic.jpeg", width = 14, height = 2, dpi = 300, plot = dispers)

dispers

# ============================================================
# 5. DIFFERENTIAL ABUNDANCE
# ============================================================

# phyloseq object considering structual zeros
SAM <- meta %>%
  column_to_rownames("SampleID") %>%
    sample_data() # convert to phyloseq sample_data object
SAM$Individuo <- as.factor(SAM$Individuo)
SAM$method <- as.factor(SAM$method)

ASV <- tabla %>%
#  mutate_all(~ replace(., . == 0, 0.00001)) %>%
# mutate_all(~./sum(.)) %>%
#  as.matrix() %>%
  otu_table(., taxa_are_rows = TRUE)

txx<- rownames(tabla)

TAX <- tax %>%
  filter(Feature.ID %in% txx) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix() %>%
  tax_table() # convert to phyloseq tax_table object

phyloseq <- merge_phyloseq(ASV, TAX, SAM)

differential_abundance <- ancombc2(
  data = phyloseq,
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "method",
  rand_formula = "(1 | Individuo)",
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  trend_control = NULL
)

# Save differential abundance to data object
save(differential_abundance.p,
     differential_abundance.c,
     differential_abundance.o,
     differential_abundance.f,
     differential_abundance.g,
     file = "data/ancom.Rdata")

load("data/ancom.Rdata")

tax <- data.frame(phyloseq@tax_table) %>%
  rownames_to_column(., "taxon")

df.method.table<-differential_abundance.g$res %>%
   filter(!grepl("Family:", taxon)) %>%
     filter(!grepl("uncultured", taxon)) %>%
   mutate(taxon = gsub("Genus:", "", taxon)) %>%
  dplyr::select(taxon=taxon, lfc_methodk, p_methodk) %>%
  filter(p_methodk < 0.05) %>%
  dplyr::arrange(lfc_methodk) %>%
  left_join(tax %>% dplyr::select(-c(taxon,Kingdom,Species)) %>% base::unique(), by = join_by(taxon==Genus))

d.a.g<- df.method.table %>%
      mutate(taxon=factor(taxon,levels=df.method.table$taxon)) %>%
      ggplot(aes(x = lfc_methodk, y = forcats::fct_rev(taxon), fill = Phylum)) +
        geom_col(size = 2) +
        scale_fill_manual(values = c("#D72542", "#F3C50D", "#2D3FA6")) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 4, linetype="dashed", color = "grey", linewidth=1) +
        geom_vline(xintercept = -4, linetype="dashed", color = "grey", linewidth=1) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.text.y = element_text(face = "italic"),
          strip.background = element_blank(),
          strip.text = element_blank()
        ) +
        labs(x="Log2fold change between methods", y="Genus")

#ggsave('figures/ancom_gen.png', width = 6, height = 1.8, dpi = 300, plot =d.a.g)

d.a.g

df.method.table.f<-differential_abundance.f$res %>%
   filter(!grepl("Order:", taxon)) %>%
     filter(!grepl("uncultured", taxon)) %>%
   mutate(taxon = gsub("Family:", "", taxon)) %>%
  dplyr::select(taxon=taxon, lfc_methodk, p_methodk) %>%
  filter(p_methodk < 0.05) %>%
  dplyr::arrange(lfc_methodk) %>%
  left_join(tax %>% dplyr::select(-c(taxon,Kingdom,Species,Genus)) %>%
              base::unique(), by = join_by(taxon==Family))

d.a.f<- df.method.table.f %>%
      mutate(taxon=factor(taxon,levels=df.method.table.f$taxon)) %>%
      ggplot(aes(x = lfc_methodk, y = forcats::fct_rev(taxon), fill = Phylum)) +
        geom_col(size = 2) +
        scale_fill_manual(values = c("#D72542", "#F3C50D", "#2D3FA6")) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 4, linetype="dashed", color = "grey", linewidth=1) +
        geom_vline(xintercept = -4, linetype="dashed", color = "grey", linewidth=1) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.text.y = element_text(),
          strip.background = element_blank(),
          strip.text = element_blank()
        ) +
        labs(x="Log2fold change between methods", y="Family")
#ggsave('figures/ancom_fam.png', width = 6, height = 1.6, dpi = 300, plot =d.a.f)

d.a.f

df.method.table.o<-differential_abundance.o$res %>%
   filter(!grepl("Class:", taxon)) %>%
     filter(!grepl("uncultured", taxon)) %>%
   mutate(taxon = gsub("Order:", "", taxon)) %>%
  dplyr::select(taxon=taxon, lfc_methodk, p_methodk) %>%
  filter(p_methodk < 0.05) %>%
  dplyr::arrange(lfc_methodk) %>%
  left_join(tax %>% dplyr::select(-c(taxon,Kingdom,Species, Genus, Family)) %>%
              base::unique(), by = join_by(taxon==Order))

d.a.o<- df.method.table.o %>%
      mutate(taxon=factor(taxon,levels=df.method.table$taxon)) %>%
      ggplot(aes(x = lfc_methodk, y = forcats::fct_rev(taxon), fill = Phylum)) +
        geom_col(size = 2) +
        scale_fill_manual(values = c("#D72542", "#F3C50D", "#2D3FA6")) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 4, linetype="dashed", color = "grey", linewidth=1) +
        geom_vline(xintercept = -4, linetype="dashed", color = "grey", linewidth=1) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.text.y = element_text(),
          strip.background = element_blank(),
          strip.text = element_blank()
        ) +
        labs(x="Log2fold change between methods", y="Order")

#ggsave('figures/ancom_ord.png', width = 6, height = 1.5, dpi = 300, plot =d.a.o)

d.a.o

df.method.table.c<-differential_abundance.c$res %>%
   filter(!grepl("Phylum:", taxon)) %>%
     filter(!grepl("uncultured", taxon)) %>%
   mutate(taxon = gsub("Class:", "", taxon)) %>%
  dplyr::select(taxon=taxon, lfc_methodk, p_methodk) %>%
  filter(p_methodk < 0.05) %>%
  dplyr::arrange(lfc_methodk) %>%
  left_join(tax %>% dplyr::select(-c(taxon,Kingdom,Species, Genus, Family, Order)) %>% base::unique(), by = join_by(taxon==Class))

df.method.table.p<-differential_abundance.p$res %>%
  dplyr::select(taxon=taxon, lfc_methodk, p_methodk) %>%
  filter(p_methodk < 0.05) %>%
  dplyr::arrange(lfc_methodk) %>% as_tibble()
