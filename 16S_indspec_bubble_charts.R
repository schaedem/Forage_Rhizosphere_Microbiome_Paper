library(tidyverse)
library(rebus)
library(stats) #for cor
library(Hmisc) #for rcorr
library(ggcorrplot)
library(ggpubr)
library(pBrackets) #for adding brackets to tile plot
library(grid) #for helping with brackets
library(ggtext) #for element_markdown

#### Formatting OTU count data ####

#Read in metadata
meta <- read_csv("final_samps.csv") %>%
  mutate(treatment = as.character(treatment),
         treatment = str_replace_all(treatment, "�", ""),
         species = as.factor(species),
         treatment = ifelse(grepl("Brachiaria v. Mulato II", treatment), "Brachiaria cv. Mulato II",
                            ifelse(grepl("Maize", treatment), "Maize monocrop", treatment)),
         treatment = as.character(treatment),
         intercrop = ifelse(grepl(PLUS, treatment), "yes", "no"))

indspec <- read_csv("indicator_species_taxonomy_2.csv")
class <- read_rds("16S_trophic_class.rds")

rel_abund <- readRDS("otu_rare_relabund_16S.rds") %>%
  dplyr::select(sample, otu, rel_abund, location, 
                treatment, species, timepoint, block)
glimpse(rel_abund)
#joining only within species and location
indspec_full <- 
  semi_join(rel_abund, indspec, by=c("otu", "location", "species")) %>%
  arrange(sample) %>%
  dplyr::select(sample, otu, rel_abund, species, location, timepoint)

#### Formatting otu labels ####
graph_abund_df <- indspec_full %>%
  #Maize
  mutate(tax_otu_label = 
           ifelse(otu == "Otu00381", "Unclassified *Bacteria* (OTU 381)", #Burera
          ifelse(otu == "Otu00571", "Unclassified *Bacteria* (OTU 571)", #Rubona
           ifelse(otu == "Otu01065", "Unclassified *Sphingomonadaceae*<br>(OTU 1065)", #Rubona
          #Napier
          ifelse(otu == "Otu01452", "*Phenylobacterium* (OTU 1452)", #Burera; aerobic soil-associated
          ifelse(otu == "Otu01953", "*Gemmatimonas* (OTU 1953)", #Burera
          ifelse(otu == "Otu00316", "Unclassified *Micrococcaceae*<br>(OTU 316)", #Karama
          ifelse(otu == "Otu00597", "*Verrucomicrobia* Subdivision 3<br>(OTU 597)", #K
          ifelse(otu == "Otu01297", "Unclassified *Betaproteobacteria*<br>(OTU 1297)", #K
          ifelse(otu == "Otu01653", "Unclassified *Bacteria* (OTU 1653)", #Karama
          ifelse(otu == "Otu00516", "*Gemmata* (OTU 516)", #Rubona
          ifelse(otu == "Otu00690", "Unclassified *Bacteria* (OTU 690)",#Rubona
          ifelse(otu == "Otu01116", "Unclassified *Planctomycetaceae*<br>(OTU 1116)",#Rubona
          ifelse(otu == "Otu01444", "Unclassified *Actinomycetales*<br>(OTU 1444)",#Rubona
           #Brachiaria
          ifelse(otu == "Otu00584", "*Acidobacteria* Gp6 (OTU 584)", #burera
          ifelse(otu == "Otu01574", "*Paenibacillus* (OTU 1574)", #Burera
          ifelse(otu == "Otu02740", "Unclassified *Ktedonobacteria*<br>(OTU 2740)", #Burera
          ifelse(otu == "Otu03228", "Unclassified *Bacteria* (OTU 3228)", #Burera
          ifelse(otu == "Otu03700", "*Singulisphaera* (OTU 3700)", #Burera
          ifelse(otu == "Otu00138", "*Conexibacter* (OTU 138)",#Karama;
          ifelse(otu == "Otu00364", "Unclassified *Bacteria* (OTU 364)", #Karama
          ifelse(otu == "Otu00861", "Unclassified *Clostridiales*<br>(OTU 861)", #Rubona
          ifelse(otu == "Otu01256", "*Aciditerrimonas* (OTU 1256)", #Rubona
          ifelse(otu == "Otu01388", "Unclassified *Bacteria* (OTU 1388)", #Rubona
         #Desmodium
          ifelse(otu == "Otu00472", "*Hyphomicrobium* (OTU 472)", 
          NA)))))))))))))))))))))))))

graph_abund_df$tax_otu_label <- factor(graph_abund_df$tax_otu_label, levels =c(
  #Maize
  "Unclassified *Bacteria* (OTU 381)","Unclassified *Bacteria* (OTU 571)","Unclassified *Sphingomonadaceae*<br>(OTU 1065)",
  #Napier
  "*Phenylobacterium* (OTU 1452)", "*Gemmatimonas* (OTU 1953)", 
  "Unclassified *Micrococcaceae*<br>(OTU 316)", "*Verrucomicrobia* Subdivision 3<br>(OTU 597)",
  "Unclassified *Betaproteobacteria*<br>(OTU 1297)","Unclassified *Bacteria* (OTU 1653)",
  "*Gemmata* (OTU 516)", "Unclassified *Bacteria* (OTU 690)", 
  "Unclassified *Planctomycetaceae*<br>(OTU 1116)","Unclassified *Actinomycetales*<br>(OTU 1444)",
  #Brachiaria
  "*Acidobacteria* Gp6 (OTU 584)", "*Paenibacillus* (OTU 1574)",
  "Unclassified *Ktedonobacteria*<br>(OTU 2740)","Unclassified *Bacteria* (OTU 3228)",
  "*Singulisphaera* (OTU 3700)","*Conexibacter* (OTU 138)", "Unclassified *Bacteria* (OTU 364)",
  "Unclassified *Clostridiales*<br>(OTU 861)", "*Aciditerrimonas* (OTU 1256)","Unclassified *Bacteria* (OTU 1388)",
  #Desmodium
  "*Hyphomicrobium* (OTU 472)"
))

final_graph_abund <- graph_abund_df %>%
  group_by(timepoint, otu) %>%
  summarise(mean_rel_abund = mean(rel_abund)) %>%
  ungroup() %>%
  inner_join(graph_abund_df, by= c("otu", "timepoint")) %>%
  dplyr::select(timepoint, location, species, mean_rel_abund, location, tax_otu_label) %>%
  unique() %>%
  mutate(rel_abund_trans = 100*(mean_rel_abund + 1/4000))


#### Graphing #####
pal <- c("#1F78B4", "#6A3D9A", "#33A02C", "#E31A1C")

bubble_chart <- ggplot(final_graph_abund, 
                aes(y=tax_otu_label, x=timepoint, color=species, size=rel_abund_trans)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values = pal, guide = "none") +
  theme_bw() +
  theme(axis.text.y = element_markdown(size=15),
        legend.position = "right",
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(face = "bold", size=15),
        ) +
  ylab("") +
  xlab("Collection Time") +
  guides(size = guide_legend(title.position = "right", 
                             title = "Mean Relative\nAbundance (%)",
                             override.aes = list(color="black", shape=21),
                             title.theme = element_text(
                               size=15, face="bold", angle=-90),
                             label.theme = element_text(size=12)))

  
bubble_chart


setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/04_Station_BNF/final_data/seq_data/16S/Figures")
ggsave("partial_cor_graph_indspec_abundance.png",  dpi=300, height = 14, width=7, units="in")
