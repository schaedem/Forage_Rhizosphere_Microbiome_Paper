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

indspec <- read_csv("indicator_species_taxonomy_its2.csv")
class <- read_rds("rare_its2_guild.rds")
glimpse(class)
#joining only within species and location
indspec_full <- 
  semi_join(class, indspec, by=c("otu", "location", "species")) %>%
  arrange(sample) %>%
  dplyr::select(sample, otu, rel_abund, species, location, 
                timepoint, guild, growth, trait, confidence) %>%
  mutate(otu = ifelse(otu == "Otu00043" & species == "Napier", "Otu00043_N",
                      ifelse(otu == "Otu00043" & species == "Maize", "Otu00043_M",
                             otu)))

#### Formatting otu labels ####
graph_abund_df <- indspec_full %>%
  mutate(tax_otu_label = as.character(
    #Maize
    ifelse(otu == "Otu00022", "Unclassified *Fungi* (OTU 22)", #Burera
           ifelse(otu == "Otu00322", "Unclassified *Fungi* (OTU 322)", #Burera
                  ifelse(otu == "Otu00043_M", "*Acrocalymma* (OTU 43) ", #Karama
                         ifelse(otu == "Otu00044", "*Neoroussoella* (OTU 44)", #Karama
                                ifelse(otu == "Otu00234", "Unclassified *Fungi* (OTU 234)",  #Karama
                                       #Napier
                                       ifelse(otu == "Otu00043_N", "*Acrocalymma* (OTU 43)", #Karama
                                              ifelse(otu == "Otu00266", "Unclassified *Fungi* (OTU 266)", #Karama
                                                     ifelse(otu == "Otu00133", "*Penicillifer* (OTU 133)", #Rubona
                                                            ifelse(otu == "Otu00200", "*Microascus* (OTU 200)", #Rubona
                                                                   ifelse(otu == "Otu00276", "*Chaetomium* (OTU 276)", #Rubona
                                                                          ifelse(otu == "Otu00499", "Unclassified *Phaeosphaeriaceae*<br>(OTU 499)", #Rubona
                                                                                 ifelse(otu == "Otu00828", "*Dentiscutata* (OTU 828)", #Rubona
                                                                                        #Brachiaria
                                                                                        ifelse(otu == "Otu00073", "*Chordomyces* (OTU 73)", #Karama
                                                                                               ifelse(otu == "Otu00729", "Unclassified *Hypocreales*<br>(OTU 729)", #Karama
                                                                                                      #Desmodium
                                                                                                      ifelse(otu == "Otu00086", "*Thermomyces* (OTU 86)", #Burera
                                                                                                             ifelse(otu == "Otu00106", "Unclassified *Orbiliales*<br>(OTU 106)", #Karama
                                                                                                                    ifelse(otu == "Otu00865", "Unclassified *Chytridiomycota*<br>(OTU 865)", #Karama
                                                                                                                           ifelse(otu == "Otu00020", "*Pyrenochaetopsis* (OTU 20)", #Rubona
                                                                                                                                  ifelse(otu == "Otu00047", "*Striaticonidium* (OTU 47)", #Rubona
                                                                                                                                         ifelse(otu == "Otu00148", "*Pyrenochaetopsis* (OTU 148)", NA)))))))))))))))))))))) 


graph_abund_df$tax_otu_label <- factor(graph_abund_df$tax_otu_label, levels = c(
  #Maize
  "Unclassified *Fungi* (OTU 22)", "Otu00322", "Unclassified *Fungi* (OTU 322)", 
  "*Acrocalymma* (OTU 43) ", "*Neoroussoella* (OTU 44)", "Unclassified *Fungi* (OTU 234)",
  #Napier
  "*Acrocalymma* (OTU 43)",  "Unclassified *Fungi* (OTU 266)", "*Penicillifer* (OTU 133)", 
  "*Microascus* (OTU 200)", "*Chaetomium* (OTU 276)", "Unclassified *Phaeosphaeriaceae*<br>(OTU 499)", 
  "*Dentiscutata* (OTU 828)", 
  #Brachiaria
  "*Chordomyces* (OTU 73)", "Unclassified *Hypocreales*<br>(OTU 729)", 
  #Desmodium
  "*Thermomyces* (OTU 86)", "Unclassified *Orbiliales*<br>(OTU 106)", 
  "Unclassified *Chytridiomycota*<br>(OTU 865)", "*Pyrenochaetopsis* (OTU 20)", 
  "*Striaticonidium* (OTU 47)", "*Pyrenochaetopsis* (OTU 148)"))

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
  guides(size = guide_legend(title.position = "top", 
                             title = "Mean Relative\nAbundance (%)",
                             override.aes = list(color="black", shape=21),
                             title.theme = element_text(
                               size=15, face="bold"),
                             label.theme = element_text(size=12)))


bubble_chart


setwd("/Volumes/Backup_1 1/Marie/Thesis/Rwanda/Lab work/04_Station_BNF/final_data/seq_data/ITS2/Figures")
ggsave("partial_cor_graph_indspec_abundance_its2.png",  dpi=300, height = 14, width=7, units="in")
