library(tidyverse)
library(rebus)
library(stats) #for cor
library(Hmisc) #for rcorr
library(ggcorrplot)
library(pBrackets) #for adding brackets to tile plot
library(grid) #for helping with brackets
library(ggtext)

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

indspec_class <- inner_join(indspec, class, tax, by=c("otu", "location", "species")) %>%
  group_by(location, species, otu) %>%
  count(guild)

# rel_abund <- readRDS("otu_rare_relabund_16S.rds") %>%
#   dplyr::select(sample, otu, rel_abund, location, 
#                 treatment, species, timepoint, block)
#glimpse(rel_abund)

#joining only within species and location
indspec_full <- 
  inner_join(indspec, class, by=c("otu", "location", "species")) %>%
  mutate(otu = ifelse(otu == "Otu00043" & species == "Napier", "Otu00043_N",
               ifelse(otu == "Otu00043" & species == "Maize", "Otu00043_M",
                      otu))) %>%
  pivot_wider(names_from = otu, values_from=rel_abund, 
              values_fill = 0, id_cols=sample) %>%
  inner_join(meta, by="sample") %>%
  mutate(sample = str_trim(str_replace_all(sample, "R", ""))) %>%
  arrange(sample) %>%
  mutate(sample = ifelse(sample == "6221" | sample == "6222", 622, sample)) %>%
  unique() #values for the two reps of 622 are identical


#### Read in NIR data ####
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Thesis/Rwanda/Lab work/04_Station_BNF/final_data/NIR")

nir <- read_csv("nir_final.csv") %>%
  select(sample, protein, adf, ndf, lignin) %>%
  mutate(sample = str_trim(str_replace_all(sample, "L", ""))) %>%
  arrange(sample)

#for partial matrix
nir_partial <- semi_join(nir, indspec_full, by="sample") %>%
  arrange(sample)

missing_nir <- anti_join(indspec_full, nir_partial, by="sample")

#### Read in IRMS data ####
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Thesis/Rwanda/Lab work/04_Station_BNF/final_data/IRMS")

irms <- read_csv("irms_final_full.csv") %>%
  dplyr::select(sample, N.corrected, d15N.corrected) %>%
  mutate(sample = as.character(sample))

ndfa <- read_csv("final_good_leg_ndfa.csv") %>%
  dplyr::select(sample, ndfa)%>%
  mutate(sample = as.character(sample))


#Join to nir dataframe
nir_irms_partial <- inner_join(nir_partial, irms, by="sample") %>%
  filter(!is.na(sample)) %>%
  group_by(sample) %>%
  mutate(across(.cols=everything(), ~mean(.x))) %>% #sample 314 was doubled
  ungroup() %>%
  unique()

#Find missing samples
missing_irms <- anti_join(nir_partial, nir_irms_partial) %>%
  mutate(N.corrected = NA, 
         d15N.corrected = NA)

nir_irms_partial_final <- rbind(nir_irms_partial, missing_irms) 

nir_irms_partial_final_mat <- nir_irms_partial_final %>%
  arrange(sample) %>%
  column_to_rownames("sample") %>%
  as.matrix()


##### Make Correlation Matrix ####
#For partial dataset, within species and location
#in which the indicator species was identified
#stats::cor(x,y, method="spearman")

#need to calculate correlations for otu 43 separately for maize and napier

indspec_full_mat <- indspec_full %>%
  semi_join(nir_irms_partial_final, by="sample") %>%
  arrange(sample) %>%
  column_to_rownames("sample") %>%
  select(-c(location, block, timepoint, species, treatment, intercrop)) %>%
  as.matrix()

#cor_partial <- cor(indspec_full_mat, nir_partial_mat, method="spearman")
cor_partial_p <- rcorr(indspec_full_mat, nir_irms_partial_final_mat, type="spearman")

##### Visualize Matrices ####

#### 1. Partial Matrix ####
pal2 <- c("#A50026", "white","#313695")
dim(cor_partial_p$r)
r_mat <- cor_partial_p$r[21:26, 1:20] 
p_mat <- cor_partial_p$P[21:26, 1:20]

#### graphing with geom_tile to rearrange Otu order according to species and location
r_tibble <- as_tibble(r_mat, rownames="quality") %>%
  pivot_longer(!quality, names_to = "otu", values_to = "spearman")

p_tibble <- as_tibble(p_mat, rownames="quality") %>% 
  pivot_longer(!quality, names_to = "otu", values_to = "p.value")

graph_cor_df <- inner_join(r_tibble, p_tibble, by=c("otu", "quality")) %>%
  mutate(adj.p.value = round(p.adjust(p.value, method="BH"), 3),
         spearman = ifelse(adj.p.value > 0.05, NA, spearman))
graph_cor_df$quality <- factor(graph_cor_df$quality, levels = c(
  "N.corrected", "d15N.corrected","adf", "ndf", "lignin", "protein"
))

##### Formatting OTU labels ####
graph_cor_df <- graph_cor_df %>%
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


graph_cor_df$tax_otu_label <- factor(graph_cor_df$tax_otu_label, levels = c(
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
    
##### Partial correlation plot ####

partial_cor_graph <- ggplot(graph_cor_df, aes(y=tax_otu_label, x=quality)) +
  geom_tile(aes(fill=spearman), color="lightgrey") +
  scale_fill_gradient2(low = "#fcfdbf", high= "#b73779", mid="lightgrey",
                       na.value = "#FFFFFF", name = "Spearman's\nCorrelation") +
  theme_bw() +
  scale_x_discrete(position ="top",labels = c("N.corrected" = "N",
                                              "d15N.corrected" = expression(delta~""^15~"N"),
                                              "adf" = "ADF", "ndf" = "NDF",
                                              "lignin" = "Lignin",
                                              "protein" = "Protein")) +
  theme(legend.position="right",
        #plot.margin = margin(r=5, unit="cm"),
        axis.text.x = element_text(hjust=-0.01, angle=45, size=15, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size=15, face="bold"),
        axis.text.y = element_markdown(size=15)) 
partial_cor_graph + theme(axis.text.y = element_blank())


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Thesis/Rwanda/Lab work/04_Station_BNF/final_data/seq_data/ITS2/Figures")
ggsave("partial_cor_graph_indspec_ITS2.png", dpi=300, height = 14, width=5, units="in")

#### Ndfa in legumes ###

ndfa_full <- ndfa %>%
  group_by(sample) %>%
  mutate(ndfa = mean(ndfa)) %>%
  ungroup() %>%
  unique() %>%
  inner_join(indspec_full, by="sample") %>%
  arrange(sample)

ndfa_ndfa <- ndfa_full %>%
  dplyr::select(sample, ndfa) %>%
  column_to_rownames("sample") %>%
  as.matrix()

ndfa_indspec <- ndfa_full %>%
  dplyr::select(-c(ndfa, location, species, treatment, timepoint, block, intercrop)) %>%
  column_to_rownames("sample") %>%
  as.matrix()

ndfa_cor <- rcorr(ndfa_indspec, ndfa_ndfa, type="spearman")

ndfa_r_mat <- ndfa_cor$r[21, 1:20] 
ndfa_p_mat <- ndfa_cor$P[21, 1:20]


ndfa_r_tibble <- as_tibble(ndfa_r_mat, rownames="otu") %>%
  rename(spearman = value)

ndfa_p_tibble <- as_tibble(ndfa_p_mat, rownames="otu") %>% 
 rename(p.val = value)

graph_labs <- graph_cor_df %>%
  dplyr::select(otu, tax_otu_label)

ndfa_graph_cor_df <- inner_join(ndfa_r_tibble, ndfa_p_tibble, by=c("otu")) %>%
  filter(!is.nan(spearman)) %>%
  mutate(adj.p.val = round(p.adjust(p.val, method="BH"), 3),
         spearman = ifelse(adj.p.val > 0.05, NA, spearman),
         quality = "Ndfa") %>%
  inner_join(graph_labs, by="otu") %>%
  unique()

#### Graph Ndfa-indspec correlations ####
ndfa_cor_graph <- ggplot(ndfa_graph_cor_df, aes(y=tax_otu_label, x=quality)) +
  geom_tile(aes(fill=spearman), color="lightgrey") +
  scale_fill_gradient2(low = "#fcfdbf", high= "#b73779", mid="lightgrey",
                       na.value = "#FFFFFF", name = "Spearman's\nCorrelation") +
  theme_bw() +
  scale_x_discrete(position ="top") +
  theme(legend.position="right",
        #plot.margin = margin(r=5, unit="cm"),
        axis.text.x = element_text(hjust=-0.01, angle=45, size=15, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size=15, face="bold"),
        axis.text.y = element_markdown(size=15)) 
ndfa_cor_graph

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Thesis/Rwanda/Lab work/04_Station_BNF/final_data/seq_data/ITS2/Figures")
ggsave("ndfa_indspec_cor.png", dpi = 300, height = 5, width=7, units = "in")
