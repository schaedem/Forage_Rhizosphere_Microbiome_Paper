library(tidyverse)
library(rebus)
library(stats) #for cor
library(Hmisc) #for rcorr
library(ggcorrplot)
library(ggpubr)
library(pBrackets) #for adding brackets to tile plot
library(grid) #for helping with brackets

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

indspec_class <- inner_join(indspec, class, tax, by="otu") %>%
  group_by(location, species, trophicMode) %>%
  count(trophicMode)

rel_abund <- readRDS("otu_rare_relabund_16S.rds") %>%
  dplyr::select(sample, otu, rel_abund, location, 
                treatment, species, timepoint, block)

#joining only within species and location
indspec_full <- 
  inner_join(indspec, rel_abund, by=c("otu", "location", "species")) %>%
  pivot_wider(names_from = otu, values_from=rel_abund, 
              values_fill = 0, id_cols=sample) %>%
  inner_join(meta, by="sample")%>%
  mutate(sample = str_trim(str_replace_all(sample, "R", ""))) %>%
  arrange(sample)
  
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

nir_irms_partial_final <- rbind(nir_irms_partial, missing_irms) %>%
  arrange(sample) %>%
  column_to_rownames("sample") %>%
  as.matrix()


##### Make Correlation Matrix ####
#For partial dataset, within species and location
#in which the indicator species was identified
#stats::cor(x,y, method="spearman")

indspec_full_mat <- indspec_full %>%
  semi_join(nir_irms_partial_final, by="sample") 

indspec_full_mat <- indspec_full_mat %>%
  arrange(sample) %>%
  column_to_rownames("sample") %>%
  select(-c(location, block, timepoint, species, treatment, intercrop)) %>%
  as.matrix()

#cor_partial <- cor(indspec_full_mat, nir_partial_mat, method="spearman")
cor_partial_p <- rcorr(indspec_full_mat, nir_irms_partial_final, type="spearman")


##### Visualize Matrices ####

#### 1. Partial Matrix ####
pal2 <- c("#A50026", "white","#313695")
r_mat <- cor_partial_p$r[25:30, 1:24]
p_mat <- cor_partial_p$P[25:30, 1:24]

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
  mutate(tax_otu_label = ifelse(#Maize
                                otu == "Otu00381", "Unclassified *Bacteria* (OTU 381)", #Burera
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

graph_cor_df$tax_otu_label <- factor(graph_cor_df$tax_otu_label, levels =c(
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


##### Partial correlation plot ####
RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
library(ggtext) #element_markdown

partial_cor_graph <- ggplot(graph_cor_df, aes(y=tax_otu_label, x=quality)) +
  geom_tile(aes(fill=spearman), color="lightgrey") +
  scale_fill_gradient2(low = "#fcfdbf", high= "#b73779", mid="lightgrey",
                       na.value = "#FFFFFF", name = "Spearman's\nCorrelation") +
  theme_bw() +
  scale_x_discrete(position ="top",labels = c(
            "N.corrected" = "N", "d15N.corrected" = expression(delta~""^15~"N"),
            "adf" = "ADF", "ndf" = "NDF", "lignin" = "Lignin", "protein" = "Protein")) +
  theme(legend.position="right",
        #plot.margin = margin(r=5, unit="cm"),
        axis.text.x = element_text(hjust=-0.01, angle=45, size=15, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size=15, face="bold"),
        axis.text.y = element_markdown(size=15)) 
partial_cor_graph + theme(axis.text.y = element_blank())

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Thesis/Rwanda/Lab work/04_Station_BNF/final_data/seq_data/16S")
ggsave("16S_partial_cor_graph_indspec.png", dpi=300, height = 14, width=5, units="in")

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
