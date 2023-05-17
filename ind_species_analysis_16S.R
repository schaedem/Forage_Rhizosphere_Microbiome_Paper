library(tidyverse)
library(rebus) #constants for special characters
library(vegan)
library(broom)
library(ggtext)
library(ggpubr) #for adding brackets to ggplots
library(cowplot)


#read in taxonomy file
tax <- read_rds("tidy_taxonomy_16S.rds")
# glimpse(tax)
# # 
# #read in indicator results files
karama <- read_csv("karama_indicator_species_2.csv")
rubona <- read_csv("rubona_indicator_species_2.csv")
burera <- read_csv("burera_indicator_spcecies_2.csv")
#
indspec <- rbind(karama, rubona, burera) %>%
  mutate(p.value.bh = stats::p.adjust(p.value, method="BH")) #Benjamini & Hochberg p-value correction; doesn't change results

indspec_tax <- inner_join(indspec, tax, by="otu") %>%
  filter(s.yes == 1) #filtering for otus that are positively associated with intercropping

unclassified <- indspec_tax %>%
  filter(grepl("unclassified", order) | grepl("classified", family) |
           grepl("classified", kingdom) | grepl("classified", phylum))

write_csv(indspec_tax, "indicator_species_taxonomy_2.csv")
write_csv(unclassified, "indicator_16S_species_unclassified_2.csv")

#### Get trophic and N function classifications ####
indspec_tax <- read_csv("indicator_species_taxonomy_2.csv")
class <- read_rds("16S_trophic_class.rds")
glimpse(class)

indspec_class <- inner_join(indspec_tax, class, tax, by="otu")

#Use Wilcox test to validate sig differences using relative abundance
data <- read_rds("otu_rare_relabund_16S.rds")
#glimpse(data)
data <- data %>%
  mutate(treatment=as.character(treatment),
         intercrop = ifelse(grepl(PLUS, treatment), "yes", "no")) 

sum_data <- inner_join(data, indspec_class, by=c("otu", "location", "species"))

sum_data %>%
  filter(rel_abund > 0) %>%
  summarize(av_abund = mean(rel_abund) * 100)

##separate OTUs by location
burera_otu <- data %>% filter(location=="Burera")

rubona_otu <- data %>% filter(location=="Rubona")

karama_otu <- data %>% filter(location=="Karama")
#which(karama_otu$sample == "222R")

#Join location dfs to indicator species classification and perform Wilcox test
#Karama
k_indspec <- inner_join(karama_otu, indspec_class, by=c("otu", "location", "species")) %>%
  group_by(species, otu) %>%
  nest() %>%
  mutate(test = map(.x=data, ~wilcox.test(rel_abund~intercrop, data=.x) %>% tidy)) %>%
  unnest() %>%
  rename(p.value.wilcox = p.value1) %>%
  mutate(p.wilcox.adjust = stats::p.adjust(p.value.wilcox, method="BH")) %>%#adjustment does not change p-value
  dplyr::select(otu, p.wilcox.adjust) %>%
  unique() %>%
  mutate(location = "Karama")

#Rubona
r_indspec <- inner_join(rubona_otu, indspec_class, by=c("otu", "location", "species")) %>%
  group_by(species, otu) %>%
  nest() %>%
  mutate(test = map(.x=data, ~wilcox.test(rel_abund~intercrop, data=.x) %>% tidy)) %>%
  unnest() %>%
  rename(p.value.wilcox = p.value1) %>%
  mutate(p.wilcox.adjust = stats::p.adjust(p.value.wilcox, method="BH")) %>%#adjustment does not change p-value
  dplyr::select(otu, p.wilcox.adjust) %>%
  unique() %>%
  mutate(location = "Rubona")

#Burera
b_indspec <- inner_join(burera_otu, indspec_class, by=c("otu", "location", "species")) %>%
  group_by(species, otu) %>%
  nest() %>%
  mutate(test = map(.x=data, ~wilcox.test(rel_abund~intercrop, data=.x) %>% tidy)) %>%
  unnest() %>%
  rename(p.value.wilcox = p.value1) %>%
  mutate(p.wilcox.adjust = stats::p.adjust(p.value.wilcox, method="BH")) %>%#adjustment does not change p-value
  dplyr::select(otu, p.wilcox.adjust) %>%
  unique() %>%
  mutate(location = "Burera")

#join wilcox-tested OTUs for all locations and merge back with taxonomic and rel-abund data
#24 OTUs
indspec_wilcox <- rbind(k_indspec, r_indspec, b_indspec) %>%
  inner_join(indspec_class, by=c("location", "otu", "species")) %>%
  dplyr::select(otu, species, location, kingdom, phylum, class, order, 
         family, genus, trophicMode, n_function)

#clean up 16S dataframe and transform rel abund to % rel abund
#need to add a tiny fraction to % rel abund to avoid issues with taking the log
data <- data %>%
  dplyr::select(sample, otu, rel_abund, location, treatment, species, timepoint, block, intercrop) %>%
  mutate(rel_abund_trans = 100*(rel_abund + 1/4000))

check <- data %>%
  group_by(sample) %>%
  summarize(N = sum(rel_abund))

indspec_join <- inner_join(indspec_wilcox, data, by=c("location", "otu", "species"))

#separate by species and format/reorder OTU labels
####Maize formatting ####
maize <- indspec_join %>% filter(species=="Maize") %>%
  mutate(tax_otu_label = ifelse(otu == "Otu00381", "Unclassified *Bacteria* (OTU 381)", #Burera
                        ifelse(otu == "Otu00571", "Unclassified *Bacteria* (OTU 571)", #Rubona
                        ifelse(otu == "Otu01065", "Unclassified *Sphingomonadaceae*<br>(OTU 1065)", #Rubona
                        NA))))
maize$tax_otu_label <- factor(maize$tax_otu_label,
                              levels = c("Unclassified *Bacteria* (OTU 381)",
                                         "Unclassified *Bacteria* (OTU 571)",
                                         "Unclassified *Sphingomonadaceae*<br>(OTU 1065)"))  
  

#   mutate(tax_otu_label = ifelse(otu == "Otu00759", "Unclassified *Deltaproteobacteria*<br>(OTU 759)",
#                          ifelse(otu == "Otu01236", "*Acidobacteria* Gp1 (OTU 1236)", #Burera
#                          ifelse(otu == "Otu00104", "*Acidobacteria* Gp3 (OTU 104)", #Karama
#                          ifelse(otu == "Otu00341", "*Singulisphaera* (OTU 341)", #Karama; degrades biopolymers under acidic conditions
#                          ifelse(otu == "Otu00468", "*Micromonospora* (OTU 468)", #Karama; spore-forming saprotrophs with antibiotic and antifungal compounds
#                          ifelse(otu == "Otu00827", "Unclassified *Bacillales*<br>(OTU 827)", #K; beneficial endophytes
#                          ifelse(otu == "Otu00934", "Unclassified *Rhizobiales*<br>(OTU 934)", #K
#                          ifelse(otu == "Otu01356", "*Acidobacteria* Gp5 (OTU 1356)", #K
#                          ifelse(otu == "Otu01535", "Unclassified *Bacteria* (OTU 1535)", #K
#                          ifelse(otu == "Otu01648", "Unlassified *Bacteria* (OTU 1648)", #k
#                          ifelse(otu == "Otu00363", "Unclassified *Rhizobiales*<br>(OTU 363)", #Rubona; beneficial endophytes
#                          ifelse(otu == "Otu00576", "Unclassified *Actinobacteria*<br>(OTU 576)", #Rubona; decompose organic matter, form nets similar to fungi
#                          ifelse(otu == "Otu01065", "Unclassified *Sphingomonadaceae*<br>(OTU 1065)", #Rubona; mostly heterotrophic
#                          ifelse(otu == "Otu01202", "Unclassified *Acetobacteraceae*<br>(OTU 1202)", #Rubona; acidophilic
#                          ifelse(otu == "Otu01857", "*Flavisolibacter* (OTU 1857)", NA)))))))))))))))) #Rubona; some can degrade chitin or cellulose
# 
# maize$tax_otu_label <- factor(maize$tax_otu_label, levels=c(
#   "Unclassified *Deltaproteobacteria*<br>(OTU 759)","*Acidobacteria* Gp1 (OTU 1236)",
#   "*Acidobacteria* Gp3 (OTU 104)","*Singulisphaera* (OTU 341)","*Micromonospora* (OTU 468)",
#   "Unclassified *Bacillales*<br>(OTU 827)","Unclassified *Rhizobiales*<br>(OTU 934)",
#   "*Acidobacteria* Gp5 (OTU 1356)","Unclassified *Bacteria* (OTU 1535)","Unlassified *Bacteria* (OTU 1648)",
#   "Unclassified *Rhizobiales*<br>(OTU 363)","Unclassified *Actinobacteria*<br>(OTU 576)",
#   "Unclassified *Sphingomonadaceae*<br>(OTU 1065)","Unclassified *Acetobacteraceae*<br>(OTU 1202)",
#   "*Flavisolibacter* (OTU 1857)"
# ))

####Napier formatting####
napier <- indspec_join %>% filter(species=="Napier") %>%
  mutate(tax_otu_label = ifelse(otu == "Otu01452", "*Phenylobacterium* (OTU 1452)", #Burera; aerobic soil-associated
                                ifelse(otu == "Otu01953", "*Gemmatimonas* (OTU 1953)", #gram-negative, non-spore forming, potentially beneficial (Li et a, 2022)
                                ifelse(otu == "Otu00316", "Unclassified *Micrococcaceae*<br>(OTU 316)", #Karama
                                ifelse(otu == "Otu00597", "*Verrucomicrobia* Subdivision 3<br>(OTU 597)",
                                ifelse(otu == "Otu01297", "Unclassified *Betaproteobacteria*<br>(OTU 1297)",
                                ifelse(otu == "Otu01653", "Unclassified *Bacteria* (OTU 1653)", #Karama
                                ifelse(otu == "Otu00516", "*Gemmata* (OTU 516)", #Rubona; aerobic chemoheterotrophs
                                ifelse(otu == "Otu00690", "Unclassified *Bacteria* (OTU 690)",
                                ifelse(otu == "Otu01116", "Unclassified *Planctomycetaceae*<br>(OTU 1116)",
                                 ifelse(otu == "Otu01444", "Unclassified *Actinomycetales*<br>(OTU 1444)",
                                        NA)))))))))))
napier$tax_otu_label <- factor(napier$tax_otu_label, levels=c(
  "*Phenylobacterium* (OTU 1452)", "*Gemmatimonas* (OTU 1953)", 
  "Unclassified *Micrococcaceae*<br>(OTU 316)", "*Verrucomicrobia* Subdivision 3<br>(OTU 597)",
  "Unclassified *Betaproteobacteria*<br>(OTU 1297)","Unclassified *Bacteria* (OTU 1653)",
  "*Gemmata* (OTU 516)", "Unclassified *Bacteria* (OTU 690)", 
  "Unclassified *Planctomycetaceae*<br>(OTU 1116)","Unclassified *Actinomycetales*<br>(OTU 1444)"
))

#   mutate(tax_otu_label = ifelse(otu == "Otu00662", "*Aciditerrimonas* (OTU 662)", #Burera; acidophilic, thermophilic Actinobacteria
#                          ifelse(otu == "Otu01631", "Unclassified *Bacteria* (OTU 1631)", #burera
#                          ifelse(otu == "Otu01972", "Unclassified *Bacteria* (OTU 1972)", #burera
#                          ifelse(otu == "Otu00987", "Unclassified *Rhizobiales*<br>(OTU 987)", #Karama
#                          ifelse(otu == "Otu01916", "Unclassified *Rhodospirillales*<br>(OTU 1916)", #Karama
#                          ifelse(otu == "Otu00185", "Unclassified *Actinomycetales*<br>(OTU 185)", #Rubona; facultatively anaerobic with branching growth; antimicrobial properties
#                          ifelse(otu == "Otu00643", "*Paenibacillus* (OTU 634)", #R; facultative anaerobes associated with rhizosphere
#                          ifelse(otu == "Otu00690", "Unclassified *Bacteria* (OTU 690)", #R
#                          ifelse(otu == "Otu01057", "*Singulisphaera* (OTU 1057)", #R; acidophilic genus
#                          ifelse(otu == "Otu01830", "Unclassified *Verrucomicrobia* Subdivision 3<br>(OTU 1830)", #R
#                          ifelse(otu == "Otu02018", "*Acidobacteria* Gp1 (OTU 2018)", NA)))))))))))) #R
# 
# napier$tax_otu_label <- factor(napier$tax_otu_label, levels=c(
#   "*Aciditerrimonas* (OTU 662)", "Unclassified *Bacteria* (OTU 1631)","Unclassified *Bacteria* (OTU 1972)",
#   "Unclassified *Rhizobiales*<br>(OTU 987)","Unclassified *Rhodospirillales*<br>(OTU 1916)",
#   "Unclassified *Actinomycetales*<br>(OTU 185)", "*Paenibacillus* (OTU 634)",
#   "Unclassified *Bacteria* (OTU 690)", "*Singulisphaera* (OTU 1057)",
#   "Unclassified *Verrucomicrobia* Subdivision 3<br>(OTU 1830)",
#   "*Acidobacteria* Gp1 (OTU 2018)"
# ))

####Brachiaria formatting ####
brach <- indspec_join %>% filter(species=="Brachiaria") %>%
  mutate(tax_otu_label = ifelse(otu == "Otu00584", "*Acidobacteria* Gp6 (OTU 584)", #burera
                         ifelse(otu == "Otu01574", "*Paenibacillus* (OTU 1574)", #Burera
                         ifelse(otu == "Otu02740", "Unclassified *Ktedonobacteria*<br>(OTU 2740)",
                         ifelse(otu == "Otu03228", "Unclassified *Bacteria* (OTU 3228)",
                         ifelse(otu == "Otu03700", "*Singulisphaera* (OTU 3700)",
                         ifelse(otu == "Otu00138", "*Conexibacter* (OTU 138)",#karama; gram-pos enriched in maize rhizosphere (Deng et al, 2021)
                         ifelse(otu == "Otu00364", "Unclassified *Bacteria* (OTU 364)",
                         ifelse(otu == "Otu00861", "Unclassified *Clostridiales*<br>(OTU 861)", #Rubona
                         ifelse(otu == "Otu01256", "*Aciditerrimonas* (OTU 1256)", #rubona; Fe-reducing thermoacidiphilic actinobacterium
                         ifelse(otu == "Otu01388", "Unclassified *Bacteria* (OTU 1388)", NA)))))))))))

brach$tax_otu_label <- factor(brach$tax_otu_label, levels = c(
  "*Acidobacteria* Gp6 (OTU 584)", "*Paenibacillus* (OTU 1574)",
  "Unclassified *Ktedonobacteria*<br>(OTU 2740)","Unclassified *Bacteria* (OTU 3228)",
  "*Singulisphaera* (OTU 3700)","*Conexibacter* (OTU 138)", "Unclassified *Bacteria* (OTU 364)",
  "Unclassified *Clostridiales*<br>(OTU 861)", "*Aciditerrimonas* (OTU 1256)","Unclassified *Bacteria* (OTU 1388)"
))

#   mutate(tax_otu_label = ifelse(otu == "Otu02595", "Unclassified *Bacteria* (OTU 2595)", #Burera
#                          ifelse(otu == "Otu00381", "Unclassified *Bacteria* (OTU 381)", #Karama
#                          ifelse(otu == "Otu00569", "Unclassified *Bacteria* (OTU 569)",#k
#                          ifelse(otu == "Otu01608", "Unclassified *Bacteria* (OTU 1608)",#k
#                          ifelse(otu == "Otu00415", "Unclassified *Proteobacteria* (OTU 415)", #Rubona
#                          ifelse(otu == "Otu01622", "Unclassified *Bacteria* (OTU 1622)", NA)))))))
# 
# brach$tax_otu_label <- factor(brach$tax_otu_label, levels=c(
#   "Unclassified *Bacteria* (OTU 2595)", "Unclassified *Bacteria* (OTU 381)",
#   "Unclassified *Bacteria* (OTU 569)", "Unclassified *Bacteria* (OTU 1608)",
#   "Unclassified *Proteobacteria* (OTU 415)","Unclassified *Bacteria* (OTU 1622)"
# ))

####Desmodium formatting ####
desmodium <- filter(indspec_join, species=="Desmodium") %>%
  mutate(tax_otu_label = ifelse(otu == "Otu00472", "*Hyphomicrobium* (OTU 472)", NA)) #denitrifier (Yang & Crowley, 2000)
desmodium$tax_otu_label <- factor(desmodium$tax_otu_label)


#   mutate(tax_otu_label = ifelse(otu == "Otu00615", "*Thermopolyspora* (OTU 615)", #Burera
#                          ifelse(otu == "Otu01053", "Unclassified *Bacteria* (OTU 1053)",#Rubona
#                          NA)))
# desmodium$tax_otu_label <- factor(desmodium$tax_otu_label, levels=c(
#   "*Thermopolyspora* (OTU 615)", "Unclassified *Bacteria* (OTU 1053)"
# ))

#### Plotting ####

#### Plots by species #####

#### Maize plot ####
maize_plot <- 
  ggplot(aes(y=rel_abund_trans, x=tax_otu_label, color=intercrop, fill=intercrop), data=maize) +
  geom_jitter(position=position_jitterdodge(dodge.width = 0.8, 
                                            jitter.width = 0.3), shape=20, size=3, alpha=0.8) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.8),
               color="black", show.legend = FALSE) +
  facet_wrap(~species) +
  theme_classic() +
  xlab("") +
  ylab("Relative abundance (%)") +
  coord_flip() +
  geom_hline(yintercept = 21/2097, size=0.5, color="gray") + #detection limit
  theme(legend.title = element_blank(),
        strip.text = element_text(size=15),
        axis.title.x = element_text(size=12),
        axis.text.y = element_markdown(size=10),
        legend.position = "top") +
  scale_color_manual(values = c("grey","#33A02C"), labels=c("Monocropped", "Intercropped")) +
  scale_fill_manual(values = c("grey","#33A02C"), labels=c("Monocropped", "Intercropped")) +
  geom_bracket(xmin = c(0.8, 1.8), xmax = c(1.4, 3.4), 
               y.position=max(maize$rel_abund_trans), label = c("Burera", "Rubona"), 
               coord.flip=TRUE, inherit.aes=FALSE, label.size=3.5) +
  annotate("text", x=c(1:3), y=0.6,
           label="p<0.01", fontface='italic', size=3.5, alpha=0.85) +
   scale_y_log10() 

maize_plot  



#### Napier plot ####
napier_plot <- ggplot(aes(y=rel_abund_trans, x=tax_otu_label, 
                          color=intercrop, fill=intercrop), data=napier) +
  geom_jitter(position=position_jitterdodge(dodge.width = 0.8, 
                                            jitter.width = 0.3), shape=20, size=3, alpha=0.8) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.8),
               color="black", show.legend = FALSE) +
  scale_y_log10() +
  facet_wrap(~species) +
  theme_classic() +
  xlab("") +
  ylab("Relative abundance (%)") +
  coord_flip() +
  geom_hline(yintercept = 21/2097, size=0.5, color="gray") + #detection limit
  theme(legend.title = element_blank(),
        strip.text = element_text(size=15),
        axis.title.x = element_text(size=12),
        axis.text.y = element_markdown(size=10),
        legend.position = "top") +
  scale_color_manual(values = c("grey","#E31A1C"), labels=c("Monocropped", "Intercropped")) +
  scale_fill_manual(values = c("grey","#E31A1C"), labels=c("Monocropped", "Intercropped")) +
  geom_bracket(xmin=c(0.8, 2.8, 6.8), xmax=c(2.4,6.4, 10.4), 
               y.position = log10(max(napier$rel_abund_trans)+0.7),
               label=c("Burera", "Karama", "Rubona"), coord.flip=TRUE, inherit.aes=FALSE, label.size=3.5) +
  annotate("text", x=c(1:10), y=0.6,
           label = "p<0.01", fontface='italic',size=3.5, alpha=0.85)

napier_plot

#### Brachiaria plot ####
brach_plot <- 
  ggplot(aes(y=rel_abund_trans, x=tax_otu_label, 
             color=intercrop, fill=intercrop), data=brach) +
  geom_jitter(position=position_jitterdodge(dodge.width = 0.8, 
                                            jitter.width = 0.3), shape=20, size=3, alpha=0.8) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.8),
               color="black", show.legend = FALSE) +
  scale_y_log10() +
  facet_wrap(~species) +
  theme_classic() +
  xlab("") +
  ylab("Relative abundance (%)") +
  coord_flip() +
  geom_hline(yintercept = 21/2097, size=0.5, color="gray") + #detection limit
  theme(legend.title = element_blank(),
        strip.text = element_text(size=15),
        axis.title.x = element_text(size=12),
        axis.text.y = element_markdown(size=10),
        legend.position = "top") +
  scale_color_manual(values = c("grey","#1F78B4"), labels=c("Monocropped", "Intercropped")) +
  scale_fill_manual(values = c("grey","#1F78B4"), labels=c("Monocropped", "Intercropped")) +
  geom_bracket(xmin = c(0.8, 5.8, 7.8), xmax=c(5.4, 7.4, 10.4),
               y.position=log10(max(brach$rel_abund_trans)+1.5),
               label=c("Burera", "Karama", "Rubona"), coord.flip = TRUE, inherit.aes = FALSE, label.size=3.5) +
  annotate("text", x=c(1:3,5:8,10), y=0.8, label="p<0.01", fontface='italic', size=3.5, alpha=0.85) +
  annotate("text", x=c(4,9), y=0.8, label="p=0.01", fontface='italic', size=3.5, alpha=0.85)



  # geom_bracket(xmin=c(0.8, 1.8, 4.8), xmax=c(1.4, 4.4,6.4),
  #              y.position=log10(max(brach$rel_abund_trans)+1.5),
  #              label=c("Burera", "Karama", "Rubona"), coord.flip=TRUE, inherit.aes=FALSE, label.size=3.5) +
  # annotate("text", x=c(3,4), y=0.8, label="p=0.01", fontface='italic', 
  #          size=3.5, alpha=0.85) +
  # annotate("text", x=c(1:2, 5:6), y=0.8, label="p<0.01", fontface='italic',
  #          size=3.5, alpha=0.85) +
  # annotate("rect", alpha=0.1, fill="grey", xmin=1.6, xmax=4.45, ymin=0, ymax=1.8)

brach_plot  

#### Desmodium plot ####
desmodium_plot <- 
  ggplot(aes(y=rel_abund_trans, x=tax_otu_label, 
             color=intercrop, fill=intercrop), data=desmodium) +
  geom_jitter(position=position_jitterdodge(dodge.width = 0.8, 
                                            jitter.width = 0.3), shape=20, size=3, alpha=0.8) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.8),
               color="black", show.legend = FALSE) +
  scale_y_log10() +
  facet_wrap(~species) +
  theme_classic() +
  xlab("") +
  ylab("Relative abundance (%)") +
  coord_flip() +
  geom_hline(yintercept = 21/2097, size=0.5, color="gray") + #detection limit
  theme(legend.title = element_blank(),
        strip.text = element_text(size=15),
        axis.title.x = element_text(size=12),
        axis.text.y = element_markdown(size=10),
        legend.position = "top") +
  scale_color_manual(values = c("grey","#6A3D9A"), labels=c("Monocropped", "Intercropped")) +
  scale_fill_manual(values = c("grey","#6A3D9A"), labels=c("Monocropped", "Intercropped")) +
  geom_bracket(xmin=c(0.8), xmax=c(1.4),
               y.position=log10(max(brach$rel_abund_trans)+1),
               label=c("Burera"), coord.flip=TRUE, inherit.aes=FALSE, label.size=3.5) +
  annotate("text", x=c(1), y=0.6,
           label = "p<0.01", fontface='italic', size=3.5, alpha=0.85)
desmodium_plot

#### Assembling plots ####

maize_nap <-
  plot_grid(maize_plot + ylab("") + theme(legend.position="none"),
            napier_plot + theme(legend.position="none"), align="v", ncol=1,
            labels=c("A","C"))
maize_nap

brach_des <-
  plot_grid(brach_plot+ylab("")+theme(legend.position="none"),
            desmodium_plot + theme(legend.position="none"), align="v", ncol=1,
            labels=c("B", "D"))

brach_des

final <- plot_grid(maize_nap, brach_des, ncol=2, align=c("v", "h"), rel_widths=c(1,1))
#final

legend <- get_legend(desmodium_plot +
                       scale_color_manual(values = c("grey","#4C4E52"), labels=c("Monocropped", "Intercropped")) +
                       scale_fill_manual(values = c("grey","#4C4E52"), labels=c("Monocropped", "Intercropped")) +
                       theme(legend.position="bottom",
                             legend.text=element_text(size=10)))
final_with_legend <- plot_grid(
  final, legend, ncol=1, rel_heights=c(4,0.3)
)

final_with_legend

ggsave("16S_indspec_high_res.tiff", final_with_legend, 
       dpi=300, width=14, height=10, units=c("in"), bg="white")

#any shared OTUs?
shared_indspec <- indspec_wilcox %>%
  group_by(otu) %>%
  nest() %>%
  mutate(num = map_dbl(data, nrow)) %>%
  filter(num > 1)
#none shared
