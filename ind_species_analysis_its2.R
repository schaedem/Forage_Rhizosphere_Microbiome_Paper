library(tidyverse)
library(ggVennDiagram)
library(rebus) #constants for special characters
library(vegan)
library(broom)
library(ggtext)
library(ggpubr) #for adding brackets to ggplots

#read in taxonomy file
tax <- read_rds("tidy_taxonomy_ITS2.rds")
# glimpse(tax)
# 
# #read in indicator results files
# karama <- read_csv("karama_its_indicator_species.csv")
# rubona <- read_csv("rubona_its_indicator_species.csv")
# burera <- read_csv("burera_its2_indicator_species.csv")
# 
# indspec <- rbind(karama, rubona, burera) %>%
#   mutate(p.value.bh = stats::p.adjust(p.value, method="BH")) #Benjamini & Hochberg p-value correction; doesn't change results

# indspec_tax <- inner_join(indspec, tax, by="otu") %>%
#   filter(s.yes == 1) #filtering for otus that are positively associated with intercropping

# unclassified <- indspec_tax %>%
#   filter(grepl("unclassified", order) | grepl("incertae_sedis", order))

# write_csv(indspec_tax, "indicator_species_taxonomy.csv")
# write_csv(unclassified, "indicator_its2_species_unclassified.csv")

#### Get trophic and guild classifications ####
indspec <- read_csv("indicator_species_taxonomy_its2.csv")
class <- read_rds("rare_its2_guild.rds") %>%
  dplyr::select(otu, trophic_mode, guild, growth, trait, confidence, notes) %>%
  group_by(otu) %>%
  unique()
#glimpse(class)

#Join indspec and funguild
indspec_class <- inner_join(indspec, class, by="otu") %>%
  filter(grepl("Probable", confidence) | grepl("Highly Probable", confidence))
#all associated with brachiaria are unassigned.
test <- indspec_class[3,]
test$notes

#Use Wilcox test to validate sig differences using relative abundance
data <- read_rds("full_rare_rel_abund_ITS2.rds")
#glimpse(data)
data <- data %>%
  rename(otu = name) %>%
  mutate(treatment=as.character(treatment),
         intercrop = ifelse(grepl(PLUS, treatment), "yes", "no")) 

sum_data <- data %>%
  inner_join(indspec, by=c("otu", "location", "species"))

sum_data %>%
  filter(rel_abund > 0) %>%
  summarize(av_abund = mean(rel_abund) *100)

##separate OTUs by location
burera <- data %>% filter(location=="Burera")

rubona <- data %>% filter(location=="Rubona")

karama <- data %>% filter(location=="Karama")
which(karama$sample == "222R")

#Join location dfs to indicator species classification and perform Wilcox test
#Karama
k_indspec <- inner_join(karama, indspec_class, by=c("otu", "location", "species")) %>%
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
r_indspec <- inner_join(rubona, indspec_class, by=c("otu", "location", "species")) %>%
  group_by(species, otu) %>%
  nest() %>%
  mutate(test = map(.x=data, ~wilcox.test(rel_abund~intercrop, data=.x) %>% tidy)) %>%
  unnest() %>%
  rename(p.value.wilcox = p.value1) %>%
  mutate(p.wilcox.adjust = stats::p.adjust(p.value.wilcox, method="BH")) %>%#adjustment does not change p-value
  select(otu, p.wilcox.adjust) %>%
  unique() %>%
  mutate(location = "Rubona")

#Burera
b_indspec <- inner_join(burera, indspec_class, by=c("otu", "location", "species")) %>%
  group_by(species, otu) %>%
  nest() %>%
  mutate(test = map(.x=data, ~wilcox.test(rel_abund~intercrop, data=.x) %>% tidy)) %>%
  unnest() %>%
  rename(p.value.wilcox = p.value1) %>%
  mutate(p.wilcox.adjust = stats::p.adjust(p.value.wilcox, method="BH")) %>%#adjustment does not change p-value
  select(otu, p.wilcox.adjust) %>%
  unique() %>%
  mutate(location = "Burera")

#join wilcox-tested OTUs for all locations and merge back with taxonomic and rel-abund data
indspec_wilcox <- rbind(k_indspec, r_indspec, b_indspec) %>%
  inner_join(indspec_class, by=c("location", "otu", "species")) %>%
  select(otu, species, location,  kingdom, phylum, class, order, 
         family, genus,trophic_mode, guild, confidence)

#clean up ITS2 dataframe and transform rel abund to % rel abund
#need to add a tiny fraction to % rel abund to avoid issues with taking the log
data <- data %>%
  select(sample, otu, rel_abund, location, treatment, species, timepoint, block, intercrop) %>%
  mutate(rel_abund_trans = 100*(rel_abund + 1/4000))

indspec_join <- inner_join(indspec_wilcox, data, by=c("location", "otu", "species"))

#finally, join with taxonomy information
indspec_join <- indspec_join %>%
  inner_join(tax, by=c("otu", "kingdom", "phylum", "class", "order", "family", "genus")) %>%
  mutate(location2 = as.numeric(ifelse(location == "Burera", 1, 
                            ifelse(location == "Karama", 2,
                                   ifelse(location == "Rubona", 3, location2)))))

#separate by location
# b_join <- indspec_join %>% filter(location=="Burera")
# k_join <- indspec_join %>% filter(location=="Karama")
# r_join <- indspec_join %>% filter(location == "Rubona")

#separate by species and format/reorder OTU labels
####Maize formatting ####
maize <- indspec_join %>% filter(species=="Maize") %>%
  mutate(tax_otu_label = ifelse(grepl("OTU 322", tax_otu_label), "Unclassified *Fungi* (OTU 322)",
                        ifelse(grepl("OTU 234", tax_otu_label), "Unclassified *Fungi* (OTU 234)", 
                        ifelse(grepl("OTU 22", tax_otu_label), "Unclassified *Fungi* (OTU 22)", 
                        ifelse(grepl("OTU 44", tax_otu_label), "*Neoroussoella* (OTU 44)",
                       ifelse(grepl("OTU 43", tax_otu_label), "*Acrocalymma* (OTU 43)", tax_otu_label)))))) 
maize$tax_otu_label <- factor(maize$tax_otu_label, levels=c("Unclassified *Fungi* (OTU 22)",
                                                        "Unclassified *Fungi* (OTU 322)",
                                                        "*Acrocalymma* (OTU 43)",
                                                        "*Neoroussoella* (OTU 44)",
                                                        "Unclassified *Fungi* (OTU 234)"))

       #  tax_otu_label = fct_reorder(tax_otu_label, location2, .desc=TRUE)) #order by location for graphing

####Napier formatting ####
napier <- indspec_join %>% filter(species=="Napier") %>%
  mutate(tax_otu_label = ifelse(grepl("OTU 43", tax_otu_label), "*Acrocalymma* (OTU 43)",
                          ifelse(grepl("OTU 266", tax_otu_label), "Unclassified *Fungi* (OTU 266)",
                          ifelse(grepl("OTU 133", tax_otu_label), "*Penicillifer* (OTU 133)",
                          ifelse(grepl("OTU 200", tax_otu_label), "*Microascus* (OTU 200)",
                          ifelse(grepl("OTU 276", tax_otu_label), "*Chaetomium* (OTU 276)",
                          ifelse(grepl("OTU 499", tax_otu_label), "Unclassified *Phaeosphaeriaceae*<br> (OTU 499)",
                          ifelse(grepl("OTU 828", tax_otu_label), "*Dentiscutata* (OTU 828)", tax_otu_label))))))))

napier$tax_otu_label <- factor(napier$tax_otu_label, 
                               levels=c("*Penicillifer* (OTU 133)", "*Microascus* (OTU 200)", "*Chaetomium* (OTU 276)", 
                                        "Unclassified *Phaeosphaeriaceae*<br> (OTU 499)", "*Dentiscutata* (OTU 828)",
                                        "*Acrocalymma* (OTU 43)", "Unclassified *Fungi* (OTU 266)"))


####Brachiaria formatting ####
brach <- indspec_join %>% filter(species=="Brachiaria") %>%
  mutate(tax_otu_label = ifelse(grepl("OTU 73", tax_otu_label), "*Chordomyces* (OTU 73)",
                         ifelse(grepl("OTU 729", tax_otu_label), "Unclassified *Hypocreales*<br>(OTU 729)", 
                                      tax_otu_label)))

####Desmodium formatting ####
desmodium <- indspec_join %>% filter(species=="Desmodium") %>%
  mutate(tax_otu_label = ifelse(grepl("OTU 106", tax_otu_label), "Unclassified *Orbiliales*<br>(OTU 106)",
                          ifelse(grepl("OTU 865", tax_otu_label), "Unclassified *Chytridiomycota*<br>(OTU 865)",
                          ifelse(grepl("OTU 20", tax_otu_label), "*Pyrenochaetopsis* (OTU 20)",
                          ifelse(grepl("OTU 47", tax_otu_label), "*Striaticonidium* (OTU 47)",
                          ifelse(grepl("OTU 148", tax_otu_label), "*Pyrenochaetopsis* (OTU 148)",
                          ifelse(grepl("OTU 86", tax_otu_label), "*Thermomyces* (OTU 86)", tax_otu_label)))))))
  
desmodium$tax_otu_label <- factor(desmodium$tax_otu_label, 
                                  levels=c("*Thermomyces* (OTU 86)", #Burera
                                  "Unclassified *Orbiliales*<br>(OTU 106)", "Unclassified *Chytridiomycota*<br>(OTU 865)", #Karama
                                  "*Pyrenochaetopsis* (OTU 20)","*Striaticonidium* (OTU 47)", "*Pyrenochaetopsis* (OTU 148)")) #Rubona

# lev <- unique(desmodium$tax_otu_label)
# lev
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
  scale_y_log10() +
  facet_wrap(~species) +
  theme_classic() +
  xlab("") +
  ylab("Relative abundance (%)") +
  coord_flip() +
  geom_hline(yintercept = 20/2097, size=0.5, color="gray") + #detection limit
  theme(legend.title = element_blank(),
        strip.text = element_text(size=15),
        axis.title.x = element_text(size=12),
        axis.text.y = element_markdown(size=10),
        legend.position = "top") +
  scale_color_manual(values = c("grey","#33A02C"), labels=c("Monocropped", "Intercropped")) +
  scale_fill_manual(values = c("grey","#33A02C"), labels=c("Monocropped", "Intercropped")) +
  geom_bracket(xmin=c(0.8,2.8), xmax=c(2.4,5.4),
               y.position=log10(max(maize$rel_abund_trans)+8),
               label=c("Burera", "Karama"), coord.flip=TRUE, inherit.aes=FALSE, label.size=3.5) +
  annotate("text", x=c(1,2,3,4,5), y=5, 
           label = "italic(p<0.01)", parse=TRUE, size=3.5, alpha=0.85)
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
  geom_hline(yintercept = 20/2097, size=0.5, color="gray") + #detection limit
  theme(legend.title = element_blank(),
        strip.text = element_text(size=15),
        axis.title.x = element_text(size=12),
        axis.text.y = element_markdown(size=10),
        legend.position = "top") +
 scale_color_manual(values = c("grey","#E31A1C"), labels=c("Monocropped", "Intercropped")) +
 scale_fill_manual(values = c("grey","#E31A1C"), labels=c("Monocropped", "Intercropped")) +
  geom_bracket(xmin=c(0.8,5.8), xmax=c(5.4,7.4),
               y.position=log10(max(napier$rel_abund_trans)+8),
               label=c("Rubona", "Karama"), coord.flip=TRUE, inherit.aes=FALSE, label.size=3.5) +
  annotate("text", x=c(1,2,3,4,5, 6, 7), y=3,
           label = "italic(p<0.01)", parse=TRUE, size=3.5, alpha=0.85)

napier_plot  

RColorBrewer::display.brewer.pal(10, "Paired")
RColorBrewer::brewer.pal(10, "Paired")

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
  geom_hline(yintercept = 20/2097, size=0.5, color="gray") + #detection limit
  theme(legend.title = element_blank(),
        strip.text = element_text(size=15),
        axis.title.x = element_text(size=12),
        axis.text.y = element_markdown(size=10),
        legend.position = "top") +
  scale_color_manual(values = c("grey","#1F78B4"), labels=c("Monocropped", "Intercropped")) +
  scale_fill_manual(values = c("grey","#1F78B4"), labels=c("Monocropped", "Intercropped")) +
  geom_bracket(xmin=c(0.8), xmax=c(2.4),
               y.position=log10(max(brach$rel_abund_trans)+8),
               label=c("Karama"), coord.flip=TRUE, inherit.aes=FALSE, label.size=3.5) +
  annotate("text", x=c(1,2), y=3,
           label = "italic(p<0.01)", parse=TRUE, size=3.5, alpha=0.85)
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
  geom_hline(yintercept = 20/2097, size=0.5, color="gray") + #detection limit
  theme(legend.title = element_blank(),
        strip.text = element_text(size=15),
        axis.title.x = element_text(size=12),
        axis.text.y = element_markdown(size=10),
        legend.position = "top") +
  scale_color_manual(values = c("grey","#6A3D9A"), labels=c("Monocropped", "Intercropped")) +
  scale_fill_manual(values = c("grey","#6A3D9A"), labels=c("Monocropped", "Intercropped")) +
  geom_bracket(xmin=c(0.8, 1.8, 3.8), xmax=c(1.4, 3.4, 6.4),
               y.position=log10(max(brach$rel_abund_trans)+8),
               label=c("Burera", "Karama", "Rubona"), coord.flip=TRUE, inherit.aes=FALSE, label.size=3.5) +
  annotate("text", x=c(1:6), y=3,
           label = "italic(p<0.01)", parse=TRUE, size=3.5, alpha=0.85)
desmodium_plot


#### Assembling plots ####
library(gridExtra)
library(grid)
library(cowplot)

maize_nap <-
          plot_grid(maize_plot + ylab("") + theme(legend.position="none"),
          napier_plot + theme(legend.position="none"), align="v", ncol=1,
          labels=c("A","C"))

brach_des <-
          plot_grid(brach_plot+ylab("")+theme(legend.position="none"),
                    desmodium_plot + theme(legend.position="none"), align="v", ncol=1,
                    labels=c("B", "D"))

final <- plot_grid(maize_nap, brach_des, ncol=2, align=c("v", "h"), rel_widths=c(1.05,1))
final

legend <- get_legend(desmodium_plot +
                       scale_color_manual(values = c("grey","#4C4E52"), labels=c("Monocropped", "Intercropped")) +
                       scale_fill_manual(values = c("grey","#4C4E52"), labels=c("Monocropped", "Intercropped")) +
                       theme(legend.position="bottom",
                             legend.text=element_text(size=10)))
final_with_legend <- plot_grid(
  final, legend, ncol=1, rel_heights=c(4,0.3)
)

final_with_legend
ggsave("its2_indspec_high_res.tiff", final_with_legend, 
       dpi=300, width=12.5, height=8, units=c("in"), bg="white")
