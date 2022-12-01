library(tidyverse)
library(vegan)
library(RColorBrewer)
library(ggtext)

#Meta data
meta <- read_csv("final_samps.csv") %>%
  mutate(treatment = as.character(treatment),
         treatment = str_replace_all(treatment, "�", ""),
         species = as.factor(species),
         treatment = ifelse(grepl("Brachiaria v. Mulato II", treatment), "Brachiaria cv. Mulato II",
                            ifelse(grepl("Maize", treatment), "Maize monocrop", treatment)),
         treatment = as.factor(treatment)) 

#Taxonomy file
# taxonomy <- read_tsv("final.opti_mcc.0.03.cons.taxonomy") %>%
#   dplyr::select("OTU", "Taxonomy") %>%
#   rename_all(tolower) 
# head(taxonomy)
# 
# taxonomy2 <- taxonomy %>%
#   mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
#          taxonomy = str_replace(taxonomy, ";$", "")) %>%
#   separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"),
#            sep=";") %>%
#   mutate(phylum = str_replace_all(phylum, '"', ""))
# head(taxonomy2)
# saveRDS(taxonomy2, "tidy_taxonomy.rds")
taxonomy2 <- readRDS("tidy_taxonomy.rds")

#OTU file
shared_full <- readRDS("final_shared_rare.rds")

# shared_join <- shared_full %>%
#   rename(otu = name) %>%
#   inner_join(taxonomy2, by="otu") %>%
#   group_by(sample) %>%
#   mutate(rel_abund = value / sum(value)) %>% #relative abundance
#   ungroup()

#saveRDS(shared_join, file="shared_rare_taxon.rds")
shared_join <- readRDS("shared_rare_taxon.rds")
head(shared_join)

otu_rel_abund <- shared_join %>%
 # select(-c(class, order, family, genus)) %>%
  group_by(sample, phylum) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  ungroup() 
head(otu_rel_abund)

mean_rel_abund <- otu_rel_abund %>%
  inner_join(meta, by="sample") %>%
  group_by(phylum, location) %>%
  summarize(mean_rel_abund = 100*(mean(rel_abund)), .groups="drop") %>%
  mutate(phylum = str_replace(phylum,
                              "(.*)_unclassified", "Unclassified *\\1*"),
         phylum = str_replace(phylum,
                              "(.*)_Chloroplast", "*\\1* / Chloroplast"),
         phylum=str_replace(phylum,
                            "^(\\S*)$", "*\\1*"))

#look at which groups are the most abundant
#for now, will use Bacteroidetes as cutoff (max > 3.78)
mean_rel_abund %>%
  group_by(phylum) %>%
  summarize(max = max(mean_rel_abund)) %>%
  arrange(desc(max))

pool_rel_abund <- mean_rel_abund %>%
  group_by(phylum) %>%
  summarize(pool = max(mean_rel_abund) < 3.7, .groups="drop", 
            mean = mean(mean_rel_abund))
head(pool_rel_abund)

graph_rel_abund <- inner_join(mean_rel_abund, pool_rel_abund, by="phylum") %>%
  mutate(phylum = if_else(pool, "Other", phylum)) %>%
  group_by(location, phylum) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), 
            mean = sum(mean),
            .groups="drop") %>%
  mutate(phylum = factor(phylum),
         phylum = fct_reorder(phylum, mean, .desc=TRUE),
         phylum = fct_shift(phylum, n=1))
head(graph_rel_abund)

##### stacked barplots by phylum ####
loc_bar <- 
  ggplot(graph_rel_abund, aes(x=location, y=mean_rel_abund, fill=phylum)) +
  geom_col() +
  theme_classic() +
  labs(y="Mean Relative Abundance (%)") +
  scale_fill_manual(name=NULL,
                        breaks = c("*Acidobacteria*", "*Actinobacteria*",
                                   "*Bacteroidetes*", "*Chloroflexi*", 
                                   "*Firmicutes*", "*Planctomycetes*",
                                   "*Proteobacteria*", "*Verrucomicrobia*", 
                                   "Unclassified *Bacteria*", "Other"),
                    values = c(brewer.pal(9, "Paired"), "grey")) +
  theme(legend.key.size=unit(10, "pt"),
        legend.text=element_markdown(),
        legend.title = element_blank()) +
  xlab("") +
  scale_y_continuous(expand = c(0,0))
#  scale_fill_brewer(palette = "Paired")
loc_bar

