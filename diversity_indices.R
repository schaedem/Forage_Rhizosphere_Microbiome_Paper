library(tidyverse)
library(vegan)
set.seed(10011992)

shared <- readRDS("final_shared_rare.rds")
print(shared, n=1)


meta <- read_csv("final_samps.csv") %>%
  mutate(treatment = as.character(treatment),
         treatment = str_replace_all(treatment, "�", ""),
         species = as.factor(species),
         treatment = ifelse(grepl("Brachiaria v. Mulato II", treatment), "Brachiaria cv. Mulato II",
                            ifelse(grepl("Maize", treatment), "Maize monocrop", treatment)),
         #treatment = as.factor(treatment),
         intercrop = ifelse(grepl(" + ", treatment), "intercrop", "monocrop")) 

head(meta)

full_shared <- shared %>%
  inner_join(meta, by="sample")

str(full_shared)
saveRDS(full_shared, file="full_shared_rare_16S.rds")
head(full_shared_alpha)

full_shared_alpha <- full_shared %>%
  group_by(sample) %>%
  summarize(
    sobs = specnumber(value),
    shannon=diversity(value, index="shannon"),
    simpson=diversity(value, index="simpson"),
    invsimpson=1/simpson, 
    n=sum(value)
  ) %>%
  ungroup() %>%
  merge(meta, by="sample")

ggplot(full_shared_alpha, aes(x=species, y=shannon, group=interaction(timepoint, species), fill=species)) +
  geom_boxplot() +
  facet_wrap(~location) +
  theme_classic() +
  xlab("")

#richness = total number per taxa

#Shannon diversity index; more balanced accounting of diversity

#Simpson and inverse Simpson

#Beta diversity

#Bray-curtis
head(full_shared)

#All locations
all_shared <- shared %>%
  pivot_wider(names_from=name, values_from=value, values_fill=0) %>%
  column_to_rownames("sample")

#all_16S_dist <- vegdist(all_shared, method="bray")

saveRDS(all_16S_dist, file="all_16S_dist.rds")

#Burera only df

burera_shared <- full_shared %>%
  filter(location=="Burera") %>%
  dplyr::select(-c(location, treatment, species, timepoint)) %>%
  pivot_wider(names_from=name, values_from=value, values_fill=0) %>%
  column_to_rownames("sample")

#burera_16S_dist <- vegdist(burera_shared, method="bray")

saveRDS(burera_16S_dist, file="burera_16S_dist.rds")

#Karama only df

karama_shared <- full_shared %>%
  filter(location=="Karama") %>%
  dplyr::select(-c(location, treatment, species, timepoint)) %>%
  pivot_wider(names_from=name, values_from=value, values_fill=0) %>%
  column_to_rownames("sample")
head(karama_shared)

#karama_16S_dist <- vegdist(karama_shared, method="bray")

saveRDS(karama_16S_dist, file="karama_16S_dist.rds")


#Rubona only df

rubona_shared <- full_shared %>%
  filter(location=="Rubona")%>%
  dplyr::select(-c(location, treatment, species, timepoint)) %>%
  pivot_wider(names_from=name, values_from=value, values_fill=0) %>%
  column_to_rownames("sample")

#rubona_16S_dist <- vegdist(rubona_shared, method="bray")

saveRDS(rubona_16S_dist, file="rubona_16S_dist.rds")
