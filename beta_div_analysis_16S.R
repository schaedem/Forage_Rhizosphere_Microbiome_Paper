library(tidyverse)
library(vegan)
library(rebus) #for special characters
library(cowplot)
library(lme4)
library(emmeans) #means separation
library(multcomp) #for compact letter display

meta <- read_csv("final_samps.csv") %>%
  mutate(treatment = as.character(treatment),
         treatment = str_replace_all(treatment, "�", ""),
         species = as.factor(species),
         treatment = ifelse(grepl("Brachiaria v. Mulato II", treatment), "Brachiaria cv. Mulato II",
                            ifelse(grepl("Maize", treatment), "Maize monocrop", treatment)),
         treatment = as.character(treatment),
         intercrop = ifelse(grepl(PLUS, treatment), "yes", "no"))

make_dist <- function(filename, meta) {
  meta <- meta
  
  dist <- readRDS(filename)%>%
    as.matrix() %>%
    as_tibble(rownames = "sample") %>%
    pivot_longer(-sample) %>%
    filter(sample<name) %>%
    inner_join(meta, by="sample") %>%
    dplyr::select(sample, name, timepoint, intercrop, species, value, block) %>% #name is comparison sample
    rename(sample_a = sample, 
           intercrop_a = intercrop,
           species_a = species,
           timepoint_a = timepoint,
           block_a = block,
           sample = name) %>%
    inner_join(meta, by="sample") %>% #getting metadata for comparison sample
    dplyr::select(-c(treatment, block, location)) %>%
    rename(sample_b = sample,
           intercrop_b = intercrop,
           timepoint_b = timepoint,
           species_b = species) %>%
    filter(species_a == species_b) %>% #to compare within species
    filter(intercrop_a < intercrop_b) %>% #to compare between intercrop
    mutate(diff = abs(timepoint_a - timepoint_b))
  return(dist)
  
}

burera_dist <- make_dist("burera_16S_dist.rds", meta)
karama_dist <- make_dist("karama_16S_dist.rds", meta)
rubona_dist <- make_dist("rubona_16S_dist.rds", meta)  

karama_median <- karama_dist %>%
  group_by(sample_a, species_a) %>%
  summarize(median = median(value))

rubona_median <- rubona_dist %>%
  group_by(sample_a, species_a) %>%
  summarize(median = median(value))

burera_median <- burera_dist %>%
  group_by(sample_a, species_a) %>%
  summarize(median = median(value))

####test for sig differences among species ####
k_beta_mod <- lmer(value ~ species_a + (1|timepoint_a) + (1|block_a), data=karama_dist)
summary(k_beta_mod)
car::Anova(k_beta_mod) #species is significant (p<0.001)

r_beta_mod <- lmer(value ~ species_a + (1|timepoint_a) + (1|block_a), data=rubona_dist)
summary(r_beta_mod)
car::Anova(r_beta_mod) #species is significant (p<0.001)

b_beta_mod <- lmer(value ~ species_a + (1|timepoint_a) + (1|block_a), data=burera_dist)
summary(b_beta_mod)
car::Anova(b_beta_mod) #species is significant (p<0.001)

#Means separation:
#Maize had the greatest difference in beta-diversity between int and mono in Burera, Karama

#Karama
k_means <- emmeans(k_beta_mod, list(pairwise~species_a), adjust="tukey")$`pairwise differences of species_a` %>%
  as_tibble() %>%
  filter(p.value < 0.06)
k_means #brach and nap are not diff; marginal between des and nap; all others are diff

#getting letters for graphing
k_spec_means <- emmeans(k_beta_mod, specs = "species_a")
k_spec_means
k_letters <- cld(k_spec_means, adjust="Tukey", Letters=letters, alpha=0.05)
k_letters <- k_letters %>% as_tibble() %>%
  mutate(.group = str_trim(.group, side="both"))

#Rubona
r_means <- emmeans(r_beta_mod, list(pairwise~species_a), adjust="tukey")$`pairwise differences of species_a` %>%
  as_tibble() %>%
  filter(p.value < 0.06)
r_means #des/maize, brach/des, and maize/nap are sig

#letters for graphing
r_spec_means <- emmeans(r_beta_mod, specs = "species_a")
r_spec_means
r_letters <- cld(r_spec_means, adjust="Tukey", Letters=letters, alpha=0.05)
r_letters <- r_letters %>% as_tibble() %>%
  mutate(.group = str_trim(.group, side="both"))

#Burera
b_means <- emmeans(b_beta_mod, list(pairwise~species_a), adjust="tukey")$`pairwise differences of species_a` %>%
  as_tibble() %>%
  filter(p.value < 0.06)
b_means #brach/maize, des/maize, and maize/nap are sig

#letters for graphing
b_spec_means <- emmeans(b_beta_mod, specs = "species_a")
b_spec_means
b_letters <- cld(b_spec_means, adjust="Tukey", Letters=letters, alpha=0.05)
b_letters <- b_letters %>% as_tibble() %>%
  mutate(.group = str_trim(.group, side="both"))
b_letters


#### Plotting results ####
#custom color palette for plotting
pal <- c("#1F78B4", "#6A3D9A", "#33A02C", "#E31A1C")

karama_bray_plot <- 
  ggplot(aes(x=species_a, y=value, color=species_a), data=karama_dist) +
  geom_boxplot() +
  theme_classic() +
  geom_jitter(aes(x=species_a, y=value, color=species_a, group=paste0(sample_a,species_a)),
              inherit.aes=FALSE, alpha=0.3) +
  ylab("Bray-Curtis Distance") +
  xlab("") +
  theme(legend.position = "none") +
  scale_color_manual(values = pal) +
  scale_y_continuous(limits=c(0.4,0.8)) +
  geom_text(data=k_letters, aes(label=`.group`, x=species_a), y=0.8, inherit.aes=FALSE)
karama_bray_plot #maize trending higher; letters are same for median vs all distances

rubona_bray_plot <- 
  ggplot(aes(x=species_a, y=value, color=species_a), data=rubona_dist) +
  geom_boxplot() +
  theme_classic() +
  geom_jitter(aes(x=species_a, y=value, color=species_a, group=paste0(sample_a,species_a)),
              inherit.aes=FALSE, alpha=0.3) +
  ylab("Bray-Curtis Distance") +
  xlab("") +
  theme(legend.position = "none")+
  scale_color_manual(values = pal)+
  scale_y_continuous(limits=c(0.4,0.8)) +
  geom_text(data=r_letters, aes(label=`.group`, x=species_a), 
            y=0.8, inherit.aes=FALSE)
rubona_bray_plot

burera_bray_plot <- 
  ggplot(aes(x=species_a, y=value, color=species_a), data=burera_dist) +
  geom_boxplot() +
  theme_classic() +
  geom_jitter(aes(x=species_a, y=value, color=species_a, group=paste0(sample_a,species_a)),
              inherit.aes=FALSE, alpha=0.3) +
  ylab("Bray-Curtis Distance") +
  xlab("") +
  theme(legend.position = "none")+
  scale_color_manual(values = pal) +
  scale_y_continuous(limits=c(0.4,0.8)) +
  geom_text(data=b_letters, aes(label=`.group`, x=species_a),
            y=0.8, inherit.aes=FALSE)
burera_bray_plot #maize trending higher

plot_grid(burera_bray_plot, 
          karama_bray_plot, 
          rubona_bray_plot, labels="AUTO", align="v", label_size = 10)

ggsave("bray_curtis_diff_intercrop.tiff", dpi=300, height = 8, width=10, units="in", bg="white")

