library(tidyverse)
library(car)
library(lme4)
library(rebus)
library(broom.mixed)
library(broom)
library(broomExtra)
library(emmeans)

#metadata
meta <- read_csv("final_samps.csv")

#functional assignments
class <- readRDS("16S_trophic_class.rds") 

#relative abundance from rarefied dataset
rel_abund <- readRDS("otu_rare_relabund_16S.rds")
head(rel_abund)

#rarefied counts
shared_new <- rel_abund %>%
  inner_join(class, by="otu") %>%
  filter(!is.na(n_function)) 

#Test whether relative abundance of endophytic N-fixers increases in legume intercropped
#plots at each location/timepoint
endophytes <- shared_new %>%
  filter(grepl("endophyte", n_function)) %>%
  mutate(intercrop = ifelse(grepl(PLUS, treatment), "intercrop", "monocrop"))

endophyte_summary <- endophytes %>%
  group_by(location, genus) %>%
  summarize(rel_abund = mean(rel_abund)*100)

rubona <- endophytes %>%
  filter(location == "Rubona")

karama <- endophytes %>%
  filter(location == "Karama")

burera <- endophytes %>%
  filter(location == "Burera")


####ANOVA to test abundance by species within location x timepoint###
anova_test <- function(loc_df) {
  data <- loc_df
  
  model_function <- function(data) {
    fit <- lmer(rel_abund~species - 1 + (1|block), data = data)
    out <- emmeans(fit, specs="species")
    letters <- multcomp::cld(out, Letters=letters, alpha=0.059, adjust="Tukey") %>% 
      as_tibble() %>%
      mutate(.group = str_trim(.group, side="both"))
    # contrasts <- out$`pairwise differences of species` %>% as_tibble() %>%
    #   rename(contrasts = `1`)
  }
  
  data2 <- data %>%
    group_by(timepoint) %>%
    group_map(~model_function(.x))
  
  return(data2)
  
}

####within loc x time for AMF ####
#Rubona
rubona_results <- anova_test(rubona)
r1_anova <- rubona_results[[1]] %>%
  mutate(location = "Rubona", timepoint = "Harvest 1")
r2_anova <- rubona_results[[2]]%>%
  mutate(location = "Rubona", timepoint = "Harvest 2")
r3_anova <- rubona_results[[3]]%>%
  mutate(location = "Rubona", timepoint = "Harvest 3")
r4_anova <- rubona_results[[4]]%>%
  mutate(location = "Rubona", timepoint = "Harvest 4")

rubona_df <- rbind(r1_anova, r2_anova, r3_anova, r4_anova) # no differences

#Karama
karama_results <- anova_test(karama)
k1_anova <- karama_results[[1]] %>%
  mutate(location = "Karama", timepoint = "Harvest 1")
k2_anova <- karama_results[[2]]%>%
  mutate(location = "Karama", timepoint = "Harvest 2")
k3_anova <- karama_results[[3]]%>%
  mutate(location = "Karama", timepoint = "Harvest 3")
k4_anova <- karama_results[[4]]%>%
  mutate(location = "Karama", timepoint = "Harvest 4")

karama_df <- rbind(k1_anova, k2_anova, k3_anova, k4_anova) #no differences

#Burera
burera_results <- anova_test(burera)
b1_anova<- burera_results[[1]] %>%
  mutate(location = "Burera", timepoint = "Harvest 1")
b2_anova <- burera_results[[2]]%>%
  mutate(location = "Burera", timepoint = "Harvest 2")
b3_anova <- burera_results[[3]]%>%
  mutate(location = "Burera", timepoint = "Harvest 3")

burera_df <- rbind(b1_anova, b2_anova, b3_anova) #no differences

graph_amf_letters <- rbind(rubona_amf_df,
                           karama_amf_df,
                           burera_amf_df)

#####Overall mixed-effects ANOVA for endophytic N fixers ####
overall_mod <- lmer(rel_abund~species -1 + (1|block:location) + (1|timepoint), data=endophytes)
Anova(overall_mod) #species is sig p<0.001
emmeans(overall_mod, specs = "species", adjust="Tukey")

#Overall in Rubona
rubona_mod <- lmer(rel_abund~species -1 + (1|block) + (1:timepoint), data=rubona)
Anova(rubona_mod) #species is ns (0.98)
rubona_emmeans <- emmeans(rubona_mod, specs="species")
rubona_letters <- multcomp::cld(rubona_emmeans, adjust="Tukey", 
                                    Letters = letters) %>%
  mutate(location = "Rubona")

#Overall in Karama
karama_mod <- lmer(rel_abund~species -1 + (1|block) + (1:timepoint), data=karama)
Anova(karama_mod) #species is ns (p=0.69)
karama_emmeans <- emmeans(karama_mod, specs="species")

#Overall in Burera
burera_mod <- lmer(rel_abund~species -1 + (1|block) + (1:timepoint), data=burera)
Anova(burera_mod) #species is ns (p=0.59)
burera_emmeans <- emmeans(burera_mod, specs="species")

####Wilcoxon to test abundance by intercrop within location x timepoint x species####

wilcox_test_abund <- function(loc_df){
  data <- loc_df
  
  wilcox_function <- function(data) {
    fit <- wilcox.test(rel_abund~intercrop, data=data) 
  }
  
  data2 <- data %>%
    group_by(timepoint, species) %>%
    nest() %>%
    mutate(wilcox_int = map(data, ~wilcox_function(.)))
  
  data3 <- data2 
  
  p_vals <- numeric()
  
  for (i in 1:nrow(data3)) {
    p_vals <- c(p_vals, pluck(data3, 4, i, "p.value"))
  }
  
  data3$p_val <- p_vals
  
  data3 <- data3 %>%
    dplyr::select(species, timepoint, p_val)
  
  return(data3)
  
}

####AMF intercrop wilcox ####
r_int_wilcox <- wilcox_test_abund(rubona) %>% mutate(location = "Rubona")

k_int_wilcox <- wilcox_test_abund(karama) %>% mutate(location = "Karama")

b_int_wilcox <- wilcox_test_abund(burera) %>% mutate(location="Burera")

all_wilcox <- rbind(r_int_wilcox, k_int_wilcox, b_int_wilcox) %>%
  mutate(
    p_val = round(p_val, 2),
    p_val = ifelse(p_val > 0.10, NA, p_val))

