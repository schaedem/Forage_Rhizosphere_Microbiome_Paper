library(car)
library(tidyverse)
library(lme4)
library(rebus)
library(broom.mixed)
library(broom)
library(broomExtra)
library(emmeans)

#Test whether relative abundance of AMF increases in legume intercropped
#plots at each location/timepoint

#read in guild data
guild <- readRDS("rare_its2_guild.rds")
str(guild)

guild2 <- guild %>%
  filter(confidence == "Probable" | confidence == "Highly Probable") %>%
  select(sample, otu, rel_abund, location, treatment, species,
         timepoint, block, trophic_mode, guild, confidence) %>%
  mutate(intercrop = ifelse(grepl(PLUS, treatment), "intercrop", "monocrop"))

#select symbiotrophs other than AMF
symbiotrophs <- guild2 %>%
  filter(grepl("ymbiotroph", trophic_mode)) %>%
  filter(!grepl("Arbuscular Mycorrhizal", guild))

#select white rot (not symbiotrophic)
white_rot <- guild %>%
  filter(grepl("White Rot", trait)) %>%
  filter(confidence == "Probable" | confidence == "Highly Probable") %>%
  mutate(intercrop = ifelse(grepl(PLUS, treatment), "intercrop", "monocrop")) %>%
  filter(!grepl("ymbiotroph", trophic_mode))

#select AMF
amf <- guild2 %>%
  filter(grepl("Arbuscular Mycorrhizal", guild))

#### prepare dataframes ####
#AMF
rubona_amf <- amf %>%
  filter(location == "Rubona")

burera_amf <- amf %>%
  filter(location == "Burera")

karama_amf <- amf %>%
  filter(location == "Karama")

#Other symbiotrophs
rubona_sym <- symbiotrophs %>%
  filter(location == "Rubona")

burera_sym <- symbiotrophs %>%
  filter(location == "Burera")

karama_sym <- symbiotrophs %>%
  filter(location == "Karama")

#White rot
rubona_wr <- white_rot %>%
  filter(location == "Rubona")

karama_wr <- white_rot %>%
  filter(location == "Karama")

burera_wr <- white_rot %>%
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
rubona_amf_results <- anova_test(rubona_amf)
r1_anova_amf <- rubona_amf_results[[1]] %>%
  mutate(location = "Rubona", timepoint = "Harvest 1")
r2_anova_amf <- rubona_amf_results[[2]]%>%
  mutate(location = "Rubona", timepoint = "Harvest 2")
r3_anova_amf <- rubona_amf_results[[3]]%>%
  mutate(location = "Rubona", timepoint = "Harvest 3")
r4_anova_amf <- rubona_amf_results[[4]]%>%
  mutate(location = "Rubona", timepoint = "Harvest 4")

rubona_amf_df <- rbind(r1_anova_amf, r2_anova_amf, r3_anova_amf, r4_anova_amf)

#Karama
karama_amf_results <- anova_test(karama_amf)
k1_anova_amf <- karama_amf_results[[1]] %>%
  mutate(location = "Karama", timepoint = "Harvest 1")
k2_anova_amf <- karama_amf_results[[2]]%>%
  mutate(location = "Karama", timepoint = "Harvest 2")
k3_anova_amf <- karama_amf_results[[3]]%>%
  mutate(location = "Karama", timepoint = "Harvest 3")
k4_anova_amf <- karama_amf_results[[4]]%>%
  mutate(location = "Karama", timepoint = "Harvest 4")

karama_amf_df <- rbind(k1_anova_amf, k2_anova_amf, k3_anova_amf, k4_anova_amf)

#Burera
burera_amf_results <- anova_test(burera_amf)
b1_anova_amf <- burera_amf_results[[1]] %>%
  mutate(location = "Burera", timepoint = "Harvest 1")
b2_anova_amf <- burera_amf_results[[2]]%>%
  mutate(location = "Burera", timepoint = "Harvest 2")
b3_anova_amf <- burera_amf_results[[3]]%>%
  mutate(location = "Burera", timepoint = "Harvest 3")

burera_amf_df <- rbind(b1_anova_amf, b2_anova_amf, b3_anova_amf)

graph_amf_letters <- rbind(rubona_amf_df,
                           karama_amf_df,
                           burera_amf_df)

#####Overall mixed-effects ANOVA for AMF ####
overall_amf_mod <- lmer(rel_abund~species -1 + (1|block:location) + (1|timepoint), data=amf)
summary(overall_amf_mod)
Anova(overall_amf_mod) #species is sig p=0.02
#Overall in Rubona
rubona_amf_mod <- lmer(rel_abund~species -1 + (1|block) + (1:timepoint), data=rubona_amf)
Anova(rubona_amf_mod) #species is is (p=0.03)
rubona_amf_emmeans <- emmeans(rubona_amf_mod, specs="species")
rubona_amf_letters <- multcomp::cld(rubona_amf_emmeans, adjust="Tukey", 
                                    Letters = letters) %>%
  mutate(location = "Rubona")

#Overall in Karama
karama_amf_mod <- lmer(rel_abund~species -1 + (1|block) + (1:timepoint), data=karama_amf)
Anova(karama_amf_mod) #species is sig (p<0.001)
karama_amf_emmeans <- emmeans(karama_amf_mod, specs="species")
karama_amf_letters <- multcomp::cld(karama_amf_emmeans, adjust="Tukey", 
                                    Letters = letters) %>%
  mutate(location = "Karama")

#Overall in Burera
burera_amf_mod <- lmer(rel_abund~species -1 + (1|block) + (1:timepoint), data=burera_amf)
Anova(burera_amf_mod) #species is ns (p=0.26)
burera_amf_emmeans <- emmeans(burera_amf_mod, specs="species")
burera_amf_letters <- multcomp::cld(burera_amf_emmeans, adjust="Tukey", Letters = letters) %>%
  mutate(location = "Burera")

graph_loc_overall <- rbind(burera_amf_letters, karama_amf_letters, rubona_amf_letters) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group, side="both"))

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
# maize and desmodium are different
r_int_amf_wilcox <- wilcox_test_abund(rubona_amf) %>% mutate(location = "Rubona")

#napier and maize are different
k_int_amf_wilcox <- wilcox_test_abund(karama_amf) %>% mutate(location = "Karama")

#maize marginal (p=0.07)
b_int_amf_wilcox <- wilcox_test_abund(burera_amf) %>% mutate(location="Burera")

all_wilcox <- rbind(r_int_amf_wilcox, k_int_amf_wilcox, b_int_amf_wilcox) %>%
  mutate(
    p_val = round(p_val, 2),
    p_val = ifelse(p_val > 0.10, NA, p_val))

##### Symbiotroph intercrop wilcox ####
r_int_sym_wilcox <- wilcox_test_abund(rubona_sym) %>% mutate(location = "Rubona")

k_int_sym_wilcox <- wilcox_test_abund(karama_sym) %>% mutate(location = "Karama")

b_int_sym_wilcox <- wilcox_test_abund(burera_sym) %>% mutate(location = "Burera")

all_wilcox_sym <- rbind(r_int_sym_wilcox, k_int_sym_wilcox, b_int_sym_wilcox) %>%
  mutate(
    p_val = round(p_val, 2),
    p_val = ifelse(p_val > 0.10, NA, p_val))

####White rot intercrop wilcox ####
r_int_wr_wilcox <- wilcox_test_abund(rubona_wr) %>% mutate(location = "Rubona")

k_int_wr_wilcox <- wilcox_test_abund(karama_wr) %>% mutate(location = "Karama")

b_int_wr_wilcox <- wilcox_test_abund(burera_wr) %>% mutate(location = "Burera")

all_wilcox_wr <- rbind(r_int_wr_wilcox, k_int_wr_wilcox, b_int_wr_wilcox) %>%
  mutate(
    p_val = round(p_val, 2),
    p_val = ifelse(p_val > 0.10, NA, p_val))

#### Graphing ####
#need to fix scale for tile plots so is uniform between plots
#tile plot for wilcox intercrop testing
wr_wilcox <- 
  ggplot(all_wilcox_wr, aes(x=species, y=timepoint, fill=p_val)) +
  geom_tile() +
  theme_test() +
  facet_wrap(~location) +
  coord_flip() +
  viridis::scale_fill_viridis(option="G",na.value = "white", name="Wilcoxon Rank Sum\np-value",
                              guide=guide_legend(direction="horizontal"), limits=c(0,0.1)) +
  xlab("") +
  ylab("Timepoint") +
  theme(legend.position = 'top', legend.title = element_text(size=9))
wr_wilcox
ggsave("wilcoxon_intercrop_white_rot_time_loc.png", dpi=300, width=12, height=5, unit="in")

amf_wilcox <- 
  ggplot(all_wilcox, aes(x=species, y=timepoint, fill=p_val)) +
  geom_tile() +
  theme_test() +
  facet_wrap(~location) +
  coord_flip() +
  viridis::scale_fill_viridis(option="G",na.value = "white", name="Wilcoxon Rank Sum\np-value",
                       guide=guide_legend(direction="horizontal"), limits=c(0,0.1)) +
  xlab("") +
  ylab("Timepoint") +
  theme(legend.position = 'top', legend.title = element_text(size=9))
amf_wilcox

ggsave("wilcoxon_intercrop_amf_time_loc.png", dpi=300, width=12, height=5, unit="in")

symb_wilcox <-
  ggplot(all_wilcox_sym, aes(x=species, y=timepoint, fill=p_val)) +
  geom_tile() +
  theme_test() +
  facet_wrap(~location) +
  coord_flip() +
  viridis::scale_fill_viridis(option="G",na.value = "white", 
                      name="Wilcoxon Rank Sum\np-value", labels=c("<0.01", "<0.02",
                                                                  "<0.05", "<0.08", "<0.10"),
                       guide=guide_legend(direction="horizontal"), limits=c(0,0.1)) +
  xlab("") +
  ylab("Timepoint") +
  theme(legend.position = 'top', legend.title = element_text(size=9))
symb_wilcox
ggsave("wilcoxon_intercrop_sym_time_loc.png", dpi=300, width=12, height=5, unit="in")

#putting tile plots together
library(cowplot)
legend <- get_legend(symb_wilcox)
both_tile <- plot_grid(amf_wilcox + theme(legend.position = "none") + ylab(""), 
                       symb_wilcox + theme(legend.position = "none") + ylab(""),
                       wr_wilcox + theme(legend.position = "none"), 
                       align="h",
                       ncol=1, labels="AUTO", label_size=12)
both_tile
tile_with_legend <- plot_grid(legend, 
                              both_tile, 
                              ncol=1, rel_heights = c(0.25,5))
tile_with_legend

ggsave("amf_sym_wr_wilcox_tile.png", dpi=300, height=10, width=10, unit="in", bg="white")

#Relative abundance plots
pal <- c("#1F78B4", "#6A3D9A", "#33A02C", "#E31A1C")
amf <- amf %>%
  mutate(timepoint = case_when(timepoint == 1 ~ "Harvest 1",
                               timepoint == 2 ~ "Harvest 2",
                               timepoint == 3 ~ "Harvest 3",
                               timepoint == 4 ~ "Harvest 4"))

#graphing rel abundance of AMF by species and location
ggplot(amf, aes(x=species, y=100*(rel_abund), color=species)) +
  #geom_boxplot() +
  theme_test() +
  theme(legend.position = "none") +
  geom_jitter(stat="identity", alpha=0.3) +
  facet_grid(~location) +
  xlab("") +
  ylab("Relative Abundance (%)") +
  scale_color_manual(values=pal) +
  scale_y_continuous(limits=c(-0.05, 2.5), expand=c(0,0)) +
  geom_text(data=graph_loc_overall, aes(x=species, y=2.4, label=.group), inherit.aes=FALSE)
ggsave("amf_rel_abund_loc.png", dpi=300, height=6, width=8.8, units="in")


#graphing rel abundance of AMF by species and location x timepoint
ggplot(amf, aes(x=species, y=100*rel_abund, color=species)) +
  #geom_boxplot() +
  theme_test() +
  theme(legend.position = "none") +
  geom_jitter(stat="identity", alpha=0.3) +
  facet_grid(location~timepoint) +
  xlab("") +
  ylab("Relative Abundance (%)") +
  scale_color_manual(values=pal) +
  scale_y_continuous(limits=c(-0.05, 2.5), expand=c(0,0))+
  geom_text(data=graph_amf_letters, aes(x=species, y=2.4, label=.group), 
            inherit.aes=FALSE)
ggsave("amf_rel_abund_loc_time.png", dpi=300, height=8, width=14, units="in")

#graphing rel abundance of symbiotrophs by species and location x timepoint
ggplot(symbiotrophs, aes(x=species, y=100*rel_abund, color=species)) +
  #geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none") +
  geom_jitter(stat="identity", alpha=0.3) +
  facet_wrap(~location) +
  xlab("") +
  ylab("Relative Abundance (%)") +
  scale_color_manual(values=pal)
ggsave("symbiotrophs_rel_abund_loc.png", dpi=300, height=8, width=14, units="in")
