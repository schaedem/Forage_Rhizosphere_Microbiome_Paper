library(tidyverse)
library(lme4)
library(emmeans)

data <- read_csv("final_good_leg_ndfa.csv")  %>%
  dplyr::select(sample, ndfa, rep, location, treatment, timepoint, intercrop) %>%
  mutate(date = case_when(timepoint ==1 ~ "3/1/2021",
                          timepoint ==2 ~ "6/1/2021", 
                          timepoint ==3 ~ "9/1/2021",
                          timepoint ==4 ~ "12/1/2021"),
         date=lubridate::mdy(date))
head(data)

##### Overall models(all locations/timepoints) #####

mod_all <- lmer(ndfa ~ intercrop + (1|location:rep) + (1|timepoint), data=data)
summary(mod_all)
car::Anova(mod_all) #intercrop highly sig (p<0.001)

mod_all.2 <- lmer(ndfa ~ intercrop + treatment + (1|location:rep) + (1|timepoint), data=data)
car::Anova(mod_all.2) #both intercrop and treatment are highly sig (p<0.001)

mod_all.3 <- lmer(ndfa ~ treatment + (1|location:rep) + (1|timepoint), data=data)
car::Anova(mod_all.3) #treatment highly sig (p<0.001)

anova(mod_all, mod_all.2, mod_all.3) #mod_all.2 or mod_all.3 are preferred

mod_all.4 <- lmer(ndfa ~ location + (1|rep) + (1|timepoint), data=data) 
car::Anova(mod_all.4) #location is significant

mod_all.5 <- lmer(ndfa~as.factor(timepoint) + (1|location:rep), data=data)
car::Anova(mod_all.5) #sampling timepoint is significant

#overall, how does the non-legume species affect Ndfa in Desmodium?
post_hoc <- emmeans(mod_all.3, list(pairwise~treatment), adjust="tukey")
multcomp::cld(post_hoc)
#how does Ndfa vary by location?
post_hoc_loc <- emmeans(mod_all.4, list(pairwise~location), adjust="tukey")

#how does Ndfa vary by time?
sum_df <- data %>%
  group_by(location, date) %>%
  summarise(ndfa = mean(ndfa))

time_plot <- ggplot(sum_df, aes(x=date, y=ndfa, color=location, group=interaction(location))) +
  #geom_boxplot() +
  geom_line() +
  geom_point(stat="identity") +
  theme_classic() +
  ylab("Ndfa (%)") +
  xlab("")
time_plot
#appears higher in timepoint 3 

##### By location #####
karama <- data %>% filter(location == "Karama") #45 obs
rubona <- data %>% filter(location == "Rubona") #51 obs
burera <- data %>% filter(location == "Burera") #27 obs

#Karama
k.mod <- lmer(ndfa~intercrop + treatment + (1|rep) + (1|timepoint), data=karama)
summary(k.mod)
car::Anova(k.mod) #both are sig (p<0.001)
k_post_hoc <- emmeans(k.mod, list(pairwise~treatment), adjust="tukey")
k_cld <- multcomp::cld(k_post_hoc) %>% mutate(location = "Karama")
# treatment              intercrop emmean    SE   df lower.CL upper.CL
# Desmodium              no          33.2  8.97 4.95    10.07     56.3
# Desmodium + Brachiaria yes         40.9  9.48 5.94    17.62     64.1
# Desmodium + maize      yes         27.1 10.03 7.29     3.58     50.6
# Desmodium + Napier     yes         85.5  9.52 5.99    62.15    108.8 *

#Rubona
r.mod <- lmer(ndfa~intercrop + treatment + (1|rep) + (1|timepoint), data=rubona)
summary(r.mod)
car::Anova(r.mod) #both are sig (p<0.001)
r_post_hoc <- emmeans(r.mod, list(pairwise~treatment), adjust="tukey")
r_cld <- multcomp::cld(r_post_hoc) %>% mutate(location = "Rubona")
# treatment              intercrop emmean   SE   df lower.CL upper.CL
# Desmodium              no          44.1 16.0 4.34    0.958     87.3
# Desmodium + Brachiaria yes         91.2 15.4 3.72   47.143    135.2*
# Desmodium + maize      yes         36.8 15.7 4.00   -6.750     80.4
# Desmodium + Napier     yes        101.9 16.0 4.32   58.691    145.1*

#Burera
b.mod <- lmer(ndfa~intercrop + treatment + (1|rep) + (1|timepoint), data=burera)
summary(b.mod)
car::Anova(b.mod) #both are sig (p=0.01)
b_post_hoc <- emmeans(b.mod, list(pairwise~treatment), adjust="tukey")
b_cld <- multcomp::cld(b_post_hoc) %>% mutate(location = "Burera")


### graphing ####

pal <- c("#6A3D9A", "#1F78B4", "#B2DF8A", "#FB9A99")

#collate cld info
cld_all <- rbind(k_cld, r_cld, b_cld) %>%
  mutate(.group = as.character(.group),
         .group = str_trim(.group, side = "both"),
         label = case_when(
           .group == "1" ~ "a",
           .group == "2" ~ "b",
           .group == "12" ~ "ab"
         ),
         treatment = as.factor(treatment))

data <- data %>%
  mutate(treatment = as.factor(treatment))

loc_box_plot <- ggplot(data, aes(x=treatment, y=ndfa, color=treatment)) +
  geom_boxplot() +
  geom_jitter(alpha=0.8) +
  scale_color_manual(values = pal) +
  theme_test() +
  facet_wrap(~location) +
  theme(axis.text.x = element_blank(),
        axis.ticks.length.x = unit(0, units="cm"),
        legend.title = element_blank(),
        ) +
  xlab("")+
  ylab("Ndfa (%)") +
  geom_text(data=cld_all, aes(x=treatment, label=label, y=165),
            color="black")

loc_box_plot
ggsave("ndfa_loc_boxplot.png", dpi=300, width = 8, 
       height = 4, units="in", bg="white")
