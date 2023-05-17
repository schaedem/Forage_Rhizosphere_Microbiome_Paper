library(tidyverse)
library(lme4)
library(emmeans)
library(cowplot)
source("ndfa_lmer.R")

data <- read_csv("final_good_leg_ndfa.csv")  %>%
  dplyr::select(sample, ndfa, rep, location, treatment, timepoint, intercrop) %>%
  mutate(date = case_when(timepoint ==1 ~ "3/1/2021",
                          timepoint ==2 ~ "6/1/2021", 
                          timepoint ==3 ~ "9/1/2021",
                          timepoint ==4 ~ "12/1/2021"),
         date=lubridate::mdy(date),
         timepoint = ifelse(timepoint=="1", "Harvest 1",
                     ifelse(timepoint=="2", "Harvest 2",
                    ifelse(timepoint == "3", "Harvest 3",
                    ifelse(timepoint == "4", "Harvest 4", timepoint)))))
head(data)

#read in soil data
soil <- read_csv("final_minN_for_analysis.csv")

#merge soil and ndfa data
data_full <- inner_join(data, soil, by=c("location", "rep", "treatment", "timepoint"))

karama <- data_full %>% filter(location=="Karama")
rubona <- data_full %>% filter(location == "Rubona")
burera <- data_full %>% filter(location == "Burera")

#calculate correlations overall (Spearmans)

overall_cor_minN <- cor.test(x=data_full$minN, y=data_full$ndfa,
                             method="spearman")
#rho = -0.395, p<0.001

overall_cor_gwc <- cor.test(x=data_full$gwc, y=data_full$ndfa,
                            method="spearman")
#rho = -0.214, p=0.02

overall_cor_ph <- cor.test(x=data_full$pH, y=data_full$ndfa,
                            method="spearman")
#rho = 0.06, p=0.48

#correlations within location (Spearmans)

#graphing
pal <- c("#6A3D9A", "#1F78B4", "#B2DF8A", "#FB9A99")
loc_pal <- c("#FC8D62","#E78AC3","#E5C494")

minN_plot <- ggplot(data_full, aes(x=minN, y=ndfa, color=location)) +
  geom_point() +
  geom_smooth(method="gam", se=F, color="black") +
  theme_test() +
  ylab("Ndfa (%)") +
  theme(legend.title = element_blank()) +
  xlab(expression("minN (mg"~kg^-1~ "soil)")) +
  #facet_wrap(~location) #+
  scale_color_manual(values = loc_pal)

gwc_plot <- ggplot(data_full, aes(x=100*gwc, y=ndfa, color=location)) +
  geom_point() +
  geom_smooth(method="gam", se=F, color="black") +
  theme_test() +
  ylab("Ndfa (%)") +
  theme(legend.title = element_blank()) +
  xlab(expression("Gravimetric water content (%)")) +
  scale_color_manual(values = loc_pal) 

#arrange plots
legend.1 <- get_legend(loc_box_plot)
legend.2 <- get_legend(gwc_plot)

both <- plot_grid(minN_plot + theme(legend.position="none"), 
                  gwc_plot + theme(legend.position = "none"), 
                  labels=c("B", "C"), nrow=1, align="v")
both
all <- plot_grid(loc_box_plot + theme(legend.position="none"), 
                 both, nrow=2, align="v",
                 labels=c("A",""))

all

legends <- plot_grid(legend.1, legend.2, ncol=1, align="v")

all_final <- plot_grid(all, legends, nrow=1, rel_widths = c(6, 1))
all_final

ggsave("ndfa_plots.png", dpi=300, bg="white", 
       width=12, height=8, units="in")

