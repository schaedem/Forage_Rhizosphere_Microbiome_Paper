library(tidyverse)
library(RColorBrewer)
source("calc_lrr.R")
display.brewer.pal(10, "Paired")
brewer.pal(10, "Paired")

leg_pal <- c("#A6CEE3", "#B2DF8A", "#FB9A99")
nonleg_pal <- c("#1F78B4", "#33A02C", "#E31A1C")

#### N LRR ####
legume_N_lrr <- ggplot(leg_merge, aes(x=treatment, y=LRR, color=treatment)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="darkred", size=0.5) +
  facet_grid(location ~ timepoint) +
  xlab("") +
  ylab("LRR (%N)") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = leg_pal) +
  geom_jitter()
legume_N_lrr


label_df <- cbind(species = "Maize", label = "**", location = "Karama") %>% as_tibble()

nonleg_N_lrr <- ggplot(non_leg_merge, aes(x=species, y=LRR, color=species)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  scale_color_manual(values = nonleg_pal) +
  geom_jitter() +
  facet_grid(~location) +
  xlab("") +
  ylab("LRR (%N)") +
  theme(legend.title = element_blank()) +
  geom_text(data=label_df, aes(x=species, y=1, label=label), 
            nudge_x = 0.25, size=8, inherit.aes=F, color="black", fontface="italic")
  
nonleg_N_lrr
ggsave("LRR_N.png", dpi=300, width=8, height=2, units="in")

##### Protein LRR #####

legume_protein_lrr <- ggplot(leg_merge_nir, aes(x=treatment, y=protein_LRR, color=treatment)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="darkred", size=0.5) +
  facet_grid(location ~ timepoint) +
  xlab("") +
  ylab("LRR (%Protein)") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = leg_pal) +
  geom_jitter() 
legume_protein_lrr

protein_label_df <- cbind(species = "Maize", label = ".", location = "Rubona") %>% as_tibble()

nonleg_protein_lrr <- ggplot(non_leg_merge_nir, aes(x=species, y=protein_LRR, color=species)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  scale_color_manual(values = nonleg_pal) +
  geom_jitter() +
  facet_grid(~location) +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  geom_text(data=protein_label_df, aes(x=species, y=0.8, label=label), 
            size=8, color="black", nudge_x = 0.5)

nonleg_protein_lrr

##### Lignin LRR #####
legume_lignin_lrr <- ggplot(leg_merge_nir, aes(x=treatment, y=lignin_LRR, color=treatment)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  facet_grid(location ~ timepoint) +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = leg_pal) +
  geom_jitter()
legume_lignin_lrr


lignin_label_df <- cbind(species = c("Napier", "Brachiaria"), 
                         label = c(".", "*"), 
                         location = c("Rubona", "Burera")) %>% as_tibble()
nonleg_lignin_lrr <- ggplot(non_leg_merge_nir, aes(x=species, y=lignin_LRR, color=species)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  scale_color_manual(values = nonleg_pal) +
  geom_jitter() +
  facet_grid(~location) +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  geom_text(data=lignin_label_df, aes(x=species, y=1.3, label=label), 
            size=8, color="black", nudge_x = 0.25)


nonleg_lignin_lrr

##### ADF LRR #####
legume_adf_lrr <- ggplot(leg_merge_nir, aes(x=treatment, y=adf_LRR, color=treatment)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  facet_grid(~location) +
  xlab("") +
  #ylab("LRR (%ADF)") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = leg_pal) +
  geom_jitter()
legume_adf_lrr

nonleg_adf_lrr <- ggplot(non_leg_merge_nir, aes(x=species, y=adf_LRR, color=species)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  scale_color_manual(values = nonleg_pal) +
  geom_jitter() +
  facet_grid(~location) +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) 

nonleg_adf_lrr

##### NDF LRR #####
legume_ndf_lrr <- ggplot(leg_merge_nir, aes(x=treatment, y=ndf_LRR, color=treatment)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  facet_grid(location ~ timepoint) +
  xlab("") +
  ylab("LRR (%NDF)") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = leg_pal) +
  geom_jitter()
legume_ndf_lrr

#NDF plot
ndf_label_df <- cbind(species = "Maize", label = ".", location = "Karama") %>% as_tibble()

ndf_rect <- cbind(species = "Maize", location = "Karama", 
                  ymin = -Inf, ymax = Inf, xmin="Brachiaria", xmax="Napier") %>% as_tibble()

nonleg_ndf_lrr <- ggplot(non_leg_merge_nir, aes(x=species, y=ndf_LRR, color=species)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +
  theme_test() +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  scale_color_manual(values = nonleg_pal) +
  geom_jitter() +
  facet_grid(~location) +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  geom_text(data=ndf_label_df, aes(x=species, y=0.18, label=label), 
            size=8, color="black", nudge_x = 0.25) 

nonleg_ndf_lrr

####Collate nonleg plots ####

legend_nonleg <- get_legend(nonleg_ndf_lrr)

nonleg_lrr <- plot_grid(nonleg_N_lrr + theme(legend.position = "none",
                                             axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank()),
                        nonleg_protein_lrr+ theme(legend.position = "none",
                                                  axis.text.y = element_blank(),
                                                  axis.ticks.y = element_blank()),
                        nonleg_lignin_lrr+ theme(legend.position = "none",
                                                 axis.text.y = element_blank(),
                                                 axis.ticks.y = element_blank()),
                        nonleg_adf_lrr+ theme(legend.position = "none",
                                              axis.text.y = element_blank(),
                                              axis.ticks.y = element_blank()),
                        nonleg_ndf_lrr+ theme(legend.position = "none",
                                              axis.text.y = element_blank(),
                                              axis.ticks.y = element_blank())
                                      + ylab("LRR"),
                        labels="AUTO",
                        ncol=1, align = "h")
nonleg_lrr

nonleg_lrr_legend <- plot_grid(nonleg_lrr, legend_nonleg, ncol=2,
                               rel_widths = c(4,0.5))
nonleg_lrr_legend

ggsave("nonleg_lrr_plots.png", dpi=300, bg="white",
       width=10, height=11, units="in")
