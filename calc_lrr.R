library(tidyverse)
#source("calc_ndfa.R")

#### LRR : %N ####
#read in irms data
irms_full <- read_csv("irms_final_full.csv")

#Getting legume dfs
leg_mono <- irms_full %>% 
  filter(intercrop != "yes") %>%
  filter(grepl("Desmodium", species)) %>%
  dplyr::select(timepoint, species, location, rep, N.corrected) %>%
  rename(mono_N = N.corrected)
         #legume = species)

leg_int <- irms_full %>%
  filter(species == "Desmodium")%>%
  filter(intercrop == "yes") %>%
  mutate(treatment = as.character(treatment)) %>%
  dplyr::select(timepoint, species, location, rep, N.corrected, treatment) %>%
  rename(int_N = N.corrected)

leg_merge <- inner_join(leg_mono, leg_int, by=c("rep", "timepoint", "location", "species")) %>%
  mutate(LRR = log(int_N/mono_N)) #response ratio for magnitude diff of intercropping

#Getting nonlegume dfs
non_leg_mono <- irms_full %>%
  filter(intercrop != "yes") %>%
  filter(!grepl("Desmodium", species)) %>%
  dplyr::select(timepoint, species, location, rep, N.corrected, treatment) %>%
  rename(mono_N = N.corrected)

non_leg_mix <- irms_full %>%
  filter(intercrop == "yes") %>%
  dplyr::select(timepoint, species, location, rep, N.corrected, treatment) %>%
  rename(int_N = N.corrected)

non_leg_merge <- inner_join(non_leg_mono, non_leg_mix,
                            by=c("rep", "timepoint", "location", "species")) %>%
  dplyr::select(-treatment.x) %>%
  rename(treatment = treatment.y) %>%
  mutate(LRR = log(int_N/mono_N)) #log response ratio 



#getting total dataframe for N LRR

nler_total <- rbind(leg_merge, non_leg_merge)
summary(nler_total$LRR)


#### LRR: protein ####
nir_full <- read_csv("nir_final.csv") %>%
  dplyr::select(sample, protein, adf, ndf, lignin, timepoint) %>%
  inner_join(irms_full, by=c("sample", "timepoint")) %>%
  dplyr::select(sample, protein, adf, ndf, lignin, intercrop, species, treatment, location, timepoint, rep)

#Getting legume dfs
leg_mono_nir <- nir_full %>% 
  filter(intercrop != "yes") %>%
  filter(grepl("Desmodium", species)) %>%
  rename(mono_protein = protein,
         mono_lignin = lignin,
         mono_ndf = ndf,
         mono_adf = adf)

leg_int_nir <- nir_full %>%
  filter(species == "Desmodium")%>%
  filter(intercrop == "yes") %>%
  rename(int_protein = protein,
         int_lignin = lignin,
         int_ndf = ndf,
         int_adf = adf)

leg_merge_nir <- inner_join(leg_mono_nir, leg_int_nir,
                            by=c("rep", "timepoint", "location", "species")) %>%
  mutate(protein_LRR = log(int_protein/mono_protein),
         lignin_LRR = log(int_lignin/mono_lignin),
         adf_LRR = log(int_adf/mono_adf),
         ndf_LRR = log(int_ndf/mono_ndf)) %>% #response ratio for magnitude diff of intercropping
  dplyr::select(-c(intercrop.x, intercrop.y, 
                   sample.x, sample.y, treatment.x)) %>%
  dplyr::select(species, protein_LRR, lignin_LRR, adf_LRR, ndf_LRR,
                timepoint, rep, location, treatment.y) %>%
  rename(treatment = treatment.y) %>%
  mutate(lignin_LRR = ifelse(is.nan(lignin_LRR), 0, lignin_LRR))

#Getting nonlegume dfs
non_leg_mono_nir <- nir_full %>%
  filter(intercrop != "yes") %>%
  filter(!grepl("Desmodium", species)) %>%
  dplyr::select(timepoint, species, location, rep, 
                protein, lignin, adf, ndf, treatment) %>%
  rename(mono_protein = protein,
         mono_lignin = lignin,
         mono_adf = adf, 
         mono_ndf = ndf)

non_leg_mix_nir <- nir_full %>%
  filter(intercrop == "yes") %>%
  dplyr::select(timepoint, species, location, rep,treatment,
                protein, lignin, adf, ndf) %>%
  rename(int_protein = protein,
         int_lignin = lignin,
         int_adf = adf, 
         int_ndf = ndf)

non_leg_merge_nir <- inner_join(non_leg_mono_nir, non_leg_mix_nir,
                            by=c("rep", "timepoint", "location", "species")) %>%
  dplyr::select(-treatment.x) %>%
  rename(treatment = treatment.y) %>%
  mutate(protein_LRR = log(int_protein/mono_protein),
         lignin_LRR = log(int_lignin/mono_lignin),
         adf_LRR = log(int_adf/mono_adf),
         ndf_LRR = log(int_ndf/mono_ndf)) %>% #log response ratio 
  dplyr::select(species,  protein_LRR, lignin_LRR, adf_LRR, ndf_LRR,
                rep, timepoint, location, species, treatment)


#getting total dataframe for N LRR

lrr_total_nir <- rbind(leg_merge_nir, non_leg_merge_nir) %>%
  mutate(lignin_LRR = ifelse(is.nan(lignin_LRR), 0, lignin_LRR))


#### Land equivalency ratio ####
#calculating LER (land equivalency ratio)
#LER = (Nlm / Nl) + (Ngm/Ng)
#Nlm = N leg mixed
#Nl = N leg mono
#Ngm = N grass mixed
#Ng = N grass mono

# nler_full <- inner_join(leg_merge, non_leg_merge, by=c("rep", "timepoint", "location", "treatment")) %>%
#   rename(legume = species.x,
#          non_leg = species.y) %>%
#   mutate(nler = (leg_int_N / leg_mono_N) + (N.nonleg.mix/N.nonleg.mono)) %>%
#   mutate(non_leg = as.character(non_leg)) %>%
#   filter(non_leg != "Desmodium")
# 
# head(nler_full)
# 
# #### Graphing LER ####
# pal <- c("#E31A1C", "#1F78B4","#33A02C")
# levels(nler_full$non_leg)
# 
# ggplot(nler_full, aes(x=non_leg, y=nler, color=non_leg)) +
#   geom_boxplot() +
#   #geom_jitter(stat="identity", alpha=0.4) +
#   facet_grid(location~timepoint) +
#   theme_classic() +
#   geom_hline(yintercept = 1, color="darkred", alpha=0.5, size=1, linetype=2) + 
# #  scale_fill_brewer(palette = "Paired") +
#   scale_fill_manual(values = pal) +
#   ylab("N Land Use Efficiency") +
#   theme(legend.position="none") +
#   xlab("")
# 
# #graphing N benefit from intercropping
# t <- rep(1:4, times=9)
# sp <- rep(c("Brachiaria", "Maize", "Napier"), times=12)
# l <-rep(c("Burera", "Karama", "Rubona"), each=12)
# 
# #creating sig dataframe to annotate facets
# sig <- cbind(t, sp, l) %>%
#   as_tibble() %>%
#   rename(timepoint = t, species=sp, location=l) %>%
#   mutate(drop = ifelse(timepoint==4 & location=="Burera", "y",
#                        ifelse(timepoint==4 & location=="Karama", "y", "n"))) %>%
#     filter(drop=="n") %>%
#   select(-drop) %>%
#   mutate(label = ifelse(species=="Maize" & location=="Karama" & timepoint==3,
#                         "*", 
#                         ifelse(species=="Brachiaria" & location=="Burera" & timepoint==2, "*", ""))) %>%
#   mutate(ypos = ifelse(species=="Maize" & location=="Karama" & timepoint==3,
#                         250, 
#                         ifelse(species=="Brachiaria" & location=="Burera" & timepoint==2, 100, 0)))
# #results from n_benefit_mods.R, lme comparing trt means to zero
# 
