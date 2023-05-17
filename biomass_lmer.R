library(tidyverse)
library(emmeans)
library(multcomp)
library(lme4)
source("biomass_graphing.R")

#Separate data frames by location
rubona <- final %>% #160 obs
  filter(location == "Rubona")

karama <- final %>% #157 obs
  filter(location == "Karama")

burera <- final %>% #114 obs
  filter(location == "Burera")

#### Overall by Location - Biomass ####

#### Rubona ####
#N
lmer1 <- lmer(N.corrected ~ species-1 + (1|timepoint) + (1|rep), data=rubona)
summary(lmer1)  
car::Anova(lmer1) #sig (p<0.001)
rubona_N <- cld(emmeans(lmer1, list(pairwise~species), adjust="tukey"),Letters=letters )

#d15N
lmer20 <- lmer(d15N.corrected ~ species -1 + (1|timepoint) + (1|rep), data=rubona)
summary(lmer20)
car::Anova(lmer20) #p<0.001
rubona_d15N <- cld(emmeans(lmer20, list(pairwise~species), adjust="tukey"), Letters=letters)

#protein
lmer4 <- lmer(protein ~ species -1 + (1|timepoint) + (1|rep), data=rubona)
summary(lmer4) 
car::Anova(lmer4) #sig (p<0.001)
rubona_protein <- cld(emmeans(lmer4, list(pairwise~species), 
                              adjust="tukey"), Letters=letters) 

#lignin
lmer7 <- lmer(lignin ~ species -1 + (1|timepoint) + (1|rep), data=rubona)
summary(lmer7) 
car::Anova(lmer7) #sig (p<0.001)
rubona_lignin <- cld(emmeans(lmer7, list(pairwise~species),
                             adjust="tukey"), Letters=letters)

#adf
lmer8 <- lmer(adf ~ species -1 + (1|timepoint) + (1|rep), data=rubona)
summary(lmer8) 
car::Anova(lmer8) #p<0.001
rubona_adf <- cld(emmeans(lmer8, list(pairwise~species),
                          adjust="tukey"), Letters=letters)

#ndf
lmer9 <- lmer(ndf ~ species -1 + (1|timepoint) + (1|rep), data=rubona)
summary(lmer9) 
car::Anova(lmer9) #p<0.001
rubona_ndf <- cld(emmeans(lmer9, list(pairwise~species),
                          adjust="tukey"), Letters=letters)

rubona_cld <- rbind(rubona_N %>% mutate(quality= "N"), 
                    rubona_d15N %>% mutate(quality="d15N"),
                    rubona_protein %>% mutate(quality = "protein"),
                    rubona_adf %>% mutate(quality = "adf"),
                    rubona_ndf %>% mutate(quality = "ndf"),
                    rubona_lignin %>% mutate(quality = "lignin")) %>%
  dplyr::select(species, .group, quality) %>%
  mutate(location = "Rubona",
         .group = str_trim(.group, side = "both"))

##### Karama ####
#N
lmer10 <- lmer(N.corrected ~ species-1 + (1|timepoint) + (1|rep), data=karama)
summary(lmer10)  
car::Anova(lmer10) #sig (p<0.001)
karama_N <- cld(emmeans(lmer10, list(pairwise~species), adjust="tukey"),
                Letters=letters ) #Des is sig higher than others

#d15N
lmer11 <- lmer(d15N.corrected ~ species -1 + (1|timepoint) + (1|rep), data=karama)
summary(lmer11)
car::Anova(lmer11) #p<0.001
karama_d15N <- cld(emmeans(lmer11, list(pairwise~species), adjust="tukey"), 
                   Letters=letters) #desmodium sig lower than others

#protein
lmer12 <- lmer(protein ~ species -1 + (1|timepoint) + (1|rep), data=karama)
summary(lmer12) 
car::Anova(lmer12) #sig (p<0.001)
karama_protein <- cld(emmeans(lmer12, list(pairwise~species), 
                              adjust="tukey"), Letters=letters) #Desmodium highest, napier lowest

#lignin
lmer13 <- lmer(lignin ~ species -1 + (1|timepoint) + (1|rep), data=karama)
summary(lmer13) 
car::Anova(lmer13) #sig (p<0.001)
karama_lignin <- cld(emmeans(lmer13, list(pairwise~species),
                             adjust="tukey"), Letters=letters) #Des highest

#adf
lmer14 <- lmer(adf ~ species -1 + (1|timepoint) + (1|rep), data=karama)
summary(lmer14) 
car::Anova(lmer14) #p<0.001
karama_adf <- cld(emmeans(lmer14, list(pairwise~species),
                          adjust="tukey"), Letters=letters) #maize highest

#ndf
lmer15 <- lmer(ndf ~ species -1 + (1|timepoint) + (1|rep), data=karama)
summary(lmer15) 
car::Anova(lmer15) #p<0.001
karama_ndf <- cld(emmeans(lmer15, list(pairwise~species),
                          adjust="tukey"), Letters=letters)

karama_cld <- rbind(karama_N %>% mutate(quality= "N"), 
                    karama_d15N %>% mutate(quality="d15N"),
                    karama_protein %>% mutate(quality = "protein"),
                    karama_adf %>% mutate(quality = "adf"),
                    karama_ndf %>% mutate(quality = "ndf"),
                    karama_lignin %>% mutate(quality = "lignin")) %>%
  dplyr::select(species, .group, quality) %>%
  mutate(location = "Karama",
         .group = str_trim(.group, side = "both"))

##### Burera ####
#N
lmer16 <- lmer(N.corrected ~ species-1 + (1|timepoint) + (1|rep), data=burera)
summary(lmer16)  
car::Anova(lmer16) #sig (p<0.001)
burera_N <- cld(emmeans(lmer16, list(pairwise~species), adjust="tukey"),
                Letters=letters ) #Des is sig higher than others

#d15N
lmer17 <- lmer(d15N.corrected ~ species -1 + (1|timepoint) + (1|rep), data=burera)
summary(lmer17)
car::Anova(lmer17) #p<0.001
burera_d15N <- cld(emmeans(lmer17, list(pairwise~species), adjust="tukey"), 
                   Letters=letters) #desmodium sig lower than others

#protein
lmer18 <- lmer(protein ~ species -1 + (1|timepoint) + (1|rep), data=burera)
summary(lmer18) 
car::Anova(lmer18) #sig (p<0.001)
burera_protein <- cld(emmeans(lmer18, list(pairwise~species), 
                              adjust="tukey"), Letters=letters) #Desmodium highest, maize lowest

#lignin
lmer19 <- lmer(lignin ~ species -1 + (1|timepoint) + (1|rep), data=burera)
summary(lmer19) 
car::Anova(lmer19) #sig (p<0.001)
burera_lignin <- cld(emmeans(lmer19, list(pairwise~species),
                             adjust="tukey"), Letters=letters) #Des & maize highest

#adf
lmer20 <- lmer(adf ~ species -1 + (1|timepoint) + (1|rep), data=burera)
summary(lmer20) 
car::Anova(lmer20) #p<0.001
burera_adf <- cld(emmeans(lmer20, list(pairwise~species),
                          adjust="tukey"), Letters=letters) #maize/napier highest

#ndf
lmer21 <- lmer(ndf ~ species -1 + (1|timepoint) + (1|rep), data=burera)
summary(lmer21) 
car::Anova(lmer21) #p<0.001
burera_ndf <- cld(emmeans(lmer21, list(pairwise~species),
                          adjust="tukey"), Letters=letters) #maize/napier highest

burera_cld <- rbind(burera_N %>% mutate(quality= "N"), 
                    burera_d15N %>% mutate(quality="d15N"),
                    burera_protein %>% mutate(quality = "protein"),
                    burera_adf %>% mutate(quality = "adf"),
                    burera_ndf %>% mutate(quality = "ndf"),
                    burera_lignin %>% mutate(quality = "lignin")) %>%
  dplyr::select(species, .group, quality) %>%
  mutate(location = "Burera",
         .group = str_trim(.group, side = "both"))

all_cld <- rbind(burera_cld, karama_cld, rubona_cld)
write_csv(all_cld, "biomass_cld.csv")
