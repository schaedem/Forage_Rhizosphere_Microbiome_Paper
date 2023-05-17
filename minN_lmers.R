library(tidyverse)
library(nlme)
library(lsmeans)
library(broom)
library(broom.mixed)
library(broomExtra)

minN <- read_csv("bnf_minN_final.csv") %>%
  mutate(treatment = str_trim(treatment, side="both"))

rubona <- minN %>%
  filter(location == "Rubona")

karama <- minN %>%
  filter(location == "Karama")

burera <- minN %>%
  filter(location =="Burera") %>%
  filter(!is.na(minN))

#### 1. Overall by Location ####

#1a. Rubona
lmer1 <- lme(minN ~ treatment -1, random = c(~1|timepoint, ~1|rep), data=rubona)
summary(lmer1)

car::Anova(lmer1) #treatment highly sig (p<0.001)

#treatment is significant (p=0.046)
tuk1 <- emmeans(lmer1, pairwise~treatment, adjust="tukey") 
lmer1_letters <- multcomp::cld(object=tuk1$emmeans, Letters=letters) %>% #several differences
  mutate(location = "Rubona")

#1b. Karama
lmer2 <- lme(minN ~ treatment -1, random = c(~1|timepoint, ~1|rep), data=karama)
summary(lmer2) 
car::Anova(lmer2) #treatment highly sig
tuk2 <- emmeans(lmer2, pairwise~treatment, adjust="tukey") 
lmer2_letters <- multcomp::cld(object=tuk2$emmeans, Letters=letters) %>% #no differences
  mutate(location = "Karama")

#1c. Burera
lmer3 <- lme(minN ~ treatment -1, random = c(~1|timepoint, ~1|rep), data=burera)
summary(lmer3)
car::Anova(lmer3) #treatment highly sig
tuk3 <-emmeans(lmer3, pairwise~treatment, adjust="tukey") #no sig diff btw treatments
lmer3_letters <- multcomp::cld(object=tuk3$emmeans, Letters = letters) %>% #no differences
  mutate(location = "Burera")


#export letters for graphing
loc_letters <- rbind(lmer1_letters, lmer2_letters, lmer3_letters)
write_csv(loc_letters, "overall_loc_minN.csv")

#### 2. Within location x timepoint ####

within_loctime <- function(data){
  data_within <- data %>%
    group_by(timepoint) %>%
    do(mod = lme(minN ~ treatment -1, random = ~1|rep, data=.)) 
  data_within$mod[[1]]
  
  #getting Anova p-vals
  sum1 <- tidy(car::Anova(data_within$mod[[1]])) %>% mutate(timepoint =1)
  sum2 <- tidy(car::Anova(data_within$mod[[2]])) %>% mutate(timepoint =2)
  sum3 <- tidy(car::Anova(data_within$mod[[3]])) %>% mutate(timepoint =3)
  #sum4 <- tidy(car::Anova(data_within$mod[[4]]))%>% mutate(timepoint =4)
  
  anova <- rbind(sum1, sum2, sum3)#, sum4)
  
  
  #getting letters for graphing
  l1 <- multcomp::cld(emmeans(data_within$mod[[1]], pairwise~treatment, adjust="tukey")$emmeans, 
                      Letters = letters) %>% mutate(timepoint = 1)

  l2 <- multcomp::cld(emmeans(data_within$mod[[2]], pairwise~treatment, adjust="tukey")$emmeans, 
                      Letters = letters) %>% mutate(timepoint = 2)
  
  l3 <- multcomp::cld(emmeans(data_within$mod[[3]], pairwise~treatment, adjust="tukey")$emmeans, 
                      Letters = letters) %>% mutate(timepoint = 3)
  
  # l4 <- multcomp::cld(emmeans(data_within$mod[[4]], pairwise~treatment, adjust="tukey")$emmeans,
  #                     Letters = letters) %>% mutate(timepoint = 4)

  letters_all <- rbind(l1, l2, l3)#, l4)
  
  return(list("anova" = anova, "letters" = letters_all))
  
}

rubona_minN <-within_loctime(rubona)
rubona_let <-rubona_minN$letters %>% #some differences
  mutate(location = "Rubona") %>%
  mutate(.group = str_trim(.group, side="both"))
rubona_minN$anova #all sig diff from zero

karama_minN <- within_loctime(karama)
karama_let <- karama_minN$letters %>% #no differences
  mutate(location = "Karama")%>%
  mutate(.group = str_trim(.group, side="both"))
karama_minN$anova #all sig diff from zero

burera_minN <- within_loctime(burera)
burera_let <- burera_minN$letters %>% #some differences in T2
  mutate(location = "Burera")
burera_minN$anova

graph_loctime <- rbind(rubona_let, karama_let, burera_let) %>%
  mutate(timepoint = paste("Harvest", timepoint, sep=" ")) %>%
  mutate(.group = ifelse(location == "Karama", "", .group),
         ypos = emmean + upper.CL,
         ypos = ifelse(location=="Burera" , ypos-7, ypos))

write_csv(graph_loctime,"minN_emmeans_treatment.csv")
