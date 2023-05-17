library(tidyverse)
library(indicspecies)
library(vegan)
library(rebus) #constants for special characters

test <- data %>%
  group_by(sample) %>%
  summarize(N = sum(value)) #checking that data is indeed rarefied

data <- read_rds("full_shared_rare_16S.rds")
print(data, n=50)
levels(data$species)

data <- data %>%
  mutate(treatment=as.character(treatment),
         intercrop = ifelse(grepl(PLUS, treatment), "yes", "no"))

##separate bact OTUs by location
burera <- data %>% filter(location=="Burera")

rubona <- data %>% filter(location=="Rubona")

karama <- data %>% filter(location=="Karama")
which(karama$sample=="222R")

#function that separates otu dataframe into separate dfs by species
#species dataframes are returned as a list

prep_otu <- function(loc_data){
  brachiaria <- loc_data %>% filter(species=="Brachiaria")
  desmodium <- loc_data %>% filter(species=="Desmodium")
  maize <- loc_data %>% filter(species=="Maize")
  napier <- loc_data %>% filter(species=="Napier")
  list <- list("b" = brachiaria, "d"= desmodium, "m"= maize,"n"= napier)
  
  return(list)
  
}

#run function on each location df
k_prep <- prep_otu(karama)
r_prep <- prep_otu(rubona)
b_prep <- prep_otu(burera)

##### Indicator Species Analysis ####

#Function that:
#1. extracts group information (intercrop = yes/no)
#2. pivots species dataframe into wide format
#3. runs multipatt to find indicator species
#4. extracts results, filters by p<0.01
#5. returns a dataframe with associated location and species

run_indic <- function(loc_sp_df, location, species) {
  df <- loc_sp_df %>%
    dplyr::select(-c(location, treatment, species, timepoint, block))
  
  group <- df %>%
    dplyr::select(sample, intercrop) %>%
    unique()
  
  df_final <- df %>%
    tidyr::pivot_wider(names_from = name, values_from = value) %>%
    column_to_rownames("sample") %>%
    dplyr::select(-intercrop)
  
  df_filter <- df_final[rowSums(df_final) > 0,]
  
  set.seed(123455)
  indic <- multipatt(df_filter, group$intercrop, 
                     control=how(nperm=1000))
  
  #extract stats table
  inval_sign <- indic$sign %>%
    rownames_to_column(var="otu") %>%
    filter(p.value < 0.01) %>%
    mutate(location = paste(location),
           species = paste(species))
  
  return(inval_sign)
  
}

#Karama
k_brach <- run_indic(k_prep$b, "Karama", "Brachiaria")
k_nap <- run_indic(k_prep$n, "Karama", "Napier") 
k_maize <- run_indic(k_prep$m, "Karama", "Maize")
k_des <- run_indic(k_prep$d, "Karama", "Desmodium") #no bacteria enriched in intercropped desmodium

k_results <- rbind(k_brach, k_nap, k_maize, k_des)
#write_csv(k_results, "karama_indicator_species.csv")
write_csv(k_results, "karama_indicator_species_2.csv")

#Rubona
r_brach <- run_indic(r_prep$b, "Rubona", "Brachiaria")
r_nap <- run_indic(r_prep$n, "Rubona", "Napier")
r_maize <- run_indic(r_prep$m, "Rubona", "Maize")
r_des <- run_indic(r_prep$d, "Rubona", "Desmodium") #1 otu is enriched in intercropped desmodium

r_results <- rbind(r_brach, r_nap, r_maize, r_des)
#write_csv(r_results, "rubona_indicator_species.csv")
write_csv(r_results, "rubona_indicator_species_2.csv")

#Burera
b_brach <- run_indic(b_prep$b, "Burera", "Brachiaria")
b_nap <- run_indic(b_prep$n, "Burera", "Napier")
b_maize <- run_indic(b_prep$m, "Burera", "Maize")
b_des <- run_indic(b_prep$d, "Burera", "Desmodium") #only 1 otu associated with intercropping

b_results <- rbind(b_brach, b_nap, b_maize, b_des)
#write_csv(b_results, "burera_indicator_species.csv")
write_csv(b_results, "burera_indicator_spcecies_2.csv")

