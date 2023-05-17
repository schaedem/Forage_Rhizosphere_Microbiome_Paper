library(tidyverse)
library(vegan)
library(rebus)

#Read in metadata
meta <- read_csv("final_samps.csv") %>%
  dplyr::select(sample, location) %>%
  mutate(sample = str_split(sample, "R", simplify=TRUE)[,1],
         sample = as.numeric(sample))

# #forage distance matrix
# forage.dist <- readRDS("forage.dist.RDS")

#forage dataframe
forage_df <- read_csv("final_forage_data.csv")  

forage_samps <- dplyr::select(forage_df, sample)

#overall mantel (all locations)
shared <- readRDS("final_shared_rare_ITS2.rds")

all_samps <- shared %>% 
  mutate(sample = str_split(sample, "R", simplify=TRUE)[,1]) %>%
 # dplyr::select(sample) %>%
  unique() %>%
  mutate(sample = as.double(sample)) %>%
  dplyr::select(sample, name, value)
str(all_samps)


#find missing samps in each dataset
missing_forage <- anti_join(all_samps, forage_samps) %>%
  dplyr::select(sample) %>%
  unique() 

missing_samps <- anti_join(forage_samps, all_samps) #6 samples from seq not present


#forage dataframe with 413 samples
final_forage_df <- forage_df %>%
  anti_join(missing_samps) %>%
  arrange(sample) %>%
  merge(meta, by="sample") %>%
  unique()

#merge with metadata to sort by location
rubona_forage <- final_forage_df %>%
  filter(location == "Rubona")
rubona_samps <- rubona_forage%>%dplyr::select(sample)

karama_forage <- final_forage_df %>%
  filter(location == "Karama")
karama_samps <- karama_forage%>%dplyr::select(sample)

burera_forage <- final_forage_df %>%
  filter(location == "Burera")
burera_samps <- burera_forage%>%dplyr::select(sample)

#### Forage distance matrices ####
#forage distance matrix
forage_mat <- forage_df %>%
  column_to_rownames("sample")

forage_dist <- vegdist(scale(forage_mat), method="euclidean")

#Burera forage distance matrix
burera_forage_mat <- burera_forage %>% #110 obs
  dplyr::select(-location) %>%
  mutate(sample = as.numeric(sample)) %>%
  arrange(sample) %>%
  column_to_rownames("sample") 
burera_forage_dist <- vegdist(scale(burera_forage_mat), method = "euclidean")

#Karama forage distance matrix
karama_forage_mat <- karama_forage %>% #156 obs
  dplyr::select(-location) %>%
  mutate(sample = as.numeric(sample)) %>%
  arrange(sample) %>%
  column_to_rownames("sample") 
karama_forage_dist <- vegdist(scale(karama_forage_mat), method = "euclidean")

#Rubona forage distance matrix
rubona_forage_mat <- rubona_forage %>% #147 obs
  dplyr::select(-location) %>%
  mutate(sample = as.numeric(sample)) %>%
  arrange(sample) %>%
  column_to_rownames("sample") 
rubona_forage_dist <- vegdist(scale(rubona_forage_mat), method = "euclidean")

#### ITS2 distance matrices ####
#ITS2 dataframe with 413 samps
wide_samps <- all_samps %>%
 # mutate(sample = as.character(sample)) %>%
  group_by(sample, name) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  unique()  #one sample was sequenced twice

#final ITS2 wide mat  
wide_2_samps <- wide_samps%>%
  pivot_wider(names_from=name, values_from=value, values_fill=0)  #430 - correct number

final_merged_samps <- semi_join(wide_2_samps, forage_samps, by="sample") %>% #413 - correct number
  arrange(sample)

final_mat <- final_merged_samps %>%
  column_to_rownames("sample")

final_dist <- vegdist(final_mat, method="bray")

#final Burera  ITS2 mat
b_wide_samps <- wide_samps %>%
  semi_join(burera_samps)%>%
  pivot_wider(names_from=name, values_from=value, values_fill=0) %>%
  arrange(sample) %>%
  column_to_rownames("sample")

burera_dist <- vegdist(b_wide_samps, method = "bray")

#final Karama ITS2 mat
k_wide_samps <- wide_samps %>% #156 obs
  semi_join(karama_samps) %>%
  pivot_wider(names_from = name, values_from = value, values_fill = 0)%>%
  arrange(sample) %>%
  column_to_rownames("sample")

karama_dist <- vegdist(k_wide_samps, method="bray")

#final Rubona ITS2 mat
r_wide_samps <- wide_samps %>% #147 obs
  semi_join(rubona_samps) %>%
  pivot_wider(names_from = name, values_from = value, values_fill = 0)%>%
  arrange(sample) %>%
  column_to_rownames("sample")

rubona_dist <- vegdist(r_wide_samps, method="bray")

#### Mantel test: all samples ####
#forage_dist, final_dist

#is forage quality related to rhizosphere fungal community?
all_its2_mantel <- mantel(final_dist, forage_dist, method="spearman")
all_its2_mantel
# Mantel statistic r: 0.08419 
# Significance: 0.001 

##Karama
k_its2_mantel <- mantel(karama_dist ,karama_forage_dist, method= "spearman")

# Mantel statistic r: 0.188 
# Significance: 0.001 

#Rubona
r_its2_mantel <- mantel(rubona_dist, rubona_forage_dist, method="spearman")
# Mantel statistic r: 0.1382 
# Significance: 0.001 

#Burera
b_its2_mantel <- mantel(burera_dist, burera_forage_dist, method="spearman")
# Mantel statistic r: 0.2196 
# Significance: 0.001 