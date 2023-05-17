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
shared <- readRDS("final_shared_rare.rds")

all_samps <- shared %>% 
  mutate(sample = str_split(sample, "R", simplify=TRUE)[,1]) %>%
  # dplyr::select(sample) %>%
  unique() %>%
  mutate(sample = as.double(sample)) %>%
  dplyr::select(sample, name, value)
str(all_samps)


#find missing samps in each dataset
missing_forage <- anti_join(all_samps, forage_samps) %>% #17 obs
  dplyr::select(sample) %>%
  unique() 

missing_samps <- anti_join(forage_samps, all_samps) #10 samples from seq not present


#forage dataframe with 409 samples
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
forage_mat <- final_forage_df %>% #409 obs
  dplyr::select(-location) %>%
  arrange(sample) %>%
  column_to_rownames("sample")

forage_dist <- vegdist(scale(forage_mat), method="euclidean")

#Burera forage distance matrix
burera_forage_mat <- burera_forage %>% #109 obs
  dplyr::select(-location) %>%
  mutate(sample = as.numeric(sample)) %>%
  arrange(sample) %>%
  column_to_rownames("sample") 
burera_forage_dist <- vegdist(scale(burera_forage_mat), method = "euclidean")

#Karama forage distance matrix
karama_forage_mat <- karama_forage %>% #155 obs
  dplyr::select(-location) %>%
  mutate(sample = as.numeric(sample)) %>%
  arrange(sample) %>%
  column_to_rownames("sample") 
karama_forage_dist <- vegdist(scale(karama_forage_mat), method = "euclidean")

#Rubona forage distance matrix
rubona_forage_mat <- rubona_forage %>% #145 obs
  dplyr::select(-location) %>%
  mutate(sample = as.numeric(sample)) %>%
  arrange(sample) %>%
  column_to_rownames("sample") 
rubona_forage_dist <- vegdist(scale(rubona_forage_mat), method = "euclidean")

#### 16s distance matrices ####
#16s dataframe with 413 samps
wide_samps <- all_samps %>%
  # mutate(sample = as.character(sample)) %>%
  group_by(sample, name) %>%
  mutate(value = mean(value)) %>%
  ungroup() %>%
  unique()  #one sample was sequenced twice

#final 16s wide mat  
wide_2_samps <- wide_samps%>% #426 obs
  pivot_wider(names_from=name, values_from=value, values_fill=0)  

final_merged_samps <- semi_join(wide_2_samps, forage_samps, by="sample") %>% #409 - correct number
  arrange(sample)

final_mat <- final_merged_samps %>%
  column_to_rownames("sample")

final_dist <- vegdist(final_mat, method="bray")

#final Burera  16s mat
b_wide_samps <- wide_samps %>% #109
  semi_join(burera_samps)%>%
  pivot_wider(names_from=name, values_from=value, values_fill=0) %>%
  arrange(sample) %>%
  column_to_rownames("sample")

burera_dist <- vegdist(b_wide_samps, method = "bray")

#final Karama 16s mat
k_wide_samps <- wide_samps %>% #155 obs
  semi_join(karama_samps) %>%
  pivot_wider(names_from = name, values_from = value, values_fill = 0)%>%
  arrange(sample) %>%
  column_to_rownames("sample")

karama_dist <- vegdist(k_wide_samps, method="bray")

#final Rubona 16s mat
r_wide_samps <- wide_samps %>% #145 obs
  semi_join(rubona_samps) %>%
  pivot_wider(names_from = name, values_from = value, values_fill = 0)%>%
  arrange(sample) %>%
  column_to_rownames("sample")

rubona_dist <- vegdist(r_wide_samps, method="bray")

#### Mantel test: all samples ####
#forage_dist, final_dist

#is forage quality related to rhizosphere fungal community?
all_16s_mantel <- mantel(final_dist, forage_dist, method="spearman")
all_16s_mantel
# Mantel statistic r: 0.03076 
# Significance: 0.026 

##Karama
k_16s_mantel <- mantel(karama_dist ,karama_forage_dist, method= "spearman")
k_16s_mantel
# Mantel statistic r: -0.005709 
# Significance: 0.566 

#Rubona
r_16s_mantel <- mantel(rubona_dist, rubona_forage_dist, method="spearman")
r_16s_mantel 
# Mantel statistic r: 0.01654 
# Significance: 0.306 

#Burera
b_16s_mantel <- mantel(burera_dist, burera_forage_dist, method="spearman")
b_16s_mantel
# Mantel statistic r: 0.1081 
# Significance: 0.022 