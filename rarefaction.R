library(tidyverse)
library(stringr)
library(vegan)
library(schtools)

count <- read_csv("count.summary.csv") %>%
  rename(sample = group)

meta <- read_csv("final_samps.csv") %>%
  mutate(treatment = as.character(treatment),
    treatment = str_replace_all(treatment, "�", ""),
    species = as.factor(species),
    treatment = ifelse(grepl("Brachiaria v. Mulato II", treatment), "Brachiaria cv. Mulato II",
                       ifelse(grepl("Maize", treatment), "Maize monocrop", treatment)),
    treatment = as.factor(treatment)) 

head(meta)
levels(meta$treatment)

#determine sampling depth for rarefaction
head(count)
summary(count)

count_full <- count %>% merge(meta, by="sample")
head(count_full)

# ggplot(data=count, aes(x=seqs)) +
#   geom_histogram(binwidth = 1000) +
#   coord_cartesian(xlim=c(0,10000))
# #break at ~6,000 or at ~2,500
# 
# ggplot(data=count, aes(x=1, y=seqs)) +
#   geom_jitter() +
#   scale_y_log10()
# #gridline at 1,000
# 
# ggplot(data=count, aes(x=1, y=seqs)) +
#   geom_violin() 

count %>%
  arrange(seqs) %>%
#  print(n=35)
  ggplot(aes(x=1:nrow(.), y=seqs)) + 
  geom_line() +
  coord_cartesian(xlim=c(0,50), ylim=c(0,6000)) +
  geom_hline(yintercept=2097, size=1, color="red") +
  geom_hline(yintercept=1561, size=1, color="green") +
#  geom_hline(yintercept=709, size=1, color="pink") +
 # geom_hline(yintercept=5840, size=1, color="brown") +
  geom_hline(yintercept=2781, size=1, color="yellow")

#candidates: 1561 (obs 8),2097 (obs 9), 2781 (obs 14)
#stay with 2097


#chopping at 5840 would remove many biological replicates

#### Tidy shared file ####

set.seed(19921001)

# shared <- read_tsv("final.opti_mcc.shared") %>%
#   select(Group, starts_with("Otu")) %>%
#   mutate() %>%
#   rename(sample = Group) %>%
#   pivot_longer(-sample) %>%
#   group_by(sample) %>%
#   mutate(total=sum(value)) %>%
#   filter(total > 2097) %>%
#   group_by(name) %>%
#   mutate(total=sum(value)) %>%
#   filter(total!=0) %>%
#   ungroup()%>%
#   select(-total)
# 
# head(shared)
# 
# shared_full <- shared %>%
#   merge(meta, by="sample")
# 
# head(shared_full)
saveRDS(shared_full, "shared_full_16S_final.rds")
#un-rarefied clean shared file with metadata

# rarefying with vegan
#chosen rarefaction sampling depth: 2097
#shared file has 3 columns: sample, name(otu name), and value (n_seqs)
#print(shared, n=1)

#reshape shared
shared_wide <- shared %>%
  pivot_wider(names_from = "name", values_from = "value", values_fill=0) %>%
  as.data.frame()

rownames(shared_wide) <- shared_wide$sample
shared_wide <- shared_wide[,-1]
str(shared_wide) #427 obs, 84307 variables

#rarefy to 2097 sequences
shared_rare <- vegan::rarefy(shared_wide, 2097)
str(shared_rare)
#rarefied species richness per sample

#tidy version
shared_rare_tb <- shared_rare %>%
  as_tibble(rownames="sample")

str(shared_rare_tb)
print(shared_rare_tb, n=1)
saveRDS(shared_rare_tb, "rarefied_richness_final_16S.rds")

#Function rrarefy generates one randomly rarefied community data frame or vector of given sample size. 

#this is the dataframe that will be used for all downstream analyses
# shared_rrarefy <- rrarefy(shared_wide, sample=2097) %>%
#   as_tibble(rownames="sample") %>%
#   pivot_longer(-sample)
# 
# print(shared_rrarefy, n=1)
# 
# shared_rrarefy_full <- shared_rrarefy %>%
#   inner_join(meta, by="sample")
# print(shared_rrarefy_full, n=1)
saveRDS(shared_rrarefy_full, "final_shared_rare_16S.rds")

#drarefy will tell us the probability that each species occurs in a rarefied community
#given the sub-sampling depth

shared_prob <- drarefy(shared_wide, sample=2097) %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample)

head(shared_prob)

#make rarefaction curve
# rarecurve_data <- rarecurve(shared_wide, step=100)
# 
# rarecurve_df <- purrr::map_dfr(rarecurve_data, bind_rows) %>%
#   bind_cols(sample = rownames(shared_wide), .) %>%
#   pivot_longer(-sample) %>%
#   drop_na() %>%
#   mutate(n_seqs = as.numeric(str_replace(name, "N", "")))%>%
#   
# rareplot <- ggplot(data=rarecurve_df, aes(x=n_seqs, y=value, group=sample)) +
#   geom_vline(xintercept=2097, color="gray") +
#   geom_line() +
#   theme_classic()

##### Coverage calculations ####

goods_calc <- shared %>%
  group_by(sample) %>%
  summarize(n_seqs=sum(value),
            n_sings = sum(value==1),
            goods = 100*(1- n_sings/n_seqs)) %>% #number of singletons in each sample
  arrange(goods) #coverage values are all pretty mediocre; many singletons

coverage_stats <- goods_calc %>%
  filter(n_seqs > 2097)
#ggplot(aes(x=n_seqs, y=goods)) + geom_point()

coverage_stats %>%
  arrange(goods)
#for soil samples, good's above ~70% is standard
#reference: https://forum.mothur.org/t/alpha-diversity-after-remove-low-abundance/2619/3

#calculate Good's coverage
#Goods = number of OTUs with 1 sequence that have been seen once
#Goods = 100(1-(n(i)/N total))
