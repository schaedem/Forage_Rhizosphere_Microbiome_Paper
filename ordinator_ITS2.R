library(tidyverse)
library(vegan)
library(easyCODA)
#load ordinator
load("pcwOrd_0.1.Rdata")

shared <- readRDS("final_shared_rare_ITS2.rds")
k_shared <- shared %>% filter(location == "Karama")
r_shared <- shared %>% filter(location == "Rubona")
b_shared <- shared  %>% filter(location == "Burera")

wlra_fun <- function(shared) {
  #Format matrix
  shared_wide <- shared %>%
    dplyr::select(-c(location, treatment, species, timepoint, block)) %>%
    pivot_wider(names_from = name, values_from = value, values_fill = 0) %>%
    column_to_rownames("sample")
  
  shared_mat <- as.matrix(shared_wide)
  
  #Remove all columns with only 0 values
  shared_mat <- shared_mat[, colSums(shared_mat != 0) > 0]
  
  #Log ratio analysis
  #1. remove 0s, set to 0.5 pseudocounts
  cc <- shared_mat
  cc[cc==0] <- 0.5 #set 0 to 0.5 for log transformation
  
  #2. close the community so that all abundances are relative abundances
  #and rowSums(x) = 1
  closed_comm <- sweep(cc, 1, rowSums(shared_mat), '/')
  
  #3. Perform centered log-ratio analysis (CLR): 
  #values are centered on the geometric mean
  geo_mean <- apply(closed_comm, 1, 
                    function(x){
                      exp(sum(log(x))/length(x))
                    })
  lr_comm <- sweep(closed_comm, 1, geo_mean, '/')
  lr_comm <-log(lr_comm, 2)
  
  #4. Weighted LRA - successful in explaining more variance than unweighted
  wlr_comm <- CLR(closed_comm) #weighted on column geometric mean 
  
  #4. Run PCA on the log-ratio transformed matrix
  w_lra <- pcwOrd(wlr_comm)
  #w_lra_plot <- plot_ord(w_lra)
  
  #5. Extract axis 1
  #samples
  row_scores <- ord_scores(w_lra, choice='row', scaling='principle', axes=c(1:5))
  #otus
  col_scores <- ord_scores(w_lra, choice="column", scaling='principle', axes=c(1:5) )
  
  return(list("row_scores" = row_scores, "col_scores" = col_scores, "wlra" = w_lra))
  
}

# #All data
# all_wlra <- wlra_fun(shared) #axis1 = 26.7, axis2=12.4
# allscores <- all_wlra$row_scores %>% rownames_to_column("sample")
# allscores_otu <- all_wlra$col_scores %>% rownames_to_column("name")
# 
#Burera
burera_wlra <- wlra_fun(b_shared) #axis1 = 15.1%, axis2=11.4%, axis3 = 10.4%
bscores <- burera_wlra$row_scores %>% rownames_to_column("sample")
head(bscores)
bscores_otu <- burera_wlra$col_scores %>% rownames_to_column("name")
plot_ord(burera_wlra$wlra ,axes=c(1:2))
ord_scree(burera_wlra$wlra)

#Karama
karama_wlra <- wlra_fun(k_shared) #axis1 = 16.7%, axis2 = 10.1%, axis3 = 8.1%
kscores <- karama_wlra$row_scores %>% rownames_to_column("sample")
kscores_otu <- karama_wlra$col_scores %>% rownames_to_column("name")
head(kscores)
plot_ord(karama_wlra$wlra ,axes=c(1:2))

#Rubona
rubona_wlra <- wlra_fun(r_shared) #axis1 = 19% , axis2 =9.99%, axis3 = 
rscores <- rubona_wlra$row_scores %>% rownames_to_column("sample")
rscores_otu <- karama_wlra$col_scores %>% rownames_to_column("name")
head(rscores)
plot_ord(rubona_wlra$wlra ,axes=c(1:2))

#save results
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Marie/Thesis/Rwanda/Lab work/04_Station_BNF/final_data/seq_data/ITS2/cla_data")
write_csv(rscores, "rscores.csv")
write_csv(rscores_otu, "rscores_otu.csv")
write_csv(kscores, "kscores.csv")
write_csv(kscores_otu, "kscores_otu.csv")
write_csv(bscores, "bscores.csv")
write_csv(bscores_otu, "bscores_otu.csv")
# write_csv(allscores, "allscores.csv")
# write_csv(allscores_otu, "allscores_otu.csv")
