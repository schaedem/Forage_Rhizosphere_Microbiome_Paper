# Forage_Rhizosphere_Microbiome_Paper
This repository contains: 1) scripts that were used to process raw amplicon (16S, ITS2) data using mothur and 2) R scripts that were used to tidy and analyze data from sequencing, soil N concentrations, and plant biomass measurements.

## final_shared_rare_16S.rds
R data file containing rarefied (seqs=2097) 16S data that has been merged with appropriate metadata (location, treatment, species, timepoint, block). This is the rarefied data that has been used for all downstream analyses and diversity measures.

## rarefied_richness_final_16S.rds
R data file containing rarefied (seqs=2097) 16S species richness data for each sample (n=427).

## shared_full_16S_final.rds
R data file containing *non*rarefied sequence data that has been merged with appropriate metadata.
