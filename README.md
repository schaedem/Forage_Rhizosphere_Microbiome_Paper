# Forage_Rhizosphere_Microbiome_Paper
This repository contains: 1) scripts that were used to process raw amplicon (16S, ITS2) data using mothur and 2) R scripts that were used to tidy and analyze data from sequencing, soil N concentrations, and plant biomass measurements.

## 16S rRNA Amplicon Sequence Data

### final_shared_rare_16S.rds
R data file containing rarefied (seqs=2097) 16S data that has been merged with appropriate metadata (location, treatment, species, timepoint, block). This is the rarefied data that has been used for all downstream analyses and diversity measures.

### rarefied_richness_final_16S.rds
R data file containing rarefied (seqs=2097) 16S species richness data for each sample (n=427).

### goods_coverage_rare_16S.rds
R data file containing Good's Coverage values for rarefied 16S sequence data.

### tidy_taxonomy_16S.rds
R data file containing parsed taxonomy key information for 16S OTUs. Separate columns are given for kingdom, phylum, class, order, family, and genus.

### otu_rare_relabund_16S.rds
R data file containing rarefied OTU abundance data with associated metadata, parsed taxonomic information, and relative abundance. Refer to file "tidy_taxonomy.R" for relative abundance calculation.

## ITS2 rRNA Amplicon Sequence Data

### final_shared_rare_ITS2.rds
R data file containing rarefied (seqs=1923) ITS2 data that has been merged with appropriate metadata (location, treatment, species, timepoint, block). This is the rarefied data that has been used for all downstream analyses and diversity measures.

### full_rare_rel_abund_ITS2.rds
R data file containing rarefied OTU abundance data with associated metadata, parsed taxonomic information, and relative abundance. Refer to file "tidy_taxonomy.R" for relative abundance calculation.

### goods_coverage_rare_ITS2.rds
R data file containing parsed taxonomy key information for ITS2 OTUs. Separate columns are given for kingdom, phylum, class, order, family, and genus.

### rare_its2_guild.rds
R data file containing FUNGuild classification results for ITS2 OTUs and associated sample metadata. Used in downstream analyses involving fungal trophic modes. 

### rarefied_richness_final_ITS2.rds
R data file containing rarefied (seqs=2097) 16S species richness data for each sample (n=427).

### tidy_taxonomy_ITS2.rds
R data file containing parsed taxonomy key information for ITS2 OTUs. Separate columns are given for kingdom, phylum, class, order, family, and genus.

### funguild_all.txt
FunGuild classification results for all ITS2 OTUs
