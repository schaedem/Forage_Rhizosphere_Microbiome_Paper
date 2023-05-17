# Forage_Rhizosphere_Microbiome_Paper
This repository contains: 1) amplicon sequence data (16S, ITS2) processed using mothur from forage rhizosphere soil samples, 2) forage biomass data, 3) bulk soil mineral N data, and 4) R scripts that were used to tidy and analyze all datasets.

The goal of this project was to:
  1. Evaluate the impact of perennial legume intercropping on biological nitrogen fixation rates
  2. Determine whether a perennial intercropping arrangement altered root-associated microbial communities
  3. Do intercropping-driven changes to community composition have implications for forage nutrive quality?


  The data contained in this repository derives from soil (bulk and rhizosphere) and plant leaf tissue sampling
  from replicated trials in Rwanda (Nyagatare (Karama research station), Nyanza (Rubona research station), Burera (Rwerere Research Station)).

  Field trials in *Rubona and Nyagatare/Karama* were sampled for a total of 4 timepoints
  approximately every 8 weeks at forage anthesis:

    T1 (forage anthesis, late rainy season) : March 2021
    T2 (forage anthesis, early dry season) : May/June 2021
    T3 (forage anthesis, short rainy season) : August/September 2021
    T4 (forage anthesis, dry season) : November 2021

  Field trials in *Burera/Rwerere* were sampled for a total of 3 timepoints
  approximately every 8 weeks at forage anthesis:

    T1 (forage anthesis, late rainy season) : March/April 2021
    T2 (forage anthesis, dry-ish season) : June/July 2021
    T3 (forage anthesis, dry season/short rainy season): October 2021

  Treatments (7), 4 replicates:
  Maize monoculture, maize + Desmodium distortum intercrop
  Napier grass monoculture, Napier + D. distortum intercrop
  Brachiaria cv. Mulato II monoculture, Brachiaria + D. distortum intercrop
  D. distortum monocrop

  Labeling scheme (BULK SOIL) :
    T1 samples have ID numbers in the 100s.
    T2 samples have ID numbers in the 200s.
    T3 samples have ID numbers in the 300s.
    T4 samples have ID numbers in the 400s.

  Labeling scheme (RHIZOSPHERE SOIL) :
    T1 rhizosphere samples: 1R - 120R
    T2 rhizosphere samples: 201R - 320R
    T3 rhizosphere samples: 401-520R
    T4 rhizosphere samples: 601-680 (Burera not sampled)

  Labeling scheme (LEAF BIOMASS) :
    T1 leaf samples: 1L - 120L
    T2 leaf samples: 201L - 320L
    T3 leaf samples: 401L-520L
    T4 leaf samples: 601L - 680L (Burera not sampled)

  Laboratory procedures (bulk soil) :

  Fresh soil was shipped in insulated containers with icepacks and was used
  to determin total extractable mineral (NO3-N + NH4-N) nitrogen.

  Laboratory procedures (rhizosphere soil) :

  Forages, including maize, were excavated and shaken to dislodge loosely bound soil.
  Rhizosphere soil was collected into plastic bags using sterile toothbrushes.
  Fresh rhizosphere soil was shipped in insulated containers with icepacks.
  Upon arrival, it was sieved and stored at -20C prior to DNA extraction.
  DNA extracts were submitted to UMN Genomics Center for ITS2 and 16S V4 amplicon sequencing.

  Laboratory procedures (leaf biomass) :

  Total aboveground biomass (same plant as used for rhizosphere sampling) was harvested and leaf tissue was separated from the stem.
  Leaf biomass was dried at 90C for 48 hours in Rwanda, then passed through a food processor to obtain a course grind.
  Upon arrival at the UMN, forage biomass was further dried for 12 hours at 90C prior to 1mm grinding with a Foss Biomass Grinder.
  1mm ground samples were used for NIR scanning (protein, lignin, ADF, NDF).
  Samples were additionally ball-ground with a GenoGrinder prior to microscaling for EA-IRMS analysis.


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
