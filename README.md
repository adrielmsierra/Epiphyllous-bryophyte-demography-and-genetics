# Epiphyllous_bryophyte_demography_and_genetics
Pipeline for the analysis of population demography and genetics of epiphyllous bryophyte in a fragmented landscape

- Population density changes from long-term censuses of two epiphyllous bryophyte species
- Estimate population genetic diversity of bryophyte species in small and large patches
- Infer population genetic structure with individual-based clustering and pairwise comparisons
- Quantify the proportion and direction of recent migration patterns among populations in small and large patches

# DATA

Genind object with assembled Single Nucleotide Polymorphisms (SNPs) by individuals in 12 populations

- Dataset filtered the loci present in at least 15% of the individuals (filter parameter: -R=0.15), with sample size of n=105
Radula_105R15_genind.rds 
- Dataset filtered the loci present in at least 20% of the individuals (filter parameter: -R=0.20), with sample sizes of n=80
Radula_80R20_genind.rds
- Dataset filtered the loci present in at least 15% of the individuals (filter parameter: -R=0.15), with sample sizes of n=107
Cololejeunea_107R15_genind.rds
- Dataset filtered the loci present in at least 20% of the individuals (filter parameter: -R=0.20), with sample sizes of n=71
Cololejeunea_71R20_genind.rds
- Raw data of the number of colonies was estimated for each study plot
Demographic dataset.csv
- Genetic diversity parameters compiled by population
Gen div Stats.csv

# Scripts in this repository are organized in the following directories

- [Estimates of epiphyllous population size dynamics](https://github.com/adrielmsierra/Epiphyllous-bryophyte-demography-and-genetics/tree/main/Estimates%20of%20epiphyllous%20population%20size%20dynamics)
- [Genetic summary statistics](https://github.com/adrielmsierra/Epiphyllous-bryophyte-demography-and-genetics/tree/main/Genetic%20summary%20statistics)
- [Population structure with individual-based clustering](https://github.com/adrielmsierra/Epiphyllous-bryophyte-demography-and-genetics/tree/main/Population%20structure%20with%20individual-based%20clustering)
- [Spatial autocorrelation and migration patterns](https://github.com/adrielmsierra/Epiphyllous-bryophyte-demography-and-genetics/tree/main/Spatial%20autocorrelation%20and%20migration%20patterns)
