# The distribution of antimicrobial resistance genes across phylogroup, host species and geography in 16,000 publicly-available E. coli genomes

This repository contains the Snakemake pipeline used to [phylogroup](https://github.com/nickp60/EzClermont) and [detect mobile ARGs](https://github.com/ncbi/amr) in _E. coli_ genomes for the above titled manuscript.

# Snakemake pipeline structure

The pipeline is situated in workflow > Snakefile_amrfinder. This directory also includes all conda environment specifications and scripts needed for it to run. Input genomes were placed under resources > nuc_fastas as nucleotide fasta files - the resources directory also contains images used to create plots. The config file allows specification of the temp directory location.  

# Script functions

This table describes all scripts stored in the "other_scripts" directory used for analysis downstream of the Snakemake pipeline, in the order in which they were run.

|Script directory |   Function    |
|-----------------|---------------|
| data_generation | scripts to retrieve, then classify biosample data (note classifications are based on manual examination of data in manuscript)|
| data_tidying    | scripts to tidy the AMR data and metadata, then generate the dataset used for modelling |
| analysis_scripts| scripts to extract subsamples, perform descriptive analysis and all modelling and finally generate main manuscript figures and supplementary tables |
| functions       | functions called for modelling within above scripts |

# Reference 

Up-to-date preprint here
