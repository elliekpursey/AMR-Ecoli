library(tidyverse)
library(data.table)

# set path
root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

# amr df tidy
amrfinder_cols = c('filename', 'prot_identifier', 'contig_id', 
                   'start', 'stop', 'strand', 'gene_symbol', 
                   'sequence_name', 'scope', 'element_type',
                   'element_subtype', 'class', 'subclass',
                   'method', 'target_length', 'ref_seq_length',
                   "percent_cov", "percent_id", "align_length",
                   "accession_closest_seq", "name_closest_seq", 'HMM_id',
                   'HMM_description')

amrfinder <- fread("amrfinder_all.csv", sep=',', col.names=amrfinder_cols)

amrfinder_tidy <- amrfinder %>%
  separate(col=filename, into=c("results", "amrfinder", "genome"), sep="/") %>%
  separate(col="genome", into=c("GCF", "genome"), sep="_") %>%
  unite(refseq_id, GCF, genome, sep = "_", remove = FALSE) %>%
  dplyr::select(refseq_id, gene_symbol, sequence_name, element_type, element_subtype, class, subclass) %>%
  filter(element_type == "AMR" & element_subtype =="AMR")

write_csv(amrfinder_tidy, "acquired_amr.csv")

## raw counts of each gene type, class etc.

gene_counts <- amrfinder_tidy %>%
  count(gene_symbol)

class_counts <- amrfinder_tidy %>%
  count(class)

subclass_counts <- amrfinder_tidy %>%
  count(subclass)

# sample size count

total_genomes <- amrfinder_tidy %>%
  dplyr::select(refseq_id) %>%
  unique() %>%
  count()
