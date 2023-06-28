#!/usr/bin/env RScript

library(tidyverse)
library(data.table)

# set path

root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

# read in raw dataframes

metadata_phylogroups <- fread("metadata_phylogroups.csv", sep=',') %>%
  drop_na()

acq_amr <- fread("acquired_amr.csv", sep=',') 

# full datatset with all AMR genes and metadata

merge_df <- metadata_phylogroups %>%
  left_join(., acq_amr) %>%
  mutate(host = str_replace(host, "_", " ")) 

write.csv(merge_df, file = "AMR_and_metadata.csv", row.names = FALSE)

# get counts of ARGs and classes

gene_counts <- merge_df %>%
  count(gene_symbol)

class_counts <- merge_df %>%
  count(class)

subclass_counts <- merge_df %>%
  count(subclass)

# neat temp dataset for merging to get full genome and metadata list

genomes_df <- merge_df %>%
  dplyr::select(refseq_id, subregion, phylogroup, host) %>%
  unique()

### creation of count, MDR/XDR, CTX-M and carbapenemase datasets

## AMR gene counts, with removal of ubiquitous classes/genes

count_df <- merge_df %>%
  dplyr::filter(class != "EFFLUX" &
                  gene_symbol != "blaEC") %>%
  count(refseq_id, host, subregion, phylogroup, name="count_amr") %>%
  full_join(., genomes_df) %>%
  replace(is.na(.), 0) 

write.csv(count_df, file = "AMR_counts_and_metadata.csv", row.names = FALSE)

## binomial MDR/XDR dataset

count_classes <- merge_df %>%
  filter(gene_symbol != "blaEC" &
           class != "EFFLUX") %>% # remove genes that are intrinsic
  count(refseq_id, class, subregion, phylogroup, host, name="count") %>% # count numbers in each class
  mutate(count, count = ifelse(count > 0, 1, 0)) %>% # convert counts to binomial data
  pivot_wider(names_from=class, values_from=count) %>% # get column for each class
  replace(is.na(.), 0) %>% 
  rowwise() %>% 
  mutate(count_classes = sum(c_across(c(5:22)), na.rm = T)) %>% # sum number of resistant classes
  ungroup() %>%
  full_join(., genomes_df) %>% # join to full dataset to include genomes with no resistance
  replace(is.na(.), 0) 

# assign categories 

MDR_df <- count_classes %>%
  mutate(MDR_binom = ifelse(count_classes > 2, 1, 0)) %>%
  mutate(XDR_binom = ifelse(count_classes > 9, 1, 0)) 

write.csv(MDR_df, file = "MDR_XDR_binomial.csv", row.names = FALSE)

## CTX-M binomial dataset

# get list of different CTX-M genes

ctxm_genes <- acq_amr %>%
  dplyr::filter(grepl('CTX-M', gene_symbol)) %>%
  dplyr::select(gene_symbol) %>%
  count(gene_symbol)

# make presence/absence df for CTX-M resistance genes, merge with df of
# all genomes and metadata

binomial_ceph_resistance <- acq_amr %>%
  dplyr::filter(grepl('CTX-M', gene_symbol)) %>%
  mutate(ceph = 1) %>%
  dplyr::select(refseq_id, ceph) %>%
  unique() %>%
  right_join(., genomes_df) %>%
  replace_na(list(ceph = 0)) 

write.csv(binomial_ceph_resistance, file = "CTX-M_binomial.csv", row.names = FALSE)

## carbapenemase binomial dataset

# get list of different CTX-M genes

carb_genes <- acq_amr %>%
  dplyr::filter(subclass == "CARBAPENEM") %>%
  dplyr::select(gene_symbol) %>%
  count(gene_symbol)

# make presence/absence df for carbapenem resistance genes, merge with df of
# all genomes and metadata

binomial_carb_resistance <- acq_amr %>%
  dplyr::filter(subclass == "CARBAPENEM") %>%
  mutate(carb = 1) %>%
  dplyr::select(refseq_id, carb) %>%
  unique() %>%
  right_join(., genomes_df) %>%
  replace_na(list(carb = 0)) 

write.csv(binomial_carb_resistance, file = "carb_binomial.csv", row.names = FALSE)
