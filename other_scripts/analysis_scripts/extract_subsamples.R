#!/usr/bin/env RScript

library(tidyverse)
library(data.table)

# set workspace

root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

# read in metadata 

genomes_and_metadata <- fread("AMR_counts_and_metadata.csv", sep=',') 

### get minimum subsample sizes ###

## phylogroup

# sample sizes

phylogroup_sample_sizes <- genomes_and_metadata %>%
  dplyr::select(phylogroup) %>%
  count(phylogroup)

# remove groups with small sample sizes - (below 100)

phylogroups_to_keep <- phylogroup_sample_sizes %>%
  dplyr::filter(n >= 100) 

phylogroups_filtered <- genomes_and_metadata %>%
  dplyr::filter(phylogroup %in% phylogroups_to_keep$phylogroup)

# get minimum subsample size (without replacement)

phylogroup_min_sample_size <- min(phylogroups_to_keep$n)

## host

# sample sizes

host_sample_sizes <- genomes_and_metadata %>%
  dplyr::select(host) %>%
  count(host)

# get minimum subsample size (without replacement)

host_min_sample_size <- min(host_sample_sizes$n)

## subregion

# sample sizes

subregion_sample_sizes <- genomes_and_metadata %>%
  dplyr::select(subregion) %>%
  count(subregion)

# remove groups with small sample sizes - (below 100)

subregions_to_keep <- subregion_sample_sizes %>%
  dplyr::filter(n >= 100) 

subregions_filtered <- genomes_and_metadata %>%
  dplyr::filter(subregion %in% subregions_to_keep$subregion)

# get minimum subsample size (without replacement)

subregion_min_sample_size <- min(subregions_to_keep$n)


#################### Do subsampling/resampling #################################

# function to subsample all groups to size of smallest group without replacement

subsample_no_replacement <- function(subsample_size, grouping_var, in_df) {
  
  subsample_df <- in_df %>%
    group_by({{grouping_var}}) %>%
    do(sample_n(., subsample_size, replace = FALSE)) %>%
    dplyr::select({{grouping_var}}, refseq_id)

    return(subsample_df)

  }

# run function for host, phylogroup and subregion

wo_rep_subsample_phylogroup <- subsample_no_replacement(phylogroup_min_sample_size, phylogroup, phylogroups_filtered)
wo_rep_subsample_host <- subsample_no_replacement(host_min_sample_size, host, genomes_and_metadata)
wo_rep_subsample_subregion <- subsample_no_replacement(subregion_min_sample_size, subregion, subregions_filtered)

# write out subsampled dfs

write.csv(wo_rep_subsample_phylogroup, file = "subsampled_dataframes/subsampled_phylogroups_no_rep.csv", row.names = FALSE)
write.csv(wo_rep_subsample_host, file = "subsampled_dataframes/subsampled_hosts_no_rep.csv", row.names = FALSE)
write.csv(wo_rep_subsample_subregion, file = "subsampled_dataframes/subsampled_subregions_no_rep.csv", row.names = FALSE)

# function to resample all groups with replacement from subsample sizes of 100-4000 (spans sample sizes of all groups)

subsample_with_replacement <- function(grouping_var, in_df) {
  
  subsample_sizes <- c(100,500,1000,2000,4000)
  
  for(sample_size in subsample_sizes){
  
    subsample_df <- in_df %>%
    group_by({{grouping_var}}) %>%
    do(sample_n(., sample_size, replace = TRUE)) %>%
    dplyr::select({{grouping_var}}, refseq_id) %>%
    ungroup()
    
    df_name <- paste0(deparse(substitute(grouping_var)), '_resample_with_rep_', as.factor(sample_size))
    
    assign(df_name, data.frame(subsample_df), envir = .GlobalEnv)
    
    file_name <- paste0('subsampled_dataframes/', df_name, '.csv')
    write.csv(subsample_df, file_name, row.names=FALSE)
  }
  
}

# use function to subsample host, phylogroup and subregion with replacement for
# various subsample sizes - use same filtered dfs for this, note this writes out
# .csv files directly as well as assigning to global environment

subsample_with_replacement(phylogroup, phylogroups_filtered)
subsample_with_replacement(subregion, subregions_filtered)
subsample_with_replacement(host, genomes_and_metadata)



