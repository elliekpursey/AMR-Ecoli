#!/usr/bin/env RScript

library(tidyverse)
library(data.table)
library(glmmTMB)
library(DHARMa)
library(ggpubr)
library(ggeffects)

# set workspace

root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

###############################################################################

# read in raw data

ceph <- fread("CTX-M_binomial.csv", sep=',')

# analysis of proportions in raw data

raw_proportions <- ceph %>%
  group_by(phylogroup) %>%
  mutate(total = sum(ceph)) %>%
  mutate(phylo_total = n()) %>%
  mutate(prop = total / n()) %>%
  ungroup() %>%
  dplyr::select(phylogroup, prop, total,phylo_total) %>%
  unique()

# read in function 

source('../other_scripts_local/functions/model_fit_and_test.R')

############################### PHYLOGROUPS ####################################

# read in subsampled/resampled data and merge with ceph data to get relevant
# subsampled dataframes

phylogroup_min <- fread("subsampled_dataframes/subsampled_phylogroups_no_rep.csv", sep=",") %>%
  left_join(., ceph)

phylogroup_500 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_500.csv", sep=",") %>%
  left_join(., ceph)

phylogroup_1000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_1000.csv", sep=",") %>%
  left_join(., ceph)

phylogroup_2000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_2000.csv", sep=",") %>%
  left_join(., ceph)

phylogroup_4000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_4000.csv", sep=",") %>%
  left_join(., ceph)

# model without subsampling

# get list of phylogroups removed for subsampling and use full dataframe with 
# these groups removed for consistency

list_of_phylogroups <- phylogroup_min %>%
  dplyr::select(phylogroup) %>%
  unique()

phylogroup_model_full_df <- ceph %>%
  dplyr::filter(phylogroup %in% list_of_phylogroups$phylogroup) 

# run model

phylogroup_ceph_no_subsampling <- glmmTMB(ceph ~ phylogroup, data=phylogroup_model_full_df, na.action = "na.fail", family = "binomial")

saveRDS(phylogroup_ceph_no_subsampling, 'model_objects/phylogroup_ceph_no_subsampling')

phylogroup_ceph_no_subsampling_summary <- summary(phylogroup_ceph_no_subsampling)

# test fit 
testDispersion(phylogroup_ceph_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = phylogroup_ceph_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, phylogroup_model_full_df$phylogroup, quantreg = T)

# extract model coefficients

phylogroup_ceph_no_subsampling_model_coeff <- phylogroup_ceph_no_subsampling_summary$coefficients$cond
write.csv(phylogroup_ceph_no_subsampling_model_coeff, file = 'tables/phylogroup_ceph_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

phylogroup_ceph_no_subsampling_pred_df <- ggpredict(phylogroup_ceph_no_subsampling, terms = c('phylogroup'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

phylogroup_ceph_no_subsampling_pred_plot <- ggplot(phylogroup_ceph_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Prob. of CTX-M")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(phylogroup_min, ceph, phylogroup, "Phylogroup", "Prob. of CTX-M")
model_fit_and_test(phylogroup_500, ceph, phylogroup, "Phylogroup", "Prob. of CTX-M")
model_fit_and_test(phylogroup_1000, ceph, phylogroup, "Phylogroup", "Prob. of CTX-M")
model_fit_and_test(phylogroup_2000, ceph, phylogroup, "Phylogroup", "Prob. of CTX-M")
model_fit_and_test(phylogroup_4000, ceph, phylogroup, "Phylogroup", "Prob. of CTX-M")

# arrange plots

arranged_ceph_subsamples_phylogroup <- ggarrange(phylogroup_ceph_no_subsampling_pred_plot, phylogroup_min_ceph, phylogroup_500_ceph,
                                                phylogroup_1000_ceph, phylogroup_2000_ceph, phylogroup_4000_ceph,
                                                labels = c('A', 'B', 'C', 'D', 'E', 'F'))

############################### HOSTS ####################################

# read in subsampled/resampled data and merge with ceph data to get relevant
# subsampled dataframes

host_min <- fread("subsampled_dataframes/subsampled_hosts_no_rep.csv", sep=",") %>%
  left_join(., ceph)

host_500 <- fread("subsampled_dataframes/host_resample_with_rep_500.csv", sep=",") %>%
  left_join(., ceph)

host_1000 <- fread("subsampled_dataframes/host_resample_with_rep_1000.csv", sep=",") %>%
  left_join(., ceph)

host_2000 <- fread("subsampled_dataframes/host_resample_with_rep_2000.csv", sep=",") %>%
  left_join(., ceph)

host_4000 <- fread("subsampled_dataframes/host_resample_with_rep_4000.csv", sep=",") %>%
  left_join(., ceph)

# model without subsampling

# run model

host_ceph_no_subsampling <- glmmTMB(ceph ~ host, data=ceph, na.action = "na.fail", family = "binomial")

saveRDS(host_ceph_no_subsampling, 'model_objects/host_ceph_no_subsampling')

host_ceph_no_subsampling_summary <- summary(host_ceph_no_subsampling)

# test fit 
testDispersion(host_ceph_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = host_ceph_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, ceph$host, quantreg = T)

# extract model coefficients

host_ceph_no_subsampling_model_coeff <- host_ceph_no_subsampling_summary$coefficients$cond
write.csv(host_ceph_no_subsampling_model_coeff, file = 'tables/host_ceph_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

host_ceph_no_subsampling_pred_df <- ggpredict(host_ceph_no_subsampling, terms = c('host'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

host_ceph_no_subsampling_pred_plot <- ggplot(host_ceph_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Host", y = "Prob. of CTX-M")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(host_min, ceph, host, "Host", "Prob. of CTX-M")
model_fit_and_test(host_500, ceph, host, "Host", "Prob. of CTX-M")
model_fit_and_test(host_1000, ceph, host, "Host", "Prob. of CTX-M")
model_fit_and_test(host_2000, ceph, host, "Host", "Prob. of CTX-M")
model_fit_and_test(host_4000, ceph, host, "Host", "Prob. of CTX-M")

# arrange plots

arranged_ceph_subsamples_host <- ggarrange(host_ceph_no_subsampling_pred_plot, host_min_ceph, host_500_ceph,
                                                 host_1000_ceph, host_2000_ceph, host_4000_ceph,
                                                 labels = c('A', 'B', 'C', 'D', 'E', 'F'))

############################### SUBREGIONS ####################################

# read in subsampled/resampled data and merge with ceph data to get relevant
# subsampled dataframes

subregion_min <- fread("subsampled_dataframes/subsampled_subregions_no_rep.csv", sep=",") %>%
  left_join(., ceph)

subregion_500 <- fread("subsampled_dataframes/subregion_resample_with_rep_500.csv", sep=",") %>%
  left_join(., ceph)

subregion_1000 <- fread("subsampled_dataframes/subregion_resample_with_rep_1000.csv", sep=",") %>%
  left_join(., ceph)

subregion_2000 <- fread("subsampled_dataframes/subregion_resample_with_rep_2000.csv", sep=",") %>%
  left_join(., ceph)

subregion_4000 <- fread("subsampled_dataframes/subregion_resample_with_rep_4000.csv", sep=",") %>%
  left_join(., ceph)

# model without subsampling

# get list of subregions removed for subsampling and use full dataframe with 
# these groups removed for consistency

list_of_subregions <- subregion_min %>%
  dplyr::select(subregion) %>%
  unique()

subregion_model_full_df <- ceph %>%
  dplyr::filter(subregion %in% list_of_subregions$subregion) 

# run model

subregion_ceph_no_subsampling <- glmmTMB(ceph ~ subregion, data=subregion_model_full_df, na.action = "na.fail", family = "binomial")

saveRDS(subregion_ceph_no_subsampling, 'model_objects/subregion_ceph_no_subsampling')

subregion_ceph_no_subsampling_summary <- summary(subregion_ceph_no_subsampling)

# test fit 
testDispersion(subregion_ceph_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = subregion_ceph_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, subregion_model_full_df$subregion, quantreg = T)

# extract model coefficients

subregion_ceph_no_subsampling_model_coeff <- subregion_ceph_no_subsampling_summary$coefficients$cond
write.csv(subregion_ceph_no_subsampling_model_coeff, file = 'tables/subregion_ceph_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

subregion_ceph_no_subsampling_pred_df <- ggpredict(subregion_ceph_no_subsampling, terms = c('subregion'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

subregion_ceph_no_subsampling_pred_plot <- ggplot(subregion_ceph_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Subregion", y = "Prob. of CTX-M")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(subregion_min, ceph, subregion, "Subregion", "Prob. of CTX-M")
model_fit_and_test(subregion_500, ceph, subregion, "Subregion", "Prob. of CTX-M")
model_fit_and_test(subregion_1000, ceph, subregion, "Subregion", "Prob. of CTX-M")
model_fit_and_test(subregion_2000, ceph, subregion, "Subregion", "Prob. of CTX-M")
model_fit_and_test(subregion_4000, ceph, subregion, "Subregion", "Prob. of CTX-M")

# arrange plots

arranged_ceph_subsamples_subregion <- ggarrange(subregion_ceph_no_subsampling_pred_plot, subregion_min_ceph, subregion_500_ceph,
                                                 subregion_1000_ceph, subregion_2000_ceph, subregion_4000_ceph,
                                                 labels = c('A', 'B', 'C', 'D', 'E', 'F'))

############ SAVE ALL PLOTS ############

ggsave(plot = arranged_ceph_subsamples_phylogroup, "plots/CTXM_phylogroup_subsampled_pred_plots.jpg", width=30, height=20, units="cm")
ggsave(plot = arranged_ceph_subsamples_host, "plots/CTXM_host_subsampled_pred_plots.jpg", width=30, height=20, units="cm")
ggsave(plot = arranged_ceph_subsamples_subregion, "plots/CTXM_subregion_subsampled_pred_plots.jpg", width=30, height=30, units="cm")

ggsave(plot = arranged_ceph_subsamples_phylogroup, "plots/CTXM_phylogroup_subsampled_pred_plots.jpg", width=30, height=20, units="cm")
ggsave(plot = arranged_ceph_subsamples_host, "plots/CTXM_host_subsampled_pred_plots.jpg", width=30, height=20, units="cm")
ggsave(plot = arranged_ceph_subsamples_subregion, "plots/CTXM_subregion_subsampled_pred_plots.jpg", width=30, height=30, units="cm")


