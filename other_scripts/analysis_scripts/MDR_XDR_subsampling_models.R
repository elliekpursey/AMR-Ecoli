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

MDR_XDR <- fread("MDR_XDR_binomial.csv", sep=",")

# analysis of proportions in raw data

raw_proportions <- MDR_XDR %>%
  group_by(phylogroup) %>%
  mutate(total = sum(MDR_binom)) %>%
  mutate(phylo_total = n()) %>%
  mutate(prop = total / n()) %>%
  ungroup() %>%
  dplyr::select(phylogroup, prop, total,phylo_total) %>%
  unique()

# read in function 

source('../other_scripts_local/functions/model_fit_and_test.R')

############################### PHYLOGROUPS ####################################

# read in subsampled/resampled data and merge with MDR/XDR data to get relevant
# subsampled dataframes

phylogroup_min <- fread("subsampled_dataframes/subsampled_phylogroups_no_rep.csv", sep=",") %>%
  left_join(., MDR_XDR)

phylogroup_500 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_500.csv", sep=",") %>%
  left_join(., MDR_XDR)

phylogroup_1000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_1000.csv", sep=",") %>%
  left_join(., MDR_XDR)

phylogroup_2000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_2000.csv", sep=",") %>%
  left_join(., MDR_XDR)

phylogroup_4000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_4000.csv", sep=",") %>%
  left_join(., MDR_XDR)

##### MDR #####

# model without subsampling

# get list of phylogroups removed for subsampling and use full dataframe with 
# these groups removed for consistency

list_of_phylogroups <- phylogroup_min %>%
  dplyr::select(phylogroup) %>%
  unique()

phylogroup_model_full_df <- MDR_XDR %>%
  dplyr::filter(phylogroup %in% list_of_phylogroups$phylogroup) 
  
# run model

phylogroup_MDR_no_subsampling <- glmmTMB(MDR_binom ~ phylogroup, data=phylogroup_model_full_df, na.action = "na.fail", family = "binomial")

saveRDS(phylogroup_MDR_no_subsampling, 'model_objects/phylogroup_MDR_no_subsampling')

phylogroup_MDR_no_subsampling_summary <- summary(phylogroup_MDR_no_subsampling)

# test fit 
testDispersion(phylogroup_MDR_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = phylogroup_MDR_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, phylogroup_model_full_df$phylogroup, quantreg = T)

# extract model coefficients

phylogroup_MDR_no_subsampling_model_coeff <- phylogroup_MDR_no_subsampling_summary$coefficients$cond 
write.csv(phylogroup_MDR_no_subsampling_model_coeff, file = 'tables/phylogroup_MDR_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

phylogroup_MDR_no_subsampling_pred_df <- ggpredict(phylogroup_MDR_no_subsampling, terms = c('phylogroup'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

phylogroup_MDR_no_subsampling_pred_plot <- ggplot(phylogroup_MDR_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Prob. of MDR")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(phylogroup_min, MDR_binom, phylogroup, "Phylogroup", "Prob. of MDR")
model_fit_and_test(phylogroup_500, MDR_binom, phylogroup, "Phylogroup", "Prob. of MDR")
model_fit_and_test(phylogroup_1000, MDR_binom, phylogroup, "Phylogroup", "Prob. of MDR")
model_fit_and_test(phylogroup_2000, MDR_binom, phylogroup, "Phylogroup", "Prob. of MDR")
model_fit_and_test(phylogroup_4000, MDR_binom, phylogroup, "Phylogroup", "Prob. of MDR")

# arrange plots

arranged_MDR_subsamples_phylogroup <- ggarrange(phylogroup_MDR_no_subsampling_pred_plot, phylogroup_min_MDR_binom, phylogroup_500_MDR_binom,
                                     phylogroup_1000_MDR_binom, phylogroup_2000_MDR_binom, phylogroup_4000_MDR_binom,
                                     labels = c('A', 'B', 'C', 'D', 'E', 'F'))

##### XDR #####

# model without subsampling

# run model

phylogroup_XDR_no_subsampling <- glmmTMB(XDR_binom ~ phylogroup, data=phylogroup_model_full_df, na.action = "na.fail", family = "binomial")

saveRDS(phylogroup_XDR_no_subsampling, 'model_objects/phylogroup_XDR_no_subsampling')

phylogroup_XDR_no_subsampling_summary <- summary(phylogroup_XDR_no_subsampling)

# test fit 
testDispersion(phylogroup_XDR_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = phylogroup_XDR_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, phylogroup_model_full_df$phylogroup, quantreg = T)

# extract model coefficients

phylogroup_XDR_no_subsampling_model_coeff <- phylogroup_XDR_no_subsampling_summary$coefficients$cond
write.csv(phylogroup_XDR_no_subsampling_model_coeff, file = 'tables/phylogroup_XDR_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

phylogroup_XDR_no_subsampling_pred_df <- ggpredict(phylogroup_XDR_no_subsampling, terms = c('phylogroup'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

phylogroup_XDR_no_subsampling_pred_plot <- ggplot(phylogroup_XDR_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Prob. of XDR")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(phylogroup_min, XDR_binom, phylogroup, "Phylogroup", "Prob. of XDR")
model_fit_and_test(phylogroup_500, XDR_binom, phylogroup, "Phylogroup", "Prob. of XDR")
model_fit_and_test(phylogroup_1000, XDR_binom, phylogroup, "Phylogroup", "Prob. of XDR")
model_fit_and_test(phylogroup_2000, XDR_binom, phylogroup, "Phylogroup", "Prob. of XDR")
model_fit_and_test(phylogroup_4000, XDR_binom, phylogroup, "Phylogroup", "Prob. of XDR")

# arrange plots

arranged_XDR_subsamples_phylogroup <- ggarrange(phylogroup_XDR_no_subsampling_pred_plot, phylogroup_min_XDR_binom, phylogroup_500_XDR_binom,
                                     phylogroup_1000_XDR_binom, phylogroup_2000_XDR_binom, phylogroup_4000_XDR_binom,
                                     labels = c('A', 'B', 'C', 'D', 'E', 'F'))

############################### HOSTS ####################################

# read in subsampled/resampled data and merge with MDR/XDR data to get relevant
# subsampled dataframes

host_min <- fread("subsampled_dataframes/subsampled_hosts_no_rep.csv", sep=",") %>%
  left_join(., MDR_XDR)

host_500 <- fread("subsampled_dataframes/host_resample_with_rep_500.csv", sep=",") %>%
  left_join(., MDR_XDR)

host_1000 <- fread("subsampled_dataframes/host_resample_with_rep_1000.csv", sep=",") %>%
  left_join(., MDR_XDR)

host_2000 <- fread("subsampled_dataframes/host_resample_with_rep_2000.csv", sep=",") %>%
  left_join(., MDR_XDR)

host_4000 <- fread("subsampled_dataframes/host_resample_with_rep_4000.csv", sep=",") %>%
  left_join(., MDR_XDR)

##### MDR #####

# model without subsampling

# run model

host_MDR_no_subsampling <- glmmTMB(MDR_binom ~ host, data=MDR_XDR, na.action = "na.fail", family = "binomial")

saveRDS(host_MDR_no_subsampling, 'model_objects/host_MDR_no_subsampling')

host_MDR_no_subsampling_summary <- summary(host_MDR_no_subsampling)

# test fit 
testDispersion(host_MDR_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = host_MDR_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, MDR_XDR$host, quantreg = T)

# extract model coefficients

host_MDR_no_subsampling_model_coeff <- host_MDR_no_subsampling_summary$coefficients$cond
write.csv(host_MDR_no_subsampling_model_coeff, file = 'tables/host_MDR_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

host_MDR_no_subsampling_pred_df <- ggpredict(host_MDR_no_subsampling, terms = c('host'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

host_MDR_no_subsampling_pred_plot <- ggplot(host_MDR_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Host", y = "Prob. of MDR")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(host_min, MDR_binom, host, "Host", "Prob. of MDR")
model_fit_and_test(host_500, MDR_binom, host, "Host", "Prob. of MDR")
model_fit_and_test(host_1000, MDR_binom, host, "Host", "Prob. of MDR")
model_fit_and_test(host_2000, MDR_binom, host, "Host", "Prob. of MDR")
model_fit_and_test(host_4000, MDR_binom, host, "Host", "Prob. of MDR")

# arrange plots

arranged_MDR_subsamples_host <- ggarrange(host_MDR_no_subsampling_pred_plot, host_min_MDR_binom, host_500_MDR_binom,
                                     host_1000_MDR_binom, host_2000_MDR_binom, host_4000_MDR_binom,
                                     labels = c('A', 'B', 'C', 'D', 'E', 'F'))

##### XDR #####

# model without subsampling

# run model

host_XDR_no_subsampling <- glmmTMB(XDR_binom ~ host, data=MDR_XDR, na.action = "na.fail", family = "binomial")

saveRDS(host_XDR_no_subsampling, 'model_objects/host_XDR_no_subsampling')

host_XDR_no_subsampling_summary <- summary(host_XDR_no_subsampling)

# test fit 
testDispersion(host_XDR_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = host_XDR_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, MDR_XDR$host, quantreg = T)

# extract model coefficients

host_XDR_no_subsampling_model_coeff <- host_XDR_no_subsampling_summary$coefficients$cond
write.csv(host_XDR_no_subsampling_model_coeff, file = 'tables/host_XDR_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

host_XDR_no_subsampling_pred_df <- ggpredict(host_XDR_no_subsampling, terms = c('host'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

host_XDR_no_subsampling_pred_plot <- ggplot(host_XDR_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Host", y = "Prob. of XDR")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(host_min, XDR_binom, host, "Host", "Prob. of XDR")
model_fit_and_test(host_500, XDR_binom, host, "Host", "Prob. of XDR")
model_fit_and_test(host_1000, XDR_binom, host, "Host", "Prob. of XDR")
model_fit_and_test(host_2000, XDR_binom, host, "Host", "Prob. of XDR")
model_fit_and_test(host_4000, XDR_binom, host, "Host", "Prob. of XDR")

# arrange plots

arranged_XDR_subsamples_host <- ggarrange(host_XDR_no_subsampling_pred_plot, host_min_XDR_binom, host_500_XDR_binom,
                                     host_1000_XDR_binom, host_2000_XDR_binom, host_4000_XDR_binom,
                                     labels = c('A', 'B', 'C', 'D', 'E', 'F'))

############################### SUBREGIONS ####################################

# read in subsampled/resampled data and merge with MDR/XDR data to get relevant
# subsampled dataframes

subregion_min <- fread("subsampled_dataframes/subsampled_subregions_no_rep.csv", sep=",") %>%
  left_join(., MDR_XDR)

subregion_500 <- fread("subsampled_dataframes/subregion_resample_with_rep_500.csv", sep=",") %>%
  left_join(., MDR_XDR)

subregion_1000 <- fread("subsampled_dataframes/subregion_resample_with_rep_1000.csv", sep=",") %>%
  left_join(., MDR_XDR)

subregion_2000 <- fread("subsampled_dataframes/subregion_resample_with_rep_2000.csv", sep=",") %>%
  left_join(., MDR_XDR)

subregion_4000 <- fread("subsampled_dataframes/subregion_resample_with_rep_4000.csv", sep=",") %>%
  left_join(., MDR_XDR)

##### MDR #####

# model without subsampling

# get list of subregions removed for subsampling and use full dataframe with 
# these groups removed for consistency

list_of_subregions <- subregion_min %>%
  dplyr::select(subregion) %>%
  unique()

subregion_model_full_df <- MDR_XDR %>%
  dplyr::filter(subregion %in% list_of_subregions$subregion) 

# run model

subregion_MDR_no_subsampling <- glmmTMB(MDR_binom ~ subregion, data=subregion_model_full_df, na.action = "na.fail", family = "binomial")

saveRDS(subregion_MDR_no_subsampling, 'model_objects/subregion_MDR_no_subsampling')

subregion_MDR_no_subsampling_summary <- summary(subregion_MDR_no_subsampling)

# test fit 
testDispersion(subregion_MDR_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = subregion_MDR_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, subregion_model_full_df$subregion, quantreg = T)

# extract model coefficients

subregion_MDR_no_subsampling_model_coeff <- subregion_MDR_no_subsampling_summary$coefficients$cond
write.csv(subregion_MDR_no_subsampling_model_coeff, file = 'tables/subregion_MDR_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

subregion_MDR_no_subsampling_pred_df <- ggpredict(subregion_MDR_no_subsampling, terms = c('subregion'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

subregion_MDR_no_subsampling_pred_plot <- ggplot(subregion_MDR_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Subregion", y = "Prob. of MDR")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(subregion_min, MDR_binom, subregion, "Subregion", "Prob. of MDR")
model_fit_and_test(subregion_500, MDR_binom, subregion, "Subregion", "Prob. of MDR")
model_fit_and_test(subregion_1000, MDR_binom, subregion, "Subregion", "Prob. of MDR")
model_fit_and_test(subregion_2000, MDR_binom, subregion, "Subregion", "Prob. of MDR")
model_fit_and_test(subregion_4000, MDR_binom, subregion, "Subregion", "Prob. of MDR")

# arrange plots

arranged_MDR_subsamples_subregion <- ggarrange(subregion_MDR_no_subsampling_pred_plot, subregion_min_MDR_binom, subregion_500_MDR_binom,
                                          subregion_1000_MDR_binom, subregion_2000_MDR_binom, subregion_4000_MDR_binom,
                                          labels = c('A', 'B', 'C', 'D', 'E', 'F'))

##### XDR #####

# model without subsampling

# run model

subregion_XDR_no_subsampling <- glmmTMB(XDR_binom ~ subregion, data=subregion_model_full_df, na.action = "na.fail", family = "binomial")

saveRDS(subregion_XDR_no_subsampling, 'model_objects/subregion_XDR_no_subsampling')

subregion_XDR_no_subsampling_summary <- summary(subregion_XDR_no_subsampling)

# test fit 
testDispersion(subregion_XDR_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = subregion_XDR_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, subregion_model_full_df$subregion, quantreg = T)

# extract model coefficients

subregion_XDR_no_subsampling_model_coeff <- subregion_XDR_no_subsampling_summary$coefficients$cond
write.csv(subregion_XDR_no_subsampling_model_coeff, file = 'tables/subregion_XDR_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

subregion_XDR_no_subsampling_pred_df <- ggpredict(subregion_XDR_no_subsampling, terms = c('subregion'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

subregion_XDR_no_subsampling_pred_plot <- ggplot(subregion_XDR_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Subregion", y = "Prob. of XDR")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(subregion_min, XDR_binom, subregion, "Subregion", "Prob. of XDR")
model_fit_and_test(subregion_500, XDR_binom, subregion, "Subregion", "Prob. of XDR")
model_fit_and_test(subregion_1000, XDR_binom, subregion, "Subregion", "Prob. of XDR")
model_fit_and_test(subregion_2000, XDR_binom, subregion, "Subregion", "Prob. of XDR")
model_fit_and_test(subregion_4000, XDR_binom, subregion, "Subregion", "Prob. of XDR")

# arrange plots

arranged_XDR_subsamples_subregion <- ggarrange(subregion_XDR_no_subsampling_pred_plot, subregion_min_XDR_binom, subregion_500_XDR_binom,
                                          subregion_1000_XDR_binom, subregion_2000_XDR_binom, subregion_4000_XDR_binom,
                                          labels = c('A', 'B', 'C', 'D', 'E', 'F'))

############ SAVE ALL PLOTS ############

ggsave(plot = arranged_MDR_subsamples_phylogroup, "plots/MDR_phylogroup_subsampled_pred_plots.jpg", width=30, height=20, units="cm")
ggsave(plot = arranged_MDR_subsamples_host, "plots/MDR_host_subsampled_pred_plots.jpg", width=30, height=30, units="cm")
ggsave(plot = arranged_MDR_subsamples_subregion, "plots/MDR_subregion_subsampled_pred_plots.jpg", width=30, height=30, units="cm")

ggsave(plot = arranged_XDR_subsamples_phylogroup, "plots/XDR_phylogroup_subsampled_pred_plots.jpg", width=30, height=20, units="cm")
ggsave(plot = arranged_XDR_subsamples_host, "plots/XDR_host_subsampled_pred_plots.jpg", width=30, height=30, units="cm")
ggsave(plot = arranged_XDR_subsamples_subregion, "plots/XDR_subregion_subsampled_pred_plots.jpg", width=30, height=30, units="cm")

