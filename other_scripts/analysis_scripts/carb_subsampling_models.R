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

carb <- fread("carb_binomial.csv", sep=',')

# analysis of proportions in raw data

raw_proportions <- carb %>%
  group_by(phylogroup) %>%
  mutate(total = sum(carb)) %>%
  mutate(phylo_total = n()) %>%
  mutate(prop = total / n()) %>%
  ungroup() %>%
  dplyr::select(phylogroup, prop, total,phylo_total) %>%
  unique()

# read in function 

source('../other_scripts_local/functions/model_fit_and_test.R')

############################### PHYLOGROUPS ####################################

# read in subsampled/resampled data and merge with carb data to get relevant
# subsampled dataframes

phylogroup_min <- fread("subsampled_dataframes/subsampled_phylogroups_no_rep.csv", sep=",") %>%
  left_join(., carb)

phylogroup_500 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_500.csv", sep=",") %>%
  left_join(., carb)

phylogroup_1000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_1000.csv", sep=",") %>%
  left_join(., carb)

phylogroup_2000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_2000.csv", sep=",") %>%
  left_join(., carb)

phylogroup_4000 <- fread("subsampled_dataframes/phylogroup_resample_with_rep_4000.csv", sep=",") %>%
  left_join(., carb)

# model without subsampling

# get list of phylogroups removed for subsampling and use full dataframe with 
# these groups removed for consistency

list_of_phylogroups <- phylogroup_min %>%
  dplyr::select(phylogroup) %>%
  unique()

phylogroup_model_full_df <- carb %>%
  dplyr::filter(phylogroup %in% list_of_phylogroups$phylogroup) 

# run model

phylogroup_carb_no_subsampling <- glmmTMB(carb ~ phylogroup, data=phylogroup_model_full_df, na.action = "na.fail", family = "binomial")

saveRDS(phylogroup_carb_no_subsampling, 'model_objects/phylogroup_carb_no_subsampling')

phylogroup_carb_no_subsampling_summary <- summary(phylogroup_carb_no_subsampling)

# test fit 
testDispersion(phylogroup_carb_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = phylogroup_carb_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, phylogroup_model_full_df$phylogroup, quantreg = T)

# extract model coefficients

phylogroup_carb_no_subsampling_model_coeff <- phylogroup_carb_no_subsampling_summary$coefficients$cond
write.csv(phylogroup_carb_no_subsampling_model_coeff, file = 'tables/phylogroup_carb_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

phylogroup_carb_no_subsampling_pred_df <- ggpredict(phylogroup_carb_no_subsampling, terms = c('phylogroup'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

phylogroup_carb_no_subsampling_pred_plot <- ggplot(phylogroup_carb_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Prob. of carb")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(phylogroup_min, carb, phylogroup, "Phylogroup", "Prob. of carb")
model_fit_and_test(phylogroup_500, carb, phylogroup, "Phylogroup", "Prob. of carb")
model_fit_and_test(phylogroup_1000, carb, phylogroup, "Phylogroup", "Prob. of carb")
model_fit_and_test(phylogroup_2000, carb, phylogroup, "Phylogroup", "Prob. of carb")
model_fit_and_test(phylogroup_4000, carb, phylogroup, "Phylogroup", "Prob. of carb")

# arrange plots

arranged_carb_subsamples_phylogroup <- ggarrange(phylogroup_carb_no_subsampling_pred_plot, phylogroup_min_carb, phylogroup_500_carb,
                                                 phylogroup_1000_carb, phylogroup_2000_carb, phylogroup_4000_carb,
                                                 labels = c('A', 'B', 'C', 'D', 'E', 'F'))

############################### HOSTS ####################################

# read in subsampled/resampled data and merge with carb data to get relevant
# subsampled dataframes

host_min <- fread("subsampled_dataframes/subsampled_hosts_no_rep.csv", sep=",") %>%
  left_join(., carb)

host_500 <- fread("subsampled_dataframes/host_resample_with_rep_500.csv", sep=",") %>%
  left_join(., carb)

host_1000 <- fread("subsampled_dataframes/host_resample_with_rep_1000.csv", sep=",") %>%
  left_join(., carb)

host_2000 <- fread("subsampled_dataframes/host_resample_with_rep_2000.csv", sep=",") %>%
  left_join(., carb)

host_4000 <- fread("subsampled_dataframes/host_resample_with_rep_4000.csv", sep=",") %>%
  left_join(., carb)

# model without subsampling

# run model

host_carb_no_subsampling <- glmmTMB(carb ~ host, data=carb, na.action = "na.fail", family = "binomial")

saveRDS(host_carb_no_subsampling, 'model_objects/host_carb_no_subsampling')

host_carb_no_subsampling_summary <- summary(host_carb_no_subsampling)

# test fit 
testDispersion(host_carb_no_subsampling)
simulationOutput <- simulateResiduals(fittedModel = host_carb_no_subsampling)
plot(simulationOutput)

plotResiduals(simulationOutput, carb$host, quantreg = T)

# extract model coefficients

host_carb_no_subsampling_model_coeff <- host_carb_no_subsampling_summary$coefficients$cond
write.csv(host_carb_no_subsampling_model_coeff, file = 'tables/host_carb_no_subsampling_model_coeff.csv')

# create prediction dataframe and plot

host_carb_no_subsampling_pred_df <- ggpredict(host_carb_no_subsampling, terms = c('host'), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

host_carb_no_subsampling_pred_plot <- ggplot(host_carb_no_subsampling_pred_df, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(size=3, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
  labs(x = "Host", y = "Prob. of carb")

# use function to run models on all subsampled data, write out estimates, plots etc.

model_fit_and_test(host_min, carb, host, "Host", "Prob. of carb")
model_fit_and_test(host_500, carb, host, "Host", "Prob. of carb")
model_fit_and_test(host_1000, carb, host, "Host", "Prob. of carb")
model_fit_and_test(host_2000, carb, host, "Host", "Prob. of carb")
model_fit_and_test(host_4000, carb, host, "Host", "Prob. of carb")

# arrange plots

arranged_carb_subsamples_host <- ggarrange(host_carb_no_subsampling_pred_plot, host_min_carb, host_500_carb,
                                           host_1000_carb, host_2000_carb, host_4000_carb,
                                           labels = c('A', 'B', 'C', 'D', 'E', 'F'))


############ SAVE ALL PLOTS ############

ggsave(plot = arranged_carb_subsamples_phylogroup, "plots/carb_phylogroup_subsampled_pred_plots.jpg", width=30, height=20, units="cm")
ggsave(plot = arranged_carb_subsamples_host, "plots/carb_host_subsampled_pred_plots.jpg", width=30, height=20, units="cm")

ggsave(plot = arranged_carb_subsamples_phylogroup, "plots/carb_phylogroup_subsampled_pred_plots.jpg", width=30, height=20, units="cm")
ggsave(plot = arranged_carb_subsamples_host, "plots/carb_host_subsampled_pred_plots.jpg", width=30, height=20, units="cm")


