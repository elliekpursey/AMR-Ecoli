#!/usr/bin/env RScript

library(tidyverse)
library(data.table)
library(glmmTMB)
library(DHARMa)
library(MASS)
library(MuMIn)

# set workspace

root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

# load in tidied dataframes and convert predictor/identifier columns to factors

counts <- fread("AMR_counts_and_metadata.csv", sep=',') %>%
  mutate(refseq_id = as.factor(refseq_id))

MDR_XDR <- fread("MDR_XDR_binomial.csv", sep=',') %>%
  mutate(refseq_id = as.factor(refseq_id))

ceph <- fread("CTX-M_binomial.csv", sep=',') %>%
  mutate(refseq_id = as.factor(refseq_id))

carb <- fread("carb_binomial.csv", sep=',') %>%
  mutate(refseq_id = as.factor(refseq_id))

# function to remove NAs from factor columns and replace with '-' for writing
# out AIC tables neatly

replace_factor_na <- function(x){
  x <- as.character(x)
  x <- if_else(is.na(x), "-", x)
  x <- as.factor(x)
}


# remove groups with sample sizes <=100 for all categories in subregion and phylogroup
# (we know host groups are all above 100)

# get sample sizes

phylogroup_sample_sizes <- counts %>%
  dplyr::select(phylogroup) %>%
  count(phylogroup)

subregion_sample_sizes <- counts %>%
  dplyr::select(subregion) %>%
  count(subregion)

# identify groups with small sample sizes - (below 100)

phylogroups_to_keep <- phylogroup_sample_sizes %>%
  dplyr::filter(n >= 100) 

subregions_to_keep <- subregion_sample_sizes %>%
  dplyr::filter(n >= 100)

# filter modelling datasets 

MDR_XDR_filtered <- MDR_XDR %>%
  dplyr::filter(phylogroup %in% phylogroups_to_keep$phylogroup) %>%
  dplyr::filter(subregion %in% subregions_to_keep$subregion) %>%
  mutate(host = as.factor(host)) %>%
  mutate(phylogroup = as.factor(phylogroup)) %>%
  mutate(subregion = as.factor(subregion))

ceph_filtered <- ceph %>%
  dplyr::filter(phylogroup %in% phylogroups_to_keep$phylogroup) %>%
  dplyr::filter(subregion %in% subregions_to_keep$subregion) %>%
  mutate(host = as.factor(host)) %>%
  mutate(phylogroup = as.factor(phylogroup)) %>%
  mutate(subregion = as.factor(subregion))

carb_filtered <- carb %>%
  dplyr::filter(phylogroup %in% phylogroups_to_keep$phylogroup) %>%
  dplyr::filter(subregion %in% subregions_to_keep$subregion)%>%
  mutate(host = as.factor(host)) %>%
  mutate(phylogroup = as.factor(phylogroup)) %>%
  mutate(subregion = as.factor(subregion))

# histogram of response variable 

ggplot(counts, aes(x=count_amr)) +
  geom_histogram()

# maximal model for MDR

maximal_MDR_model <- glmmTMB(MDR_binom ~ host + phylogroup + subregion,
                               data=MDR_XDR_filtered, na.action = "na.fail", family = "binomial")

maximal_MDR_model_summary <- summary(maximal_MDR_model)

testDispersion(maximal_MDR_model)
simulationOutput <- simulateResiduals(fittedModel = maximal_MDR_model)
plot(simulationOutput)

plotResiduals(simulationOutput, MDR_XDR_filtered$host, quantreg = T)
plotResiduals(simulationOutput, MDR_XDR_filtered$subregion, quantreg = T)
plotResiduals(simulationOutput, MDR_XDR_filtered$phylogroup, quantreg = T)

MDR_dredge <- as.data.frame(dredge(maximal_MDR_model)) %>%
  mutate_if(is.factor, replace_factor_na) 

write.csv(MDR_dredge, file = "tables/AIC_maximal_MDR_model.csv", row.names=FALSE)

maximal_MDR_model_coeff <- as.data.frame(maximal_MDR_model_summary$coefficients$cond) 

write.csv(maximal_MDR_model_coeff, file = "tables/estimates_maximal_MDR_model.csv")

# maximal model for XDR

maximal_XDR_model <- glmmTMB(XDR_binom ~ host + phylogroup + subregion,
                             data=MDR_XDR_filtered, na.action = "na.fail", family = "binomial")

maximal_XDR_model_summary <- summary(maximal_XDR_model)

testDispersion(maximal_XDR_model)
simulationOutput <- simulateResiduals(fittedModel = maximal_XDR_model)
plot(simulationOutput)

plotResiduals(simulationOutput, MDR_XDR_filtered$host, quantreg = T)
plotResiduals(simulationOutput, MDR_XDR_filtered$subregion, quantreg = T)
plotResiduals(simulationOutput, MDR_XDR_filtered$phylogroup, quantreg = T)

XDR_dredge <- as.data.frame(dredge(maximal_XDR_model)) %>%
  mutate_if(is.factor, replace_factor_na) 

write.csv(XDR_dredge, file = "tables/AIC_maximal_XDR_model.csv", row.names=FALSE)

maximal_XDR_model_coeff <- as.data.frame(maximal_XDR_model_summary$coefficients$cond) 

write.csv(maximal_XDR_model_coeff, file = "tables/estimates_maximal_XDR_model.csv")

# maximal model for CTX-M

maximal_ceph_model <- glmmTMB(ceph ~ host + phylogroup + subregion,
                             data=ceph_filtered, na.action = "na.fail", family = "binomial")

maximal_ceph_model_summary <- summary(maximal_ceph_model)

testDispersion(maximal_ceph_model)
simulationOutput <- simulateResiduals(fittedModel = maximal_ceph_model)
plot(simulationOutput)

plotResiduals(simulationOutput, ceph_filtered$host, quantreg = T)
plotResiduals(simulationOutput, ceph_filtered$subregion, quantreg = T)
plotResiduals(simulationOutput, ceph_filtered$phylogroup, quantreg = T)

ceph_dredge <- as.data.frame(dredge(maximal_ceph_model)) %>%
  mutate_if(is.factor, replace_factor_na)

write.csv(ceph_dredge, file = "tables/AIC_maximal_ceph_model.csv", row.names=FALSE)

maximal_ceph_model_coeff <- as.data.frame(maximal_ceph_model_summary$coefficients$cond) 

write.csv(maximal_ceph_model_coeff, file = "tables/estimates_maximal_ceph_model.csv")

# maximal model for carbapenemases

carb_count_subregion <- carb %>%
  count(subregion, carb)

# this model does not converge with subregion included

maximal_carb_model <- glmmTMB(carb ~ host + phylogroup,
                              data=carb_filtered, na.action = "na.fail", family = "binomial")

maximal_carb_model_summary <- summary(maximal_carb_model)

testDispersion(maximal_carb_model)
simulationOutput <- simulateResiduals(fittedModel = maximal_carb_model)
plot(simulationOutput)

plotResiduals(simulationOutput, carb_filtered$host, quantreg = T)
plotResiduals(simulationOutput, carb_filtered$phylogroup, quantreg = T)

carb_dredge <- as.data.frame(dredge(maximal_carb_model)) %>%
  mutate_if(is.factor, replace_factor_na) 

write.csv(carb_dredge, file = "tables/AIC_maximal_carb_model.csv", row.names=FALSE)

maximal_carb_model_coeff <- as.data.frame(maximal_carb_model_summary$coefficients$cond) 

write.csv(maximal_carb_model_coeff, file = "tables/estimates_maximal_carb_model.csv")



