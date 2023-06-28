#!/usr/bin/env RScript

library(tidyverse)
library(data.table)

# set workspace

root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

################### model AIC and coefficient tables ###########################

##### read in data ######

# AIC tables #

MDR_dredge <- fread("tables/AIC_maximal_MDR_model.csv", sep=",") %>%
  mutate(model = "MDR")

XDR_dredge <- fread("tables/AIC_maximal_XDR_model.csv", sep=",") %>%
  mutate(model = "XDR")

ceph_dredge <- fread("tables/AIC_maximal_ceph_model.csv", sep=",") %>%
  mutate(model = "CTX-M")

carb_dredge <- fread("tables/AIC_maximal_carb_model.csv", sep=",") %>%
  mutate(model = "carbapenemase")

# maximal model coefficients #

maximal_MDR_model_coeff <- fread("tables/estimates_maximal_MDR_model.csv", sep=",") %>%
  mutate(model = "MDR") %>%
  rename(predictor = V1) # convert rownames with predictors into column

maximal_XDR_model_coeff <- fread("tables/estimates_maximal_XDR_model.csv", sep=",") %>%
  mutate(model = "XDR") %>%
  rename(predictor = V1) # convert rownames with predictors into column

maximal_ceph_model_coeff <- fread("tables/estimates_maximal_ceph_model.csv", sep=",") %>%
  mutate(model = "CTX-M") %>%
  rename(predictor = V1) # convert rownames with predictors into column

maximal_carb_model_coeff <- fread("tables/estimates_maximal_carb_model.csv", sep=",")%>%
  mutate(model = "carbapenemase") %>%
  rename(predictor = V1) # convert rownames with predictors into column

# MDR model coefficients #

# host

host_full_MDR_model_coeff <- fread("tables/host_MDR_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "MDR full dataset") %>%
  rename(predictor = V1)

host_min_MDR_model_coeff <- fread("subsampled_model_estimates/host_min_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR min") %>%
  rename(predictor = V1)
  
host_500_MDR_model_coeff <- fread("subsampled_model_estimates/host_500_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 500") %>%
  rename(predictor = V1)

host_1000_MDR_model_coeff <- fread("subsampled_model_estimates/host_1000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 1000") %>%
  rename(predictor = V1)

host_2000_MDR_model_coeff <- fread("subsampled_model_estimates/host_2000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 2000") %>%
  rename(predictor = V1)

host_4000_MDR_model_coeff <- fread("subsampled_model_estimates/host_4000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 4000") %>%
  rename(predictor = V1)

# phylogroup

phylogroup_full_MDR_model_coeff <- fread("tables/phylogroup_MDR_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "MDR full dataset") %>%
  rename(predictor = V1)

phylogroup_min_MDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_min_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR min") %>%
  rename(predictor = V1)

phylogroup_500_MDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_500_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 500") %>%
  rename(predictor = V1)

phylogroup_1000_MDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_1000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 1000") %>%
  rename(predictor = V1)

phylogroup_2000_MDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_2000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 2000") %>%
  rename(predictor = V1)

phylogroup_4000_MDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_4000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 4000") %>%
  rename(predictor = V1)

# subregion

subregion_full_MDR_model_coeff <- fread("tables/subregion_MDR_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "MDR full dataset") %>%
  rename(predictor = V1)

subregion_min_MDR_model_coeff <- fread("subsampled_model_estimates/subregion_min_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR min") %>%
  rename(predictor = V1)

subregion_500_MDR_model_coeff <- fread("subsampled_model_estimates/subregion_500_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 500") %>%
  rename(predictor = V1)

subregion_1000_MDR_model_coeff <- fread("subsampled_model_estimates/subregion_1000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 1000") %>%
  rename(predictor = V1)

subregion_2000_MDR_model_coeff <- fread("subsampled_model_estimates/subregion_2000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 2000") %>%
  rename(predictor = V1)

subregion_4000_MDR_model_coeff <- fread("subsampled_model_estimates/subregion_4000_MDR_binom.csv", sep=",") %>%
  mutate(model = "MDR 4000") %>%
  rename(predictor = V1)

# XDR model coefficients #

# host

host_full_XDR_model_coeff <- fread("tables/host_XDR_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "XDR full dataset") %>%
  rename(predictor = V1)

host_min_XDR_model_coeff <- fread("subsampled_model_estimates/host_min_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR min") %>%
  rename(predictor = V1)

host_500_XDR_model_coeff <- fread("subsampled_model_estimates/host_500_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 500") %>%
  rename(predictor = V1)

host_1000_XDR_model_coeff <- fread("subsampled_model_estimates/host_1000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 1000") %>%
  rename(predictor = V1)

host_2000_XDR_model_coeff <- fread("subsampled_model_estimates/host_2000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 2000") %>%
  rename(predictor = V1)

host_4000_XDR_model_coeff <- fread("subsampled_model_estimates/host_4000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 4000") %>%
  rename(predictor = V1)

# phylogroup

phylogroup_full_XDR_model_coeff <- fread("tables/phylogroup_XDR_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "XDR full dataset") %>%
  rename(predictor = V1)

phylogroup_min_XDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_min_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR min") %>%
  rename(predictor = V1)

phylogroup_500_XDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_500_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 500") %>%
  rename(predictor = V1)

phylogroup_1000_XDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_1000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 1000") %>%
  rename(predictor = V1)

phylogroup_2000_XDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_2000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 2000") %>%
  rename(predictor = V1)

phylogroup_4000_XDR_model_coeff <- fread("subsampled_model_estimates/phylogroup_4000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 4000") %>%
  rename(predictor = V1)

# subregion

subregion_full_XDR_model_coeff <- fread("tables/subregion_XDR_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "XDR full dataset") %>%
  rename(predictor = V1)

subregion_min_XDR_model_coeff <- fread("subsampled_model_estimates/subregion_min_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR min") %>%
  rename(predictor = V1)

subregion_500_XDR_model_coeff <- fread("subsampled_model_estimates/subregion_500_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 500") %>%
  rename(predictor = V1)

subregion_1000_XDR_model_coeff <- fread("subsampled_model_estimates/subregion_1000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 1000") %>%
  rename(predictor = V1)

subregion_2000_XDR_model_coeff <- fread("subsampled_model_estimates/subregion_2000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 2000") %>%
  rename(predictor = V1)

subregion_4000_XDR_model_coeff <- fread("subsampled_model_estimates/subregion_4000_XDR_binom.csv", sep=",") %>%
  mutate(model = "XDR 4000") %>%
  rename(predictor = V1)

# ceph model coefficients #

# host

host_full_ceph_model_coeff <- fread("tables/host_ceph_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "CTX-M full dataset") %>%
  rename(predictor = V1)

host_min_ceph_model_coeff <- fread("subsampled_model_estimates/host_min_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M min") %>%
  rename(predictor = V1)

host_500_ceph_model_coeff <- fread("subsampled_model_estimates/host_500_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 500") %>%
  rename(predictor = V1)

host_1000_ceph_model_coeff <- fread("subsampled_model_estimates/host_1000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 1000") %>%
  rename(predictor = V1)

host_2000_ceph_model_coeff <- fread("subsampled_model_estimates/host_2000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 2000") %>%
  rename(predictor = V1)

host_4000_ceph_model_coeff <- fread("subsampled_model_estimates/host_4000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 4000") %>%
  rename(predictor = V1)

# phylogroup

phylogroup_full_ceph_model_coeff <- fread("tables/phylogroup_ceph_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "CTX-M full dataset") %>%
  rename(predictor = V1)

phylogroup_min_ceph_model_coeff <- fread("subsampled_model_estimates/phylogroup_min_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M min") %>%
  rename(predictor = V1)

phylogroup_500_ceph_model_coeff <- fread("subsampled_model_estimates/phylogroup_500_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 500") %>%
  rename(predictor = V1)

phylogroup_1000_ceph_model_coeff <- fread("subsampled_model_estimates/phylogroup_1000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 1000") %>%
  rename(predictor = V1)

phylogroup_2000_ceph_model_coeff <- fread("subsampled_model_estimates/phylogroup_2000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 2000") %>%
  rename(predictor = V1)

phylogroup_4000_ceph_model_coeff <- fread("subsampled_model_estimates/phylogroup_4000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 4000") %>%
  rename(predictor = V1)

# subregion

subregion_full_ceph_model_coeff <- fread("tables/subregion_ceph_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "CTX-M full dataset") %>%
  rename(predictor = V1)

subregion_min_ceph_model_coeff <- fread("subsampled_model_estimates/subregion_min_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M min") %>%
  rename(predictor = V1)

subregion_500_ceph_model_coeff <- fread("subsampled_model_estimates/subregion_500_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 500") %>%
  rename(predictor = V1)

subregion_1000_ceph_model_coeff <- fread("subsampled_model_estimates/subregion_1000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 1000") %>%
  rename(predictor = V1)

subregion_2000_ceph_model_coeff <- fread("subsampled_model_estimates/subregion_2000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 2000") %>%
  rename(predictor = V1)

subregion_4000_ceph_model_coeff <- fread("subsampled_model_estimates/subregion_4000_ceph.csv", sep=",") %>%
  mutate(model = "CTX-M 4000") %>%
  rename(predictor = V1)

# carb model coefficients #

# host

host_full_carb_model_coeff <- fread("tables/host_carb_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "carbapenemase full dataset") %>%
  rename(predictor = V1)

host_min_carb_model_coeff <- fread("subsampled_model_estimates/host_min_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase min") %>%
  rename(predictor = V1)

host_500_carb_model_coeff <- fread("subsampled_model_estimates/host_500_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase 500") %>%
  rename(predictor = V1)

host_1000_carb_model_coeff <- fread("subsampled_model_estimates/host_1000_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase 1000") %>%
  rename(predictor = V1)

host_2000_carb_model_coeff <- fread("subsampled_model_estimates/host_2000_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase 2000") %>%
  rename(predictor = V1)

host_4000_carb_model_coeff <- fread("subsampled_model_estimates/host_4000_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase 4000") %>%
  rename(predictor = V1)

# phylogroup

phylogroup_full_carb_model_coeff <- fread("tables/phylogroup_carb_no_subsampling_model_coeff.csv", sep=",") %>%
  mutate(model = "carbapenemase full dataset") %>%
  rename(predictor = V1)

phylogroup_min_carb_model_coeff <- fread("subsampled_model_estimates/phylogroup_min_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase min") %>%
  rename(predictor = V1)

phylogroup_500_carb_model_coeff <- fread("subsampled_model_estimates/phylogroup_500_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase 500") %>%
  rename(predictor = V1)

phylogroup_1000_carb_model_coeff <- fread("subsampled_model_estimates/phylogroup_1000_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase 1000") %>%
  rename(predictor = V1)

phylogroup_2000_carb_model_coeff <- fread("subsampled_model_estimates/phylogroup_2000_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase 2000") %>%
  rename(predictor = V1)

phylogroup_4000_carb_model_coeff <- fread("subsampled_model_estimates/phylogroup_4000_carb.csv", sep=",") %>%
  mutate(model = "carbapenemase 4000") %>%
  rename(predictor = V1)

########## make tidy tables ##########

### maximal models ###

AIC_tables <- MDR_dredge %>%
  full_join(., XDR_dredge) %>%
  full_join(., ceph_dredge) %>%
  full_join(., carb_dredge)

maximal_coefficient_tables <- maximal_MDR_model_coeff %>%
  full_join(., maximal_XDR_model_coeff) %>%
  full_join(., maximal_ceph_model_coeff) %>%
  full_join(., maximal_carb_model_coeff)

### host models ###

host_coeff_dfs <- list(host_full_MDR_model_coeff, host_min_MDR_model_coeff, host_500_MDR_model_coeff, 
                  host_1000_MDR_model_coeff, host_2000_MDR_model_coeff, host_4000_MDR_model_coeff,
                  host_full_XDR_model_coeff, host_min_XDR_model_coeff, host_500_XDR_model_coeff, 
                  host_1000_XDR_model_coeff, host_2000_XDR_model_coeff, host_4000_XDR_model_coeff,
                  host_full_ceph_model_coeff, host_min_ceph_model_coeff, host_500_ceph_model_coeff, 
                  host_1000_ceph_model_coeff, host_2000_ceph_model_coeff, host_4000_ceph_model_coeff,
                  host_full_carb_model_coeff, host_min_carb_model_coeff, host_500_carb_model_coeff, 
                  host_1000_carb_model_coeff, host_2000_carb_model_coeff, host_4000_carb_model_coeff)

host_coeffs <- host_coeff_dfs %>% 
  reduce(full_join)

### phylogroup models ###

phylogroup_coeff_dfs <- list(phylogroup_full_MDR_model_coeff, phylogroup_min_MDR_model_coeff, phylogroup_500_MDR_model_coeff, 
                       phylogroup_1000_MDR_model_coeff, phylogroup_2000_MDR_model_coeff, phylogroup_4000_MDR_model_coeff,
                       phylogroup_full_XDR_model_coeff, phylogroup_min_XDR_model_coeff, phylogroup_500_XDR_model_coeff, 
                       phylogroup_1000_XDR_model_coeff, phylogroup_2000_XDR_model_coeff, phylogroup_4000_XDR_model_coeff,
                       phylogroup_full_ceph_model_coeff, phylogroup_min_ceph_model_coeff, phylogroup_500_ceph_model_coeff, 
                       phylogroup_1000_ceph_model_coeff, phylogroup_2000_ceph_model_coeff, phylogroup_4000_ceph_model_coeff,
                       phylogroup_full_carb_model_coeff, phylogroup_min_carb_model_coeff, phylogroup_500_carb_model_coeff, 
                       phylogroup_1000_carb_model_coeff, phylogroup_2000_carb_model_coeff, phylogroup_4000_carb_model_coeff)

phylogroup_coeffs <- phylogroup_coeff_dfs %>% 
  reduce(full_join)

### subregion models ###

subregion_coeff_dfs <- list(subregion_full_MDR_model_coeff, subregion_min_MDR_model_coeff, subregion_500_MDR_model_coeff, 
                             subregion_1000_MDR_model_coeff, subregion_2000_MDR_model_coeff, subregion_4000_MDR_model_coeff,
                             subregion_full_XDR_model_coeff, subregion_min_XDR_model_coeff, subregion_500_XDR_model_coeff, 
                             subregion_1000_XDR_model_coeff, subregion_2000_XDR_model_coeff, subregion_4000_XDR_model_coeff,
                             subregion_full_ceph_model_coeff, subregion_min_ceph_model_coeff, subregion_500_ceph_model_coeff, 
                             subregion_1000_ceph_model_coeff, subregion_2000_ceph_model_coeff, subregion_4000_ceph_model_coeff)

subregion_coeffs <- subregion_coeff_dfs %>% 
  reduce(full_join)

# write out as tables

write.csv(AIC_tables, file = "tables/all_maximal_model_AICs.csv")
write.csv(maximal_coefficient_tables, file = "tables/all_maximal_model_coeffs.csv")

write.csv(host_coeffs, file = "tables/host_model_coeffs.csv")
write.csv(phylogroup_coeffs, file = "tables/phylogroup_model_coeffs.csv")
write.csv(subregion_coeffs, file = "tables/subregion_model_coeffs.csv")

