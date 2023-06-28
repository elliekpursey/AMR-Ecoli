library(tidyverse)
library(data.table)
library(MetBrewer)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# set path
root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

# read in data frames 

metadata <- fread("host_sort_EC.csv", sep=',') %>%
  dplyr::select(refseq_id, category, isolation_source, location, date)

phylogroups <-  fread("all_phylogroups.csv", sep=',', col.names = c('genome', 'phylogroup'))

# tidy metadata columns, match country names to world dataset
metadata_tidy <- metadata %>%
  separate(col=location, into=c("name_long"), sep=":") %>%
  mutate(category = str_replace(category, "_", " ")) %>%
  mutate(name_long = str_replace(name_long, "Viet Nam", "Vietnam")) %>%
  mutate(name_long = str_replace(name_long, "USA", "United States")) %>%
  mutate(name_long = str_replace(name_long, "Russia", "Russian Federation")) %>%
  mutate(name_long = str_replace(name_long, "Gambia", "The Gambia")) %>%
  mutate(name_long = str_replace(name_long, "Laos", "Lao PDR")) %>%
  mutate(name_long = str_replace(name_long, "South Korea", "Republic of Korea")) %>%
  mutate(name_long = str_replace(name_long, "Brunei", "Brunei Darussalam")) %>%
  mutate(name_long = str_replace(name_long, "Republic of the Congo", "Republic of Congo")) %>%
  mutate(name_long = str_replace(name_long, "Zaire", "Democratic Republic of the Congo")) %>%
  mutate(name_long = str_replace(name_long, "Democratic Republic of Congo", "Democratic Republic of the Congo")) %>%
  drop_na()

metadata_remove_missing <- metadata_tidy %>%
  dplyr::filter(name_long != "missing" &
                  name_long != "" & 
                  name_long != "not determined" &
                  name_long != "Na" & 
                  name_long != "not collected" &
                  name_long != "not available" &
                  name_long != "na" &
                  name_long != "Not collected" &
                  name_long != "Not applicable" &
                  name_long != "Not Collected" &
                  name_long != "Not Applicable" &
                  name_long != "not applicable" &
                  name_long != "unknown" & 
                  name_long != "Unknown" &
                  name_long != "Missing") %>%
  dplyr::filter(category != "missing" & category != "other") %>% # remove missing data for hosts
  rename(host = category) 

tidy_phylogroups <- phylogroups %>%
  separate(col="genome", into=c("GCF", "genome"), sep="_") %>%
  unite(refseq_id, GCF, genome, sep = "_", remove = FALSE) %>%
  dplyr::select(-GCF, -genome) %>%
  dplyr::filter(phylogroup != "EC_control_fail")  # remove EC_control_fail

# load in map data 

world <- ne_countries(scale = "medium", returnclass = "sf")

world <- as_tibble(world) %>%
  dplyr::select(name_long, subregion, continent, iso_a3) %>%
  rename(ISO_3DIGIT = iso_a3)

# join world data to genome data and phylogroups 

full_df <- metadata_remove_missing %>%
  left_join(., world) %>%
  left_join(., tidy_phylogroups) %>%
  drop_na() # removes 2 datapoints with countries that can't be assigned (Korea and Atlantic Ocean)

# count number of each host category, country, region and continent

hosts_count <- full_df %>%
  count(host, name="count") %>%
  mutate(total = sum(count)) 

country_count <- full_df %>%
  count(name_long, name="count") %>%
  mutate(total = sum(count)) %>%
  rename(country = name_long)

subregion_count <- full_df %>%
  count(subregion, name="count") %>%
  mutate(total = sum(count)) 

phylogroup_count <- full_df %>%
  count(phylogroup, name="count") %>%
  mutate(total = sum(count)) 

write_csv(hosts_count,"tables/hosts_count.csv")
write_csv(country_count,"tables/country_count.csv")
write_csv(subregion_count,"tables/subregion_count.csv")
write_csv(phylogroup_count,"tables/phylogroup_count.csv")

write_csv(full_df, "metadata_phylogroups.csv")
