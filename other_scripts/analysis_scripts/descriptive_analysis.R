#!/usr/bin/env RScript

library(tidyverse)
library(data.table)
library(ggpubr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(MetBrewer)

# set workspace

root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

# load in tidied dataframes

counts <- fread("AMR_counts_and_metadata.csv", sep=',') 

MDR_XDR <- fread("MDR_XDR_binomial.csv", sep=',')

ceph <- fread("CTX-M_binomial.csv", sep=',')

carb <- fread("carb_binomial.csv", sep=',')

# load in tables with group sample sizes for full dataset

phylogroup_counts <- fread("tables/phylogroup_count.csv", sep=',')
host_counts <- fread("tables/hosts_count.csv", sep=',')
subregion_counts <- fread("tables/subregion_count.csv", sep=',')

# manually line break axis labels for hosts

host_labels <- c("Agricultural/\ndomestic",
               "Human", 
               "Wild\nanimal",
               "Wild\nbird") 

# get overall mean count of AMR genes across dataset

mean_overall_amr_count <- mean(counts$count_amr)

# get median ARG counts for hosts, phylogroups and subregions

median_ARG_host <- counts %>%
  group_by(host) %>%
  summarize(median_ARG = median(count_amr, na.rm = TRUE))

median_ARG_phylogroup <- counts %>%
  group_by(phylogroup) %>%
  summarize(median_ARG = median(count_amr, na.rm = TRUE))

median_ARG_subregion <- counts %>%
  group_by(subregion) %>%
  summarize(median_ARG = median(count_amr, na.rm = TRUE))


# plot sample sizes for hosts, phylogroups and subregions

phylogroup_counts$phylogroup <- factor(phylogroup_counts$phylogroup, levels=c("A","B1", "B2", "C", "D", "E", "F", "G", "U", "U/cryptic", 'cryptic'))

phylogroup_count_plot <- ggplot(phylogroup_counts, aes(x=phylogroup, y=count)) +
  geom_col() +
  theme_light() +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(x = "Phylogroup", y = "No. of genomes") 

host_count_plot <- ggplot(host_counts, aes(x=host, y=count, fill=host)) +
  geom_col() +
  theme_light() +
  theme(text = element_text(size=10),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none") +
  scale_x_discrete(labels=host_labels) +
  scale_fill_manual(values=met.brewer("Juarez"), name="Host") +
  labs(x = "Host", y = "No. of genomes") 

subregion_count_plot <- ggplot(subregion_counts, aes(x=subregion, y=count)) +
  geom_col() +
  theme_light() +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(x = "Subregion", y = "No. of genomes")

# sample size tables for hosts, phylogroups and subregions

rename_columns_phylogroup_counts <- phylogroup_counts %>%
  rename(group = phylogroup) %>%
  mutate(variable = "phylogroup") %>%
  dplyr::select(-total)

rename_columns_host_counts <- host_counts %>%
  rename(group = host) %>%
  mutate(variable = "host") %>%
  dplyr::select(-total)

rename_columns_subregion_counts <- subregion_counts %>%
  rename(group = subregion) %>%
  mutate(variable = "subregion") %>%
  dplyr::select(-total)

sample_sizes_all <- rename_columns_phylogroup_counts %>%
  full_join(., rename_columns_host_counts) %>%
  full_join(., rename_columns_subregion_counts)

write.csv(sample_sizes_all, file = "tables/host_subregion_phylogroup_sample_sizes.csv", row.names = FALSE)

# plots of raw AMR counts per phylogroup, host & subregion

counts$phylogroup <- factor(counts$phylogroup, levels=c("A","B1", "B2", "C", "D", "E", "F", "G", "U", "U/cryptic", 'cryptic'))

raw_amr_phylogroup <- ggplot(counts, aes(x=phylogroup, y=count_amr)) +
  theme_light() +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  geom_jitter(alpha=0.05) +
  geom_boxplot(alpha=0.5) +
  labs(x = "Phylogroup", y = "AMR gene count") 

counts$host <- as.character(counts$host)

raw_amr_host <- ggplot(counts, aes(x=host, y=count_amr, fill=host)) +
  theme_light() +
  theme(text = element_text(size=10),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none") +
  geom_jitter(aes(colour=host), alpha=0.02) +
  geom_boxplot(alpha=0.5) +
  labs(x = "Host", y = "AMR gene count") +
  scale_x_discrete(labels=host_labels) +
  scale_fill_manual(values=met.brewer("Juarez"), name="Host") +
  scale_colour_manual(values=met.brewer("Juarez"), name="Host") 

raw_amr_subregion <- ggplot(counts, aes(x=subregion, y=count_amr)) +
  theme_light() +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  geom_jitter(alpha=0.05) +
  geom_boxplot(alpha=0.5) +
  labs(x = "Subregion", y = "AMR gene count") 

arrange_raw_plots <- ggarrange(raw_amr_phylogroup, raw_amr_host,
                               raw_amr_subregion,
                               labels = c("A", "B", "C"),
                               ncol=1,
                               nrow=3,
                               heights = c(1,1,2),
                               font.label = list(size = 25))

ggsave(plot=arrange_raw_plots, "plots/raw_counts_amr.tiff", width=10, height=20, dpi=300)
ggsave(plot=arrange_raw_plots, "plots/raw_counts_amr.jpg", width=10, height=20, dpi=300)

# plots of phylogroup distribution per host and subregion

# phylogroups per geographic region - proportions (plotting counts doesn't add
# much value)

total_subregions <- counts %>%
  dplyr::select(refseq_id, phylogroup, subregion) %>%
  unique() %>%
  count(subregion, name="total")

phylo_subregions <- counts %>%
  dplyr::select(refseq_id, phylogroup, subregion) %>%
  unique() %>%
  count(phylogroup, subregion) %>%
  full_join(., total_subregions) %>%
  mutate(prop = n/total)

phylo_subregion_plot_prop <- ggplot(phylo_subregions, aes(x=phylogroup, y=prop)) +
  geom_col() +
  facet_wrap(~subregion) +
  theme_light() +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(x = "Phylogroup", y = "Prop. of genomes") 

ggsave(plot=phylo_subregion_plot_prop, "plots/phylogroups_by_subregion.tiff", width=15, height=15, dpi=300)
ggsave(plot=phylo_subregion_plot_prop, "plots/phylogroups_by_subregion.jpg", width=15, height=15, dpi=300)

# phylogroups per host

total_hosts <- counts %>%
  dplyr::select(refseq_id, phylogroup, host) %>%
  unique() %>%
  count(host, name="total")

phylo_hosts <- counts %>%
  dplyr::select(refseq_id, phylogroup, host) %>%
  unique() %>%
  count(phylogroup, host) %>%
  full_join(., total_hosts) %>%
  mutate(prop = n/total)

phylo_host_plot <- ggplot(phylo_hosts, aes(x=phylogroup, y=prop)) +
  geom_col() +
  facet_wrap(~host) +
  theme_light() +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(x = "Phylogroup", y = "Prop. of genomes") 

ggsave(plot=phylo_host_plot, "plots/phylogroups_by_host.tiff", width=15, height=15, dpi=300)
ggsave(plot=phylo_host_plot, "plots/phylogroups_by_host.jpg", width=15, height=15, dpi=300)

# sample size tables for binomial categories: MDR, XDR, CTX-M & carbapenemases

sample_sizes_MDR_phylogroup <- MDR_XDR %>%
  dplyr::filter(MDR_binom == 1) %>%
  count(phylogroup, name = "count") %>%
  mutate(category = "MDR")

sample_sizes_XDR_phylogroup <- MDR_XDR %>%
  dplyr::filter(XDR_binom == 1) %>%
  count(phylogroup, name = "count") %>%
  mutate(category = "XDR")

sample_sizes_MDR_host <- MDR_XDR %>%
  dplyr::filter(MDR_binom == 1) %>%
  count(host, name = "count") %>%
  mutate(category = "MDR")

sample_sizes_XDR_host <- MDR_XDR %>%
  dplyr::filter(XDR_binom == 1) %>%
  count(host, name = "count") %>%
  mutate(category = "XDR")

sample_sizes_MDR_subregion <- MDR_XDR %>%
  dplyr::filter(MDR_binom == 1) %>%
  count(subregion, name = "count") %>%
  mutate(category = "MDR")

sample_sizes_XDR_subregion <- MDR_XDR %>%
  dplyr::filter(XDR_binom == 1) %>%
  count(subregion, name = "count") %>%
  mutate(category = "XDR")

write.csv(sample_sizes_MDR_phylogroup, file = "tables/sample_sizes_MDR_phylogroup.csv", row.names = FALSE)
write.csv(sample_sizes_XDR_phylogroup, file = "tables/sample_sizes_XDR_phylogroup.csv", row.names = FALSE)
write.csv(sample_sizes_MDR_host, file = "tables/sample_sizes_MDR_host.csv", row.names = FALSE)
write.csv(sample_sizes_XDR_host, file = "tables/sample_sizes_XDR_host.csv", row.names = FALSE)
write.csv(sample_sizes_MDR_subregion, file = "tables/sample_sizes_MDR_subregion.csv", row.names = FALSE)
write.csv(sample_sizes_XDR_subregion, file = "tables/sample_sizes_XDR_subregion.csv", row.names = FALSE)

sample_sizes_ceph_phylogroup <- ceph %>%
  dplyr::filter(ceph == 1) %>%
  count(phylogroup, name = "count") %>%
  mutate(category = "CTX-M")

sample_sizes_ceph_subregion <- ceph %>%
  dplyr::filter(ceph == 1) %>%
  count(subregion, name = "count") %>%
  mutate(category = "CTX-M")

sample_sizes_ceph_host <- ceph %>%
  dplyr::filter(ceph == 1) %>%
  count(host, name = "count") %>%
  mutate(category = "CTX-M")

write.csv(sample_sizes_ceph_phylogroup, file = "tables/sample_sizes_ctx-m_phylogroup.csv", row.names = FALSE)
write.csv(sample_sizes_ceph_subregion, file = "tables/sample_sizes_ctx-m_subregion.csv", row.names = FALSE)
write.csv(sample_sizes_ceph_host, file = "tables/sample_sizes_ctx-m_host.csv", row.names = FALSE)

sample_sizes_carb_phylogroup <- carb %>%
  dplyr::filter(carb == 1) %>%
  count(phylogroup, name = "count") %>%
  mutate(category = "carbapenemase")

sample_sizes_carb_subregion <- carb %>%
  dplyr::filter(carb == 1) %>%
  count(subregion, name = "count") %>%
  mutate(category = "carbapenemase")

sample_sizes_carb_host <- carb %>%
  dplyr::filter(carb == 1) %>%
  count(host, name = "count") %>%
  mutate(category = "carbapenemase")

write.csv(sample_sizes_carb_phylogroup, file = "tables/sample_sizes_carb_phylogroup.csv", row.names = FALSE)
write.csv(sample_sizes_carb_subregion, file = "tables/sample_sizes_carb_subregion.csv", row.names = FALSE)
write.csv(sample_sizes_carb_host, file = "tables/sample_sizes_carb_host.csv", row.names = FALSE)

##### make merged supplementary tables with full sample sizes ####

phylogroup_sample_sizes_dfs <- list(sample_sizes_MDR_phylogroup, sample_sizes_XDR_phylogroup,
                                sample_sizes_ceph_phylogroup, sample_sizes_carb_phylogroup)

phylogroup_sample_sizes <- phylogroup_sample_sizes_dfs %>%
  reduce(full_join)

write.csv(phylogroup_sample_sizes, file = "tables/all_sample_sizes_phylogroups.csv", row.names = FALSE)

host_sample_sizes_dfs <- list(sample_sizes_MDR_host, sample_sizes_XDR_host,
                                    sample_sizes_ceph_host, sample_sizes_carb_host)

host_sample_sizes <- host_sample_sizes_dfs %>%
  reduce(full_join)

write.csv(host_sample_sizes, file = "tables/all_sample_sizes_hosts.csv", row.names = FALSE)

subregion_sample_sizes_dfs <- list(sample_sizes_MDR_subregion, sample_sizes_XDR_subregion,
                              sample_sizes_ceph_subregion, sample_sizes_carb_subregion)

subregion_sample_sizes <- subregion_sample_sizes_dfs %>%
  reduce(full_join)

write.csv(subregion_sample_sizes, file = "tables/all_sample_sizes_subregions.csv", row.names = FALSE)

###### map showing location of subregions ######

world <- ne_countries(scale = "medium", returnclass = "sf")

subregion_palette <- c(rgb(160,71,102, maxColorValue=255),
             rgb(111,190,71, maxColorValue=255),
             rgb(175,92,211, maxColorValue=255),
             rgb(94,136,39, maxColorValue=255),
             rgb(99,108,217, maxColorValue=255),
             rgb(199,170,49, maxColorValue=255),
             rgb(206,79,173, maxColorValue=255),
             rgb(92,190,128, maxColorValue=255),
             rgb(217,64,126, maxColorValue=255),
             rgb(78,199,194, maxColorValue=255),
             rgb(206,60,66, maxColorValue=255),
             rgb(51,139,112, maxColorValue=255),
             rgb(222,107,48, maxColorValue=255),
             rgb(97,160,217, maxColorValue=255),
             rgb(213,149,85, maxColorValue=255),
             rgb(90,109,177, maxColorValue=255),
             rgb(172,176,98, maxColorValue=255),
             rgb(138,78,155, maxColorValue=255),
             rgb(71,126,68, maxColorValue=255),
             rgb(191,144,218, maxColorValue=255),
             rgb(129,110,43, maxColorValue=255),
             rgb(221,130,177, maxColorValue=255),
             rgb(162,82,47, maxColorValue=255),
             rgb(223,125,120, maxColorValue=255))

subregions <- ggplot(data = world) +
  geom_sf(aes(fill = subregion)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=10), 
        legend.text=element_text(size=10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = "right",
        plot.margin = unit(c(1, -5, 1, -5), "cm")) +
  scale_fill_manual(values=subregion_palette, name="Subregion") 
  

ggsave(plot=subregions, "plots/subregion_locations.jpg", dpi=200, width=35, height=10, units = "cm")

## arrange main descriptive figures into panel plot ##


arrange_sample_size_plots <- ggarrange(phylogroup_count_plot,
                                       subregion_count_plot,
                                       host_count_plot,
                                       labels = c("A", "B", "C"),
                                       ncol=1,
                                       nrow=3,
                                       heights = c(1,1.5,1),
                                       font.label = list(size = 15))



arrange_count_plots <- ggarrange(raw_amr_phylogroup,
                                 raw_amr_subregion,
                                 raw_amr_host,
                                 labels = c("D", "E", "F"),
                                 ncol=1,
                                 nrow=3,
                                 heights = c(1,1.5,1),
                                 font.label = list(size = 15))


arrange_main_descriptive_plots <- ggarrange(arrange_sample_size_plots,
                                            arrange_count_plots,
                                            phylo_host_plot,
                                           labels = c("", "", "G"),
                                           ncol=2,
                                           nrow=2,
                                           heights = c(3,1),
                                           font.label = list(size = 15))

ggsave(plot=arrange_main_descriptive_plots, "plots/main_descriptive_plots.jpg", dpi=300, width=20, height=30, units = "cm")
ggsave(plot=arrange_main_descriptive_plots, "plots/main_descriptive_plots.tiff", dpi=300, width=20, height=30, units = "cm")






