library(tidyverse)
library(data.table)
library(ggpubr)
library(ggeffects)
library(glmmTMB)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(MetBrewer)

# set workspace

root_directory <- ""

# set working dir

working_dir <- paste(root_directory, "results", sep = "/")
setwd(working_dir)

###############################################################################

# read in model objects

# MDR

min_MDR_phylogroup <- readRDS("model_objects/phylogroup_min_MDR_binom")
min_MDR_host <- readRDS("model_objects/host_min_MDR_binom")
min_MDR_subregion <- readRDS("model_objects/subregion_min_MDR_binom")
full_MDR_phylogroup <- readRDS("model_objects/phylogroup_MDR_no_subsampling")
full_MDR_host <- readRDS("model_objects/host_MDR_no_subsampling")
full_MDR_subregion <- readRDS("model_objects/subregion_MDR_no_subsampling")

# XDR

min_XDR_phylogroup <- readRDS("model_objects/phylogroup_min_XDR_binom")
min_XDR_host <- readRDS("model_objects/host_min_XDR_binom")
min_XDR_subregion <- readRDS("model_objects/subregion_min_XDR_binom")
full_XDR_phylogroup <- readRDS("model_objects/phylogroup_XDR_no_subsampling")
full_XDR_host <- readRDS("model_objects/host_XDR_no_subsampling")
full_XDR_subregion <- readRDS("model_objects/subregion_XDR_no_subsampling")

# CTX-M

min_ceph_phylogroup <- readRDS("model_objects/phylogroup_min_ceph")
min_ceph_host <- readRDS("model_objects/host_min_ceph")
min_ceph_subregion <- readRDS("model_objects/subregion_min_ceph")
full_ceph_phylogroup <- readRDS("model_objects/phylogroup_ceph_no_subsampling")
full_ceph_host <- readRDS("model_objects/host_ceph_no_subsampling")
full_ceph_subregion <- readRDS("model_objects/subregion_ceph_no_subsampling")

# carbapenem res.

min_carb_phylogroup <- readRDS("model_objects/phylogroup_min_carb")
min_carb_host <- readRDS("model_objects/host_min_carb")
full_carb_phylogroup <- readRDS("model_objects/phylogroup_carb_no_subsampling")
full_carb_host <- readRDS("model_objects/host_carb_no_subsampling")

# load world dataset

world <- ne_countries(scale = "medium", returnclass = "sf")

# pictures of hosts

host_pics <- c(
  'Agricultural/domestic' = "<img src = '../resources/cow.jpg' width = '60' height = '50' /><br>Agricultural/domestic",
  'Human' = "<img src = '../resources/people.jpg' width = '75' height = '50' /><br>Human",
  'Wild animal' = "<img src = '../resources/fox.jpg' width = '45' height = '50' /><br>Wild animal",
  'Wild bird' = "<img src = '../resources/crow.jpg' width = '45' height = '50' /><br>Wild bird")

############################# make figures ####################################

### MDR figure ###

## phylogroup

# subsampled to min group size

pred_min_MDR_phylogroup <- ggpredict(min_MDR_phylogroup, terms = c("phylogroup"), type='fixed', ci.lvl = 0.95)

plot_MDR_min_phylogroup <- ggplot(pred_min_MDR_phylogroup, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Predicted probability\n of MDR") 

# full dataset

pred_full_MDR_phylogroup <- ggpredict(full_MDR_phylogroup, terms = c("phylogroup"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

plot_MDR_full_phylogroup <- ggplot(pred_full_MDR_phylogroup, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Predicted probability\n of MDR") 

## host

# subsampled to min group size

pred_MDR_min_host <- ggpredict(min_MDR_host, terms = c("host"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x)) %>%
  mutate(x = str_replace(x, "agricultural/domestic", "Agricultural/domestic"))  %>%
  mutate(x = str_replace(x, "human", "Human")) %>%
  mutate(x = str_replace(x, "wild animal", "Wild animal")) %>%
  mutate(x = str_replace(x, "wild bird", "Wild bird"))

plot_MDR_min_host <- ggplot(pred_MDR_min_host, aes(x = x, y = predicted, colour=x))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Host species", y = "Predicted probability\n of MDR") +
  scale_colour_manual(values=met.brewer("Juarez")) +
  scale_x_discrete(name = NULL,
                   labels = host_pics) +
  theme(axis.text.x = ggtext::element_markdown(size = 12))

# full dataset

pred_MDR_full_host <- ggpredict(full_MDR_host, terms = c("host"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x)) %>%
  mutate(x = str_replace(x, "agricultural/domestic", "Agricultural/domestic"))  %>%
  mutate(x = str_replace(x, "human", "Human")) %>%
  mutate(x = str_replace(x, "wild animal", "Wild animal")) %>%
  mutate(x = str_replace(x, "wild bird", "Wild bird"))

plot_MDR_full_host <- ggplot(pred_MDR_full_host, aes(x = x, y = predicted, colour=x))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Host species", y = "Predicted probability\n of MDR") +
  scale_colour_manual(values=met.brewer("Juarez")) +
  scale_x_discrete(name = NULL,
                   labels = host_pics) +
  theme(axis.text.x = ggtext::element_markdown(size = 12))

## subregion

# subsampled to min group size

pred_MDR_min_subregion <- ggpredict(min_MDR_subregion, terms = c("subregion"), type='fixed', ci.lvl = 0.95) %>%
  rename(subregion = x) %>%
  dplyr::select(-group)

map_MDR_min_subregion <- world %>%
  full_join(pred_MDR_min_subregion)

plot_MDR_min_subregion <- ggplot(data = map_MDR_min_subregion) +
  geom_sf(aes(fill = predicted)) +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), name = "Predicted probability\n of MDR") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=15), 
        legend.text=element_text(size=10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = "bottom",
        plot.margin = unit(c(1, -5, 1, -5), "cm"))

# full dataset

pred_MDR_full_subregion <- ggpredict(full_MDR_subregion, terms = c("subregion"), type='fixed', ci.lvl = 0.95) %>%
  rename(subregion = x) %>%
  dplyr::select(-group)

map_MDR_full_subregion <- world %>%
  full_join(pred_MDR_full_subregion)

plot_MDR_full_subregion <- ggplot(data = map_MDR_full_subregion) +
  geom_sf(aes(fill = predicted)) +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), name = "Predicted probability\n of MDR") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=15), 
        legend.text=element_text(size=10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = "bottom",
        plot.margin = unit(c(1, -5, 1, -5), "cm"))

## arrange

# 1 row = 1 predictor, left column is full dataset, right column is min dataset

MDR_figure_arranged <- ggarrange(plot_MDR_min_phylogroup, plot_MDR_full_phylogroup,
                                 plot_MDR_min_host, plot_MDR_full_host,
                                 plot_MDR_min_subregion, plot_MDR_full_subregion,
                                 nrow=3, ncol=2,
                                 heights = c(1,1,1.15),
                                 labels = c('A', 'D', 'B', 'E', 'C', 'F'))

ggsave(plot = MDR_figure_arranged, "plots/MDR_final_fig.jpg", width=30, height=30, units="cm", dpi=300)
ggsave(plot = MDR_figure_arranged, "plots/MDR_final_fig.svg", width=30, height=30, units="cm", dpi=300)


### XDR figure ###

## phylogroup

# subsampled to min group size

pred_min_XDR_phylogroup <- ggpredict(min_XDR_phylogroup, terms = c("phylogroup"), type='fixed', ci.lvl = 0.95)

plot_XDR_min_phylogroup <- ggplot(pred_min_XDR_phylogroup, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  coord_cartesian(ylim = c(0.00, 0.15)) +
  labs(x = "Phylogroup", y = "Predicted probability\n of XDR") 

# full dataset

pred_full_XDR_phylogroup <- ggpredict(full_XDR_phylogroup, terms = c("phylogroup"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

plot_XDR_full_phylogroup <- ggplot(pred_full_XDR_phylogroup, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Predicted probability\n of XDR") 

## host

# subsampled to min group size

pred_XDR_min_host <- ggpredict(min_XDR_host, terms = c("host"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x)) %>%
  mutate(x = str_replace(x, "agricultural/domestic", "Agricultural/domestic"))  %>%
  mutate(x = str_replace(x, "human", "Human")) %>%
  mutate(x = str_replace(x, "wild animal", "Wild animal")) %>%
  mutate(x = str_replace(x, "wild bird", "Wild bird"))

plot_XDR_min_host <- ggplot(pred_XDR_min_host, aes(x = x, y = predicted, colour=x))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  coord_cartesian(ylim = c(0.00, 0.20)) +
  labs(x = "Host species", y = "Predicted probability\n of XDR") +
  scale_colour_manual(values=met.brewer("Juarez")) +
  scale_x_discrete(name = NULL,
                   labels = host_pics) +
  theme(axis.text.x = ggtext::element_markdown(size = 12))

# full dataset

pred_XDR_full_host <- ggpredict(full_XDR_host, terms = c("host"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x)) %>%
  mutate(x = str_replace(x, "agricultural/domestic", "Agricultural/domestic"))  %>%
  mutate(x = str_replace(x, "human", "Human")) %>%
  mutate(x = str_replace(x, "wild animal", "Wild animal")) %>%
  mutate(x = str_replace(x, "wild bird", "Wild bird"))

plot_XDR_full_host <- ggplot(pred_XDR_full_host, aes(x = x, y = predicted, colour=x))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Host species", y = "Predicted probability\n of XDR") +
  scale_colour_manual(values=met.brewer("Juarez")) +
  scale_x_discrete(name = NULL,
                   labels = host_pics) +
  theme(axis.text.x = ggtext::element_markdown(size = 12))

## subregion

# subsampled to min group size

pred_XDR_min_subregion <- ggpredict(min_XDR_subregion, terms = c("subregion"), type='fixed', ci.lvl = 0.95) %>%
  rename(subregion = x) %>%
  dplyr::select(-group) %>%
  dplyr::filter(conf.high < '1') # removes regions with error bar spanning full range, misleading as error bars not shown on map

map_XDR_min_subregion <- world %>%
  full_join(pred_XDR_min_subregion)

plot_XDR_min_subregion <- ggplot(data = map_XDR_min_subregion) +
  geom_sf(aes(fill = predicted)) +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), name = "Predicted probability\n of XDR") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=15), 
        legend.text=element_text(size=10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = "bottom",
        plot.margin = unit(c(1, -5, 1, -5), "cm"))

# full dataset

pred_XDR_full_subregion <- ggpredict(full_XDR_subregion, terms = c("subregion"), type='fixed', ci.lvl = 0.95) %>%
  rename(subregion = x) %>%
  dplyr::select(-group) %>%
  dplyr::filter(conf.high < '1') # removes regions with error bar spanning full range, misleading as error bars not shown on map

map_XDR_full_subregion <- world %>%
  full_join(pred_XDR_full_subregion)

plot_XDR_full_subregion <- ggplot(data = map_XDR_full_subregion) +
  geom_sf(aes(fill = predicted)) +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), name = "Predicted probability\n of XDR") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=15), 
        legend.text=element_text(size=10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = "bottom",
        plot.margin = unit(c(1, -5, 1, -5), "cm"))

## arrange

# 1 row = 1 predictor, left column is full dataset, right column is min dataset

XDR_figure_arranged <- ggarrange(plot_XDR_min_phylogroup, plot_XDR_full_phylogroup,
                                 plot_XDR_min_host, plot_XDR_full_host,
                                 plot_XDR_min_subregion, plot_XDR_full_subregion,
                                 nrow=3, ncol=2,
                                 heights = c(1,1,1.15),
                                 labels = c('A', 'D', 'B', 'E', 'C', 'F'))

ggsave(plot = XDR_figure_arranged, "plots/XDR_final_fig.jpg", width=30, height=30, units="cm", dpi=300)
ggsave(plot = XDR_figure_arranged, "plots/XDR_final_fig.svg", width=30, height=30, units="cm", dpi=300)

### CTX-M figure ###

## phylogroup

# subsampled to min group size

pred_min_ceph_phylogroup <- ggpredict(min_ceph_phylogroup, terms = c("phylogroup"), type='fixed', ci.lvl = 0.95)

plot_ceph_min_phylogroup <- ggplot(pred_min_ceph_phylogroup, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Predicted probability\n of CTX-M") 

# full dataset

pred_full_ceph_phylogroup <- ggpredict(full_ceph_phylogroup, terms = c("phylogroup"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

plot_ceph_full_phylogroup <- ggplot(pred_full_ceph_phylogroup, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Predicted probability\n of CTX-M") 

## host

# subsampled to min group size

pred_ceph_min_host <- ggpredict(min_ceph_host, terms = c("host"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x)) %>%
  mutate(x = str_replace(x, "agricultural/domestic", "Agricultural/domestic"))  %>%
  mutate(x = str_replace(x, "human", "Human")) %>%
  mutate(x = str_replace(x, "wild animal", "Wild animal")) %>%
  mutate(x = str_replace(x, "wild bird", "Wild bird"))

plot_ceph_min_host <- ggplot(pred_ceph_min_host, aes(x = x, y = predicted, colour=x))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Host species", y = "Predicted probability\n of CTX-M") +
  scale_colour_manual(values=met.brewer("Juarez")) +
  scale_x_discrete(name = NULL,
                   labels = host_pics) +
  theme(axis.text.x = ggtext::element_markdown(size = 12))

# full dataset

pred_ceph_full_host <- ggpredict(full_ceph_host, terms = c("host"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x)) %>%
  mutate(x = str_replace(x, "agricultural/domestic", "Agricultural/domestic"))  %>%
  mutate(x = str_replace(x, "human", "Human")) %>%
  mutate(x = str_replace(x, "wild animal", "Wild animal")) %>%
  mutate(x = str_replace(x, "wild bird", "Wild bird"))

plot_ceph_full_host <- ggplot(pred_ceph_full_host, aes(x = x, y = predicted, colour=x))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Host species", y = "Predicted probability\n of CTX-M") +
  scale_colour_manual(values=met.brewer("Juarez")) +
  scale_x_discrete(name = NULL,
                   labels = host_pics) +
  theme(axis.text.x = ggtext::element_markdown(size = 12))

## subregion

# subsampled to min group size

pred_ceph_min_subregion <- ggpredict(min_ceph_subregion, terms = c("subregion"), type='fixed', ci.lvl = 0.95) %>%
  rename(subregion = x) %>%
  dplyr::select(-group)

map_ceph_min_subregion <- world %>%
  full_join(pred_ceph_min_subregion)

plot_ceph_min_subregion <- ggplot(data = map_ceph_min_subregion) +
  geom_sf(aes(fill = predicted)) +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), name = "Predicted probability\n of CTX-M") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=15), 
        legend.text=element_text(size=10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = "bottom",
        plot.margin = unit(c(1, -5, 1, -5), "cm"))

# full dataset

pred_ceph_full_subregion <- ggpredict(full_ceph_subregion, terms = c("subregion"), type='fixed', ci.lvl = 0.95) %>%
  rename(subregion = x) %>%
  dplyr::select(-group)

map_ceph_full_subregion <- world %>%
  full_join(pred_ceph_full_subregion)

plot_ceph_full_subregion <- ggplot(data = map_ceph_full_subregion) +
  geom_sf(aes(fill = predicted)) +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), name = "Predicted probability\n of CTX-M") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=15), 
        legend.text=element_text(size=10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = "bottom",
        plot.margin = unit(c(1, -5, 1, -5), "cm"))

## arrange

# 1 row = 1 predictor, left column is full dataset, right column is min dataset

ceph_figure_arranged <- ggarrange(plot_ceph_min_phylogroup, plot_ceph_full_phylogroup,
                                 plot_ceph_min_host, plot_ceph_full_host,
                                 plot_ceph_min_subregion, plot_ceph_full_subregion,
                                 nrow=3, ncol=2,
                                 heights = c(1,1,1.15),
                                 labels = c('A', 'D', 'B', 'E', 'C', 'F'))

ggsave(plot = ceph_figure_arranged, "plots/ceph_final_fig.jpg", width=30, height=30, units="cm", dpi=300)
ggsave(plot = ceph_figure_arranged, "plots/ceph_final_fig.svg", width=30, height=30, units="cm", dpi=300)

### carbapenem res. figure ###

## phylogroup

# subsampled to min group size

pred_min_carb_phylogroup <- ggpredict(min_carb_phylogroup, terms = c("phylogroup"), type='fixed', ci.lvl = 0.95)

plot_carb_min_phylogroup <- ggplot(pred_min_carb_phylogroup, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Predicted probability\n of carbapenem res.") 

# full dataset

pred_full_carb_phylogroup <- ggpredict(full_carb_phylogroup, terms = c("phylogroup"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x))

plot_carb_full_phylogroup <- ggplot(pred_full_carb_phylogroup, aes(x = x, y = predicted))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Phylogroup", y = "Predicted probability\n of carbapenem res.") 

## host

# subsampled to min group size

pred_carb_min_host <- ggpredict(min_carb_host, terms = c("host"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x)) %>%
  mutate(x = str_replace(x, "agricultural/domestic", "Agricultural/domestic"))  %>%
  mutate(x = str_replace(x, "human", "Human")) %>%
  mutate(x = str_replace(x, "wild animal", "Wild animal")) %>%
  mutate(x = str_replace(x, "wild bird", "Wild bird"))

plot_carb_min_host <- ggplot(pred_carb_min_host, aes(x = x, y = predicted, colour=x))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Host species", y = "Predicted probability\n of carbapenem res.") +
  scale_colour_manual(values=met.brewer("Juarez")) +
  scale_x_discrete(name = NULL,
                   labels = host_pics) +
  theme(axis.text.x = ggtext::element_markdown(size = 12))

# full dataset

pred_carb_full_host <- ggpredict(full_carb_host, terms = c("host"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x = as.character(x)) %>%
  mutate(x = str_replace(x, "agricultural/domestic", "Agricultural/domestic"))  %>%
  mutate(x = str_replace(x, "human", "Human")) %>%
  mutate(x = str_replace(x, "wild animal", "Wild animal")) %>%
  mutate(x = str_replace(x, "wild bird", "Wild bird"))

plot_carb_full_host <- ggplot(pred_carb_full_host, aes(x = x, y = predicted, colour=x))+ 
  theme_light() +
  theme(text = element_text(size=15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  geom_point(size=4, position = position_dodge(width=0.8)) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.8, position = position_dodge(width=0.8)) +  
  labs(x = "Host species", y = "Predicted probability\n of carbapenem res.") +
  scale_colour_manual(values=met.brewer("Juarez")) +
  scale_x_discrete(name = NULL,
                   labels = host_pics) +
  theme(axis.text.x = ggtext::element_markdown(size = 12))

## arrange

# 1 row = 1 predictor, left column is full dataset, right column is min dataset

carb_figure_arranged <- ggarrange(plot_carb_min_phylogroup, plot_carb_full_phylogroup,
                                  plot_carb_min_host, plot_carb_full_host,
                                  nrow=2, ncol=2,
                                  labels = c('A', 'C', 'B', 'D'))

ggsave(plot = carb_figure_arranged, "plots/carb_final_fig.jpg", width=30, height=20, units="cm", dpi=300)
ggsave(plot = carb_figure_arranged, "plots/carb_final_fig.svg", width=30, height=30, units="cm", dpi=300)
