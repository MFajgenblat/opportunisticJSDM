#-------------------------------------------------------------------------------
# Setting working directory
#-------------------------------------------------------------------------------

# Choose the top-level directory
setwd("~/opportunisticJSDM")

#-------------------------------------------------------------------------------
# Load packages
#-------------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(RANN)
library(lubridate)
library(mgcv)
library(splines)
library(tidybayes)
library(posterior)
library(doParallel)
library(ape)
library(cmdstanr)
library(reshape2)
library(patchwork)
library(ggtree)
library(ggstats)
library(magick)
library(ggimage)
library(ggtext)
library(abind)

#-------------------------------------------------------------------------------
# Read data
#-------------------------------------------------------------------------------

waterbodies_metadata <- read.csv("F3_Application/Data/waterbodies_metadata.csv")
species_metadata <- read.csv("F3_Application/Data/species_metadata.csv", sep = ";") %>%
  arrange(Species) %>%
  mutate(alphabetical_id = row_number())
waterbodies <- read.csv("F3_Application/Data/waterbodies_metadata.csv", sep = ",") %>%
  st_as_sf(coords = c("X", "Y"), crs = 31370)
datalist <- readRDS("F3_Application/Data/datalist.rds")

#-------------------------------------------------------------------------------
# Generate taxonomic tree plot
#-------------------------------------------------------------------------------

transparent <- function(img) {
  magick::image_fx(img, expression = "0.75*a", channel = "alpha")
}

tree <- species_metadata %>%
  mutate(Species = factor(Species),
         Genus = factor(Genus),
         Family = factor(Family),
         Superfamily = factor(Superfamily),
         Suborder = factor(Suborder),
         Order = factor(Order)) %>%
  as.phylo(~ Order/Suborder/Superfamily/Family/Genus/Species, data = ., collapse = FALSE)
tree$edge.length = rep(1, length(tree$edge))

treeplot <- species_metadata %>%
  mutate(Species = factor(paste0("   ", Species)),
         Genus = factor(Genus),
         Family = factor(Family),
         Superfamily = factor(Superfamily),
         Suborder = factor(Suborder),
         Order = factor(Order)) %>%
  as.phylo(~ Order/Suborder/Superfamily/Family/Genus/Species, data = ., collapse = FALSE) %>%
  ggtree(ladderize = F, size = 0.2) +
  geom_image(aes(subset=(node==81), image="F3_Application/Species_silhouettes/Lestes viridis.svg"), size = .06, color = "#1B9E77", image_fun = transparent, nudge_x = -1) +
  geom_image(aes(subset=(node==73), image="F3_Application/Species_silhouettes/Ischnura elegans.svg"), size = .06, color = "#D95F02", image_fun = transparent, nudge_x = -1) +
  geom_image(aes(subset=(node==65), image="F3_Application/Species_silhouettes/Sympetrum striolatum.svg"), size = .06, color = "#7570B3", image_fun = transparent, nudge_x = -1) +
  geom_image(aes(subset=(node==62), image="F3_Application/Species_silhouettes/Somatochlora metallica.svg"), size = .06, color = "#E7298A", image_fun = transparent, nudge_x = -1) +
  geom_image(aes(subset=(node==59), image="F3_Application/Species_silhouettes/Gomphus pulchellus.svg"), size = .06, color = "#66A61E", image_fun = transparent, nudge_x = -1) +
  geom_image(aes(subset=(node==55), image="F3_Application/Species_silhouettes/Anax imperator.svg"), size = .06, color = "#E6AB02", image_fun = transparent, nudge_x = -1) +
  geom_hilight(node=81, fill="#1B9E77", alpha=0.15) +
  geom_hilight(node=73, fill="#D95F02", alpha=0.15) +
  geom_hilight(node=65, fill="#7570B3", alpha=0.15) +
  geom_hilight(node=62, fill="#E7298A", alpha=0.15) +
  geom_hilight(node=59, fill="#66A61E", alpha=0.15) +
  geom_hilight(node=55, fill="#E6AB02", alpha=0.15) +
  scale_x_continuous(limits = c(NA,NA), expand = c(0,0)) +
  scale_y_reverse("Species", expand = c(0,0.0)) +
  coord_fixed(ratio = 2, clip = "off") +
  #geom_tiplab(size = 1.65, fontface = 3, color = "grey25") +
  theme(axis.title.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))
treeplot

#-------------------------------------------------------------------------------
# Read posterior
#-------------------------------------------------------------------------------

fit_constrained <- readRDS("F3_Application/Fits/fit.rds")
#fit_unconstrained <- readRDS("F3_Application/Fits/fit_unconstrained.rds")

#-------------------------------------------------------------------------------
# Preparing metadata
#-------------------------------------------------------------------------------

species_names <- gsub("_", " ", gsub("Sp.", "", colnames(datalist$any_seen)))
day_names <- seq(as.Date("2000-01-01"), as.Date("2000-12-31"), by = 1)
year_names <- 2009:2023
psi_covariate_names <- c("Intercept", "Pond area", "Alluvial forest", "Brownfield", "Coastal habitat", "Deciduous forest", "Heathland", "Pine forest", "Plantation forest", "Semi-natural grassland", "Shrubland", "Urban area", "Wetland and open water")
p_covariate_names <- gsub("_", " ", colnames(datalist$detection_X))
trait_names <- gsub(" scaled", "", gsub("Migratory", "Good disperser", gsub("_", " ", rownames(datalist$Tr))))

#-------------------------------------------------------------------------------
# Phenological effects
#-------------------------------------------------------------------------------

Sys.setlocale("LC_ALL", "English")
expit <- function(x) {1/(1+exp(-x))}

phenology <- melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 10) == "phenology[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    week = 1:datalist$N_weeks,
                    Species = species_names)[,-1]) %>%
  mutate(phenology_prob = expit(value))

phenology %>%
  group_by(Species, draw) %>%
  summarise(peak = which.max(phenology_prob)) %>%
  group_by(Species) %>%
  summarise(peak = mean(peak)) %>%
  arrange(-peak) %>%
  pull(Species) -> phenology_rank

p_phenology_summarized <- phenology %>%
  mutate(Species = factor(Species, levels = phenology_rank)) %>%
  group_by(Species, week) %>%
  summarise(p_phenology = median(phenology_prob))

phenology_plot <- ggplot(p_phenology_summarized) +
  geom_tile(aes(x = seq(as.Date("2022-01-01"), as.Date("2022-12-31"), by = 1)[1+(week-1)*7], y = Species, fill = p_phenology), color = "white", linewidth = 0.1) +
  scale_x_date("Date", expand = c(0,0), date_breaks = "2 months", date_labels = "%b %d") +
  scale_y_discrete("Species", guide = guide_axis(n.dodge = 1), expand = c(0,0)) +
  scale_fill_gradientn("Posterior median phenological intensity", 
                       colors = c("#FEF4DE","#FDEDC8","#FDE5B2","#fcde9c","#faa476","#f0746e","#e34f6f","#dc3977","#b9257a","#7c1d6f"),
                       limits = c(0,1), guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  ggtitle("(a) Phenological detection patterns") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 6),
        axis.title.x = element_text(size = 6, face = "bold", vjust = 2),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 4, face = "italic"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5),
        legend.position = "bottom",
        aspect.ratio = 1,
        plot.margin = unit(c(0, 0.05, 0, 0), "cm"))

#-------------------------------------------------------------------------------
# Observer effects
#-------------------------------------------------------------------------------

Eta_observer <- melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 13) == "observer_Eta[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    observer = 1:datalist$N_observers,
                    dim = 1:datalist$N_observer_dims)[,-1]) %>%
  group_by(observer, dim) %>%
  summarize(value = median(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = dim, values_from = value) %>%
  select(-observer) %>%
  as.matrix()

observer_Lambda <- melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 16) == "observer_Lambda[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    dim = 1:datalist$N_observer_dims,
                    species = 1:datalist$N_species)[,-1]) %>%
  group_by(species, dim) %>%
  summarize(value = median(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = dim, values_from = value) %>%
  select(-species) %>%
  as.matrix()

observer_effects <- Eta_observer %*% t(observer_Lambda)

observer_plot <- melt(observer_effects, varnames = c("Observer", "Species")) %>%
  mutate(Species = factor(species_names[Species], levels = species_names[order(prcomp(observer_effects)$rotation[,1])])) %>%
  ggplot() +
  geom_tile(aes(x = factor(Observer, levels = order(prcomp(observer_effects)$x[,1])), y = Species, fill = value)) +
  scale_x_discrete("Observer (anonymized)", expand = c(0,0)) +
  scale_y_discrete("Species", guide = guide_axis(n.dodge = 1), expand = c(0,0)) +
  scale_fill_gradient2("Posterior median observer random effect",
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296",
                       limits = c(-max(abs(observer_effects)), max(abs(observer_effects))),
                       guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  ggtitle("(b) Observer-specific detection patterns") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(face = "bold", size = 6),
        axis.title.x = element_text(size = 6, face = "bold", vjust = 2),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 4, face = "italic"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5),
        legend.position = "bottom",
        plot.margin = unit(c(0, 0, 0, 0.05), "cm"))

#-------------------------------------------------------------------------------
# Phenological and observer effects
#-------------------------------------------------------------------------------

phenology_plot + observer_plot + plot_layout(ncol = 2)
ggsave("Figure_1.png", width = 17.3, height = 9, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Detection covariate (and trait) influences
#-------------------------------------------------------------------------------

treeplot +
  # Species-specific effects
  melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 15) == "detection_beta[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    cov = 1:datalist$N_detection_covs,
                    Species = species_names)[,-1]) %>%
  mutate(Species = factor(Species, levels = rev(tree$tip.label)),
         cov = factor(p_covariate_names[cov], levels = p_covariate_names)) %>%
  group_by(Species, cov) %>%
  summarize(important = abs(mean(value > 0)*2-1) > 0.95,
            value = median(value)) %>%
  ggplot() +
  geom_tile(aes(x = cov, y = Species, fill = value), color = "white", linewidth = 0.2) +
  geom_point(aes(x = cov, y = Species, shape = important), size = 0.75, alpha = 0.5) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2("Posterior\nmedian\neffect",
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                       guide = guide_colorbar(barwidth = unit(0.2, "cm"), barheight = unit(5.5, "cm"), title.position = "top", title.hjust = 0, title.vjust = -2)) +
  scale_shape_manual(values = c(NA, 8), guide = "none") +
  coord_equal() +
  ggtitle("(a) Detection-related\npredictor associations") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(face = "bold", size = 6),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 0),
        axis.text.y = element_text(size = 4, face = "italic"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0.05), "cm")) +
  plot_spacer() +
  melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 12) == "detection_G[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    cov = 1:datalist$N_detection_covs,
                    trait = trait_names)[,-1]) %>%
  mutate(trait = factor(trait, levels = rev(trait_names)),
         cov = factor(p_covariate_names[cov], levels = p_covariate_names)) %>%
  group_by(trait, cov) %>%
  summarize(important = abs(mean(value > 0)*2-1) > 0.95,
            value = median(value)) %>%
  ggplot() +
  geom_tile(aes(x = cov, y = trait, fill = value), color = "white", linewidth = 0.2) +
  geom_point(aes(x = cov, y = trait, shape = important), size = 0.75, alpha = 0.5) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2("Posterior\nmedian\neffect",
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                       guide = "none") +
  scale_shape_manual(values = c(NA, 8), guide = "none") +
  coord_equal() +
  ggtitle("(b) Detection-related\ntrait influences") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(face = "bold", size = 5.8),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0.05), "cm")) +
  plot_layout(ncol = 2, guide = "collect")
ggsave("Detection_process.png", width = 8.2, height = 18, units = "cm", dpi = 600, bg = "white")

#-------------------------------------------------------------------------------
# Occupancy covariate (and trait) influences
#-------------------------------------------------------------------------------

treeplot +
  # Species-specific effects
  melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 15) == "occupancy_beta[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    cov = 1:datalist$N_occupancy_covs,
                    Species = species_names)[,-1]) %>%
  mutate(Species = factor(Species, levels = rev(tree$tip.label)),
         cov = factor(gsub("Pine forest", "Coniferous forest", psi_covariate_names)[cov], levels = gsub("Pine forest", "Coniferous forest", psi_covariate_names))) %>%
  group_by(Species, cov) %>%
  summarize(important = abs(mean(value > 0)*2-1) > 0.95,
            value = median(value)) %>%
  mutate(value = case_when(value > 5 ~ 5,
                           value < -5 ~ -5,
                           T ~ value)) %>%
  ggplot() +
  geom_tile(aes(x = cov, y = Species, fill = value), color = "white", linewidth = 0.2) +
  geom_point(aes(x = cov, y = Species, shape = important), size = 0.75, alpha = 0.5) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2("Posterior\nmedian\neffect",
                       breaks = seq(-5, 5, by = 2.5), labels = c("<-5", "-2.5", "0", "2.5", ">5"),
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                       guide = guide_colorbar(barwidth = unit(0.2, "cm"), barheight = unit(5.5, "cm"), title.position = "top", title.hjust = 0, title.vjust = -2)) +
  scale_shape_manual(values = c(NA, 8), guide = "none") +
  coord_equal() +
  ggtitle("(a) Occupancy-related predictor associations") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(face = "bold", size = 6),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 0),
        axis.text.y = element_text(size = 4, face = "italic"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0.05), "cm")) +
  plot_spacer() +
  melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 12) == "occupancy_G[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    cov = 1:datalist$N_occupancy_covs,
                    trait = trait_names)[,-1]) %>%
  mutate(trait = factor(trait, levels = rev(trait_names)),
         cov = factor(psi_covariate_names[cov], levels = psi_covariate_names)) %>%
  group_by(trait, cov) %>%
  summarize(important = abs(mean(value > 0)*2-1) > 0.95,
            value = median(value)) %>%
  mutate(value = case_when(value > 2 ~ 2,
                           value < -2 ~ -2,
                           T ~ value)) %>%
  ggplot() +
  geom_tile(aes(x = cov, y = trait, fill = value), color = "white", linewidth = 0.2) +
  geom_point(aes(x = cov, y = trait, shape = important), size = 0.75, alpha = 0.5) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2("Posterior\nmedian\neffect",
                       limits = c(-5,5),
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                       guide = "none") +
  scale_shape_manual(values = c(NA, 8), guide = "none") +
  coord_equal() +
  ggtitle("(b) Occupancy-related trait influences") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(face = "bold", size = 5.8),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0.05), "cm")) +
  plot_layout(ncol = 2, guide = "collect")
ggsave("Figure_2.png", width = 8.2, height = 17.5, units = "cm", dpi = 600, bg = "white")

#-------------------------------------------------------------------------------
# Omega matrix - spatial
#-------------------------------------------------------------------------------

temp1 <- melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 15) == "spatial_Lambda[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    dim = 1:datalist$N_spatial_dims,
                    species = species_names)[,-1]) %>%
  select(-variable)
temp2 <- lapply(unique(temp1$draw), function(i) as.matrix(select(pivot_wider(filter(temp1, draw == i), names_from = dim, values_from = value), -draw, -species)))
temp3 <- apply(abind(lapply(temp2, function(i) cov2cor(i %*% t(i))), along = 3), 1:2, median)
temp4 <- melt(fit_unconstrained[,which(substr(colnames(fit_unconstrained), 1, 15) == "spatial_Lambda[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    dim = 1:datalist$N_spatial_dims,
                    species = species_names)[,-1]) %>%
  select(-variable)
temp5 <- lapply(unique(temp4$draw), function(i) as.matrix(select(pivot_wider(filter(temp4, draw == i), names_from = dim, values_from = value), -draw, -species)))
temp6 <- apply(abind(lapply(temp5, function(i) cov2cor(i %*% t(i))), along = 3), 1:2, median)
rbind(melt(temp3) %>%
        filter(Var1 <= Var2),
      melt(temp6) %>%
        filter(Var1 > Var2)) %>%
  mutate(Species_from = factor(species_names[Var1], levels = (species_names[order(prcomp(temp3)$rotation[,1])])),
         Species_to = factor(species_names[Var2], levels = rev(species_names[order(prcomp(temp3)$rotation[,1])]))) %>%
  ggplot() +
  geom_tile(aes(x = Species_from, y = Species_to, fill = value), color = "white", linewidth = 0.1) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2("Posterior\nmedian\ninterspecific\ncorrelation",
                       limits = c(-1,1),
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                       guide = guide_colorbar(barwidth = unit(0.2, "cm"), barheight = unit(5.5, "cm"), title.position = "top", title.hjust = 0)) +
  coord_equal() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(face = "bold", size = 6),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 4, face = "italic", angle = 75, hjust = 0),
        axis.text.y = element_text(size = 4, face = "italic"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5),
        legend.position = "right")
ggsave("Figure_3.png", width = 11, height = 9.5, units = "cm", dpi = 600, bg = "white")

#-------------------------------------------------------------------------------
# Omega matrix - spatiotemporal
#-------------------------------------------------------------------------------

temp1 <- melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 22) == "spatiotemporal_Lambda[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    dim = 1:datalist$N_spatiotemporal_dims,
                    species = species_names)[,-1]) %>%
  select(-variable)
temp2 <- lapply(unique(temp1$draw), function(i) as.matrix(select(pivot_wider(filter(temp1, draw == i), names_from = dim, values_from = value), -draw, -species)))
temp3 <- apply(abind(lapply(temp2, function(i) cov2cor(i %*% t(i))), along = 3), 1:2, median)
temp4 <- melt(fit_unconstrained[,which(substr(colnames(fit_unconstrained), 1, 22) == "spatiotemporal_Lambda[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    dim = 1:datalist$N_spatiotemporal_dims,
                    species = species_names)[,-1]) %>%
  select(-variable)
temp5 <- lapply(unique(temp4$draw), function(i) as.matrix(select(pivot_wider(filter(temp4, draw == i), names_from = dim, values_from = value), -draw, -species)))
temp6 <- apply(abind(lapply(temp5, function(i) cov2cor(i %*% t(i))), along = 3), 1:2, median)
rbind(melt(temp3) %>%
        filter(Var1 <= Var2),
      melt(temp6) %>%
        filter(Var1 > Var2)) %>%
  mutate(Species_from = factor(species_names[Var1], levels = (species_names[order(prcomp(temp3)$rotation[,1])])),
         Species_to = factor(species_names[Var2], levels = rev(species_names[order(prcomp(temp3)$rotation[,1])]))) %>%
  ggplot() +
  geom_tile(aes(x = Species_from, y = Species_to, fill = value), color = "white", linewidth = 0.1) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2("Posterior\nmedian\ninterspecific\ncorrelation",
                       limits = c(-1,1),
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                       guide = guide_colorbar(barwidth = unit(0.2, "cm"), barheight = unit(5.5, "cm"), title.position = "top", title.hjust = 0)) +
  coord_equal() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(face = "bold", size = 6),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 4, face = "italic", angle = 75, hjust = 0),
        axis.text.y = element_text(size = 4, face = "italic"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5),
        legend.position = "right")
ggsave("Figure_S6.png", width = 11, height = 9.5, units = "cm", dpi = 600, bg = "white")

#-------------------------------------------------------------------------------
# Add data to the Stan datalist for prediction purposes
#-------------------------------------------------------------------------------

inputdata <- fit_constrained
#inputdata[1,] <- apply(inputdata, 2, median)
#inputdata <- inputdata[1,,drop=F]

all_units <- merge(waterbodies_metadata, data.frame(year = 2009:2023))

datalist <- readRDS("F3_Application/Data/datalist.rds")
datalist_pred <- datalist
datalist_pred$N_units_pred <- nrow(all_units)
datalist_pred$N_sites_pred <- length(unique(all_units$waterbody))
datalist_pred$occupancy_X_pred <- occupancy_X
datalist_pred$spatial_bf_pred <- spatial_bf
datalist_pred$spatiotemporal_bf_spatial_pred <- spatiotemporal_bf_spatial
datalist_pred$spatiotemporal_bf_spatial_pred <- spatiotemporal_bf_spatial
datalist_pred$site_pred <- all_units$waterbody
datalist_pred$year_pred <- all_units$year - 2008
datalist_pred$included <- (waterbodies_metadata$waterbody %in% visitmetadata$waterbody)*1
datalist_pred$site_crossid <- distinct(data.frame(waterbody = visitmetadata$waterbody, site_crossid = as.numeric(droplevels(visitmetadata$waterbody)))) %>%
  tidyr::complete(waterbody, fill = list(site_crossid = 0)) %>%
  pull(site_crossid)

#-------------------------------------------------------------------------------
# Species-specific distribution
#-------------------------------------------------------------------------------

datalist_pred$species_pred <- 4
mod_gq <- cmdstan_model("F3_Application/opportunisticJSDM_generating_occupancies.stan")
fit_gq <- mod_gq$generate_quantities(inputdata, data = datalist_pred, sig_figs = 4)
fit_gqs <- read_cmdstan_csv(fit_gq$output_files(), format = "draws_matrix")$post_warmup_draws
temp <- as_draws_df(fit_gq$draws())
temp2 <- data.frame(all_units, Z = apply(temp[,1:nrow(all_units)], 2, median))
temp2 %>%
  st_as_sf(coords = c("X", "Y"), crs = 31370) %>%
  ggplot() +
  geom_sf(aes(color = Z), shape = 16, size = 0.01) +
  scale_color_gradientn("Posterior median occupancy probability",
                        limits = c(0,1), breaks = seq(0, 1, by = 0.2),
                        colors = c("#FEF4DE","#FDEDC8","#FDE5B2","#fcde9c","#faa476","#f0746e","#e34f6f","#dc3977","#b9257a","#7c1d6f"),
                        guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  facet_wrap(~ year) +
  theme_void() +
  theme(plot.title = element_markdown(size = 6, hjust = 0.08),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0))
ggsave("Figure_S7.png", width = 17.3, height = 9, units = "cm", dpi = 600, bg = "white")

# Plot B: yearly temporal trend
inputdata <- fit_constrained
datalist_pred$species_pred <- 4
mod_gq <- cmdstan_model("F3_Application/opportunisticJSDM_generating_occupancies.stan")
fit_gq <- mod_gq$generate_quantities(inputdata, data = datalist_pred, sig_figs = 4)
fit_gqs <- read_cmdstan_csv(fit_gq$output_files(), format = "draws_matrix")$post_warmup_draws
temp <- as_draws_df(fit_gq$draws())

split_data <- split(t(as.matrix(temp[,1:nrow(all_units)])), all_units$year)
result <- do.call(rbind, lapply(split_data, function(x) colMeans(matrix(x, ncol=nrow(temp)))))
melt(result) %>%
  ggplot() +
  stat_lineribbon(aes(x = Var1, y = value), .width = 0.5, color = "#274a59", fill = "#274a5950", size = 0.5) +
  scale_x_continuous("Year", expand = c(0,0), breaks = seq(2009, 2023, by = 2)) +
  scale_y_continuous("Fraction of\noccupiedponds", limits = c(0, NA), expand = c(0,0)) +
  ggtitle("<b>(b) Temporal trend of <i>Aeshna isoceles</i></b>") +
  theme(plot.title = element_markdown(size = 6, hjust = 0.08),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.title = element_text(face = "bold", size = 6),
        axis.text = element_text(size = 5))
ggsave("Occupancies_B.png", width = 10.3/2, height = 3.5, units = "cm", dpi = 600, bg = "white")

# Plot C: trend slope histogram

temp5 <- melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 15) == "temporal_trend[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    year = 1:datalist$N_years,
                    Species = species_names)[,-1]) %>%
  select(-variable) %>%
  pivot_wider(names_from = year, values_from = value)
year_range <- seq(0, 1, length.out = 15)
temp6 <- data.frame(select(temp5, draw, Species),
                    Slope = apply(as.matrix(select(temp5, -draw, -Species)), 1, function(y) lm(y ~ year_range)$coef[2]))
temp6 %>%
  ggplot() +
  stat_interval(aes(x = Slope, y = Species), .width = c(0.5,0.8,0.95,0.99)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer()

temp6 %>%
  group_by(Species) %>%
  summarize(Category = case_when(mean(Slope > 0) > 0.95 ~ "Increasing",
                                 mean(Slope < 0) > 0.95 ~ "Decreasing",
                                 T ~ "Other")) %>%
  group_by(Category) %>%
  tally()

temp6 %>%
  group_by(Species) %>%
  summarize(Slope = median(Slope)) %>%
  ggplot() +
  geom_histogram(aes(x = Slope), color = "#274a59", fill = "#274a5950", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous("Posterior median linearized trend slope", expand = c(0,0)) +
  scale_y_continuous("Fraction of\noccupiedponds", limits = c(0,NA), expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  ggtitle("(c) Temporal trend slopes across species") +
  theme(plot.title = element_text(face = "bold", size = 6, hjust = 0.08),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0.5, "mm"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_blank())
ggsave("Occupancies_C.png", width = 10.3/2, height = 3.5, units = "cm", dpi = 600, bg = "white")
  
temp6 %>%
  group_by(Species) %>%
  summarize(Slope = median(Slope)) %>%
  left_join(species_metadata %>%
              ungroup() %>%
              select(Species, Body_size_scaled, Migratory, Permanent_ponds, Nutrient_level, Antipredation_index, STI)) %>%
  pivot_longer(!c(Species, Slope)) %>%
  mutate(name = case_when(name == "Antipredation_index" ~ "Antipredation index",
                          name == "Body_size_scaled" ~ "Body size",
                          name == "Migratory" ~ "Good disperser",
                          name == "Nutrient_level" ~ "Nutrient level",
                          name == "Permanent_ponds" ~ "Permanent ponds",
                          name == "STI" ~ "STI")) %>%
  ggplot(aes(x = value, y = Slope)) +
  geom_point(shape = 16, size = 0.25, color = "#274a59") +
  geom_smooth(method = "lm", color = "#274a59", fill = "#274a5950", linewidth = 0.2) +
  scale_x_continuous("") +
  scale_y_continuous("Posterior median linearized trend slope") +
  facet_wrap(~ name, scales = "free_x", switch = "x") +
  theme(plot.title = element_text(face = "bold", size = 6, hjust = 0),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.text = element_text(size = 5),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6))
ggsave("Figure_S8.png", width = 17.3, height = 10, units = "cm", dpi = 600, bg = "white")

#-------------------------------------------------------------------------------
# Temporal trend of all species
#-------------------------------------------------------------------------------

for (i in 1:50) {
  print(i)
  inputdata <- fit_constrained
  datalist_pred$species_pred <- i
  set_cmdstan_path("C:/Users/maxim/Documents/cmdstan-2.30.1")
  mod_gq <- cmdstan_model("F3_Application/opportunisticJSDM_gq_occupancy.stan")
  fit_gq <- mod_gq$generate_quantities(inputdata, data = datalist_pred, sig_figs = 4)
  fit_gqs <- read_cmdstan_csv(fit_gq$output_files(), format = "draws_matrix")$post_warmup_draws
  temp <- as_draws_df(fit_gq$draws())
  
  split_data <- split(t(as.matrix(temp[,1:nrow(all_units)])), all_units$year)
  result <- do.call(rbind, lapply(split_data, function(x) colMeans(matrix(x, ncol=nrow(temp)))))
  write.table(data.frame(species = i, melt(result)), "species_trends.csv", row.names = F, col.names = F, sep = ",", append = T)
}

read.csv("species_trends.csv", header=FALSE) %>%
  ggplot() +
  stat_lineribbon(aes(x = V2, y = V4), .width = 0.95, color = "#274a59", fill = "#274a5950", size = 0.5) +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous("Fraction of occupied waterbodies", limits = c(NA, NA), expand = c(0,0)) +
  facet_wrap(~ species_names[V1], scales = "free", ncol = 6) +
  theme(plot.title = element_markdown(size = 6, hjust = 0.08),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93", size = 0.25),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 5),
        strip.clip = "off",
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.title = element_text(face = "bold", size = 6),
        axis.text = element_text(size = 5))
ggsave("Figure_4.png", width = 17.3, height = 18, units = "cm", dpi = 600, bg = "white")

temp15 <- read.csv("species_trends.csv", header=FALSE) %>%
  group_by(V1, V2, V3) %>%
  summarize(V4 = median(V4)) %>%
  pivot_wider(names_from = V2, values_from = V4) %>%
  ungroup()
year_range <- seq(0, 1, length.out = 15)
temp6 <- data.frame(select(temp15, V1, V3),
                    Slope = apply(as.matrix(select(temp15, -V1, -V3)), 1, function(y) lm(y ~ year_range)$coef[2]))
temp6 %>%
  ggplot() +
  stat_interval(aes(x = Slope, y = species_names[V1]), .width = c(0.5,0.8,0.95,0.99)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer()

temp6 %>%
  group_by(V1) %>%
  summarize(Category = case_when(mean(Slope > 0) > 0.95 ~ "Increasing",
                                 mean(Slope < 0) > 0.95 ~ "Decreasing",
                                 T ~ "Other")) %>%
  group_by(Category) %>%
  tally()

temp6 %>%
  mutate(Species = species_names[V1]) %>%
  group_by(Species) %>%
  summarize(Slope = median(Slope)) %>%
  left_join(species_metadata %>%
              ungroup() %>%
              select(Species, Body_size_scaled, Migratory, Permanent_ponds, Nutrient_level, Antipredation_index, STI)) %>%
  pivot_longer(!c(Species, Slope)) %>%
  mutate(name = case_when(name == "Antipredation_index" ~ "Antipredation index",
                          name == "Body_size_scaled" ~ "Body size",
                          name == "Migratory" ~ "Good disperser",
                          name == "Nutrient_level" ~ "Nutrient level",
                          name == "Permanent_ponds" ~ "Permanent ponds",
                          name == "STI" ~ "STI")) %>%
  ggplot(aes(x = value, y = Slope)) +
  geom_point(shape = 16, size = 0.25, color = "#274a59") +
  geom_smooth(method = "lm", color = "#274a59", fill = "#274a5950", linewidth = 0.2) +
  scale_x_continuous("") +
  scale_y_continuous("Posterior median linearized trend slope") +
  facet_wrap(~ name, scales = "free_x", switch = "x") +
  theme(plot.title = element_text(face = "bold", size = 6, hjust = 0),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.text = element_text(size = 5),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6))
ggsave("Figure_S8.png", width = 17.3, height = 10, units = "cm", dpi = 600, bg = "white")

#-------------------------------------------------------------------------------
# Species richness LIJKT KAPOT TE ZIJN
#-------------------------------------------------------------------------------

datalist_pred$species_pred <- 0
set_cmdstan_path("C:/Users/maxim/Documents/cmdstan-2.30.1")
mod_gq <- cmdstan_model("F3_Application/opportunisticJSDM_gq_occupancy.stan")
fit_gq <- mod_gq$generate_quantities(inputdata, data = datalist_pred, sig_figs = 4)
fit_gqs <- read_cmdstan_csv(fit_gq$output_files(), format = "draws_matrix")$post_warmup_draws
temp <- as_draws_df(fit_gq$draws())
temp2 <- data.frame(all_units, Z = apply(temp[,1:nrow(all_units)], 2, median))

# Plot A: yearly species richness maps
plotdata <- temp2 %>%
  filter(year %in% seq(2009, 2023, by = 2)) %>%
  mutate(X = X - min(X),
         Y = Y - min(Y),
         X = X - (max(X)/2),
         X = X * (1 - 0.000005*Y),
         Y = Y/3) %>%
  mutate(Z = case_when(Z < 10 ~ 10,
                       Z > 25 ~ 25,
                       T ~ Z)) %>%
  arrange(year, Z)
ggplot() +
  geom_point(data = plotdata %>%
               filter(Y >= 10000),
             aes(x = X, y = Y + 3*5850*(year - 2008), color = Z), shape = 16, size = 0.05) +
  geom_segment(aes(x = -83000, xend = -83000, y = 15000, yend = 290000),
               arrow = arrow(length = unit(1.5, "mm")), size = 0.2) +
  geom_point(data = plotdata %>%
               filter(Y < 10000),
             aes(x = X, y = Y + 3*5850*(year - 2008), color = Z), shape = 16, size = 0.05) +
  geom_text(data = data.frame(label = seq(2009, 2023, by = 2), y = 24000+3*5850*(seq(2009, 2023, by = 2)-2008)), aes(x = -95000, y = y, label = label), size = 1.5) +
  scale_color_gradientn("Posterior median species richness",
                        breaks = seq(10, 25, by = 5),
                        labels = c("<10", "15", "20", ">25"),
                        colors = c("#FEF4DE","#FDEDC8","#FDE5B2","#fcde9c","#faa476","#f0746e","#e34f6f","#dc3977","#b9257a","#7c1d6f"),
                        guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  coord_equal() +
  theme_void() +
  ggtitle("(a) Yearly species richness maps") +
  theme(plot.title = element_text(face = "bold", size = 6, hjust = 0.08),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,0,0))
ggsave("SR_A.png", width = 7, height = 10.2, units = "cm", dpi = 600, bg = "white")

# Plot B: average species richness map
temp2 %>%
  group_by(waterbody, X, Y) %>%
  summarize(Z = mean(Z)) %>%
  arrange(Z) %>%
  mutate(Z = case_when(Z < 10 ~ 10,
                       Z > 25 ~ 25,
                       T ~ Z)) %>%
  st_as_sf(coords = c("X", "Y"), crs = 31370) %>%
  ggplot() +
  geom_sf(aes(color = Z), shape = 16, size = 0.075) +
  scale_color_gradientn("Posterior median species richness",
                        breaks = seq(10, 25, by = 5),
                        labels = c("<10", "15", "20", ">25"),
                        colors = c("#FEF4DE","#FDEDC8","#FDE5B2","#fcde9c","#faa476","#f0746e","#e34f6f","#dc3977","#b9257a","#7c1d6f"),
                        guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  ggtitle("(b) Average species richness map") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 6, hjust = 0.06, vjust = -6.5),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,0,0),
        plot.margin = unit(c(-5,0,-1,0), "mm"))
ggsave("SR_B.png", width = 10.3, height = 5.1, units = "cm", dpi = 600, bg = "white")

# Plot C: temporal species richness change
year_range <- seq(0, 1, length.out = 15)
temp3 <- temp2 %>%
  pivot_wider(names_from = year, values_from = Z)
temp4 <- data.frame(select(temp3, waterbody, WVLC, Area, X, Y),
                    Slope = apply(as.matrix(select(temp3, -waterbody, -WVLC, -Area, -X, -Y)), 1, function(y) lm(y ~ year_range)$coef[2]))
temp4 %>%
  arrange(abs(Slope)) %>%
  mutate(Slope = case_when(Slope > 3 ~ 3,
                           Slope < -3 ~ -3,
                           T ~ Slope)) %>%
  st_as_sf(coords = c("X", "Y"), crs = 31370) %>%
  ggplot() +
  geom_sf(aes(color = Slope), shape = 16, size = 0.075) +
  scale_color_gradient2("Posterior median change in species richness",
                        limits = c(-3,3),
                        low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                        guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  ggtitle("(c) Temporal change in species richness") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 6, hjust = 0.06, vjust = -6.5),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,0,0),
        plot.margin = unit(c(-5,0,-1,0), "mm"))
ggsave("SR_C.png", width = 10.3, height = 5.1, units = "cm", dpi = 600, bg = "white")

#-------------------------------------------------------------------------------
# Gradient plots
#-------------------------------------------------------------------------------

datalist_gradientpred <- datalist
resolution = 11
occupancy_X_mean <- cbind(t(as.matrix(colMeans(occupancy_X)))[rep(1,resolution),], Agriculture = 1 - sum(colMeans(occupancy_X)[3:13]))
occupancy_X_gradient <- rbind(data.frame(Variable = "Pond area", occupancy_X_mean[,1,drop=F], Area = seq(min(occupancy_X[,2]), max(occupancy_X[,2]), length.out = resolution), occupancy_X_mean[,3:14,drop=F]),
                              data.frame(Variable = "Alluvial forest", occupancy_X_mean[,1:2,drop=F], Alluvial.forest = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,4:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Brownfield", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Brownfield = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,5:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Coastal habitat", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:4,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Coastal.habitat = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,6:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Deciduous forest", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:5,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Deciduous.forest = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,7:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Heathland", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:6,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Heathland = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,8:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Pine forest", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:7,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Pine.forest = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,9:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Plantation forest", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:8,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Plantation.forest = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,10:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Semi-natural grassland", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:9,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Semi.natural.grassland = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,11:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Shrubland", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:10,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Shrubland = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,12:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Urban area", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:11,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Urban.area = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,13:14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Wetland and open water", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:12,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Wetland.and.open.water = seq(0, 1, length.out = resolution), sweep(occupancy_X_mean[,14,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*")),
                              data.frame(Variable = "Agriculture", occupancy_X_mean[,1:2,drop=F], sweep(occupancy_X_mean[,3:13,drop=F], 1, seq(1/(1 - mean(occupancy_X[,3])), 0, length.out = resolution), "*"), Agriculture = seq(0, 1, length.out = resolution)))
datalist_gradientpred$N_gradientpred <- nrow(occupancy_X_gradient)
datalist_gradientpred$occupancy_X_gradient <- occupancy_X_gradient[,2:14]

mod_gq <- cmdstan_model("F3_Application/opportunisticJSDM_generating_gradients.stan")
fit_gq <- mod_gq$generate_quantities(fit_constrained, data = datalist_gradientpred, sig_figs = 4)
fit_gqs <- read_cmdstan_csv(fit_gq$output_files(), format = "draws_matrix")$post_warmup_draws
temp <- as_draws_df(fit_gq$draws())
temp2 <- cbind(merge(data.frame(draw = 1:nrow(fit_constrained)), data.frame(Variable = occupancy_X_gradient$Variable, Gradient = seq(0, 1, length.out = resolution))), Z = melt(temp[,1:nrow(occupancy_X_gradient)])) %>%
  mutate(Variable = gsub("Pine forest", "Coniferous forest", Variable))

temp2 %>%
  mutate(Variable = factor(Variable, levels = temp2 %>%
                             filter(Gradient == 1) %>%
                             group_by(Variable) %>%
                             summarize(Z.value = median(Z.value)) %>%
                             arrange(Z.value) %>%
                             pull(Variable))) %>%
  ggplot() +
  stat_lineribbon(aes(x = Gradient, y = Z.value), .width = c(0.5,0.8,0.95,0.99), color = NA) +
  scale_x_continuous("Variable gradient", breaks = seq(0, 1, by = 0.2), expand = c(0,0)) +
  scale_y_continuous("Expected species richness", expand = c(0,0), limits = c(NA, NA)) +
  scale_fill_brewer("Credible\ninterval") +
  facet_wrap(~ Variable, ncol = 4, strip.position = "bottom", scales = "free_x") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.spacing.x = unit(5, "mm"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.text = element_text(size = 5),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5))
ggsave("Figure_S5.png", width = 17.3, height = 15, units = "cm", dpi = 600, bg = "white")


fit_constrained %>%
  spread_draws(temporal_trend[year,species]) %>%
  mutate(Species = species_names[species],
         year = year + 2008) %>%
  ggplot() +
  stat_lineribbon(aes(x = year, y = temporal_trend), .width = c(0.5,0.8,0.95,0.99), color = NA) +
  scale_fill_brewer("Credible\ninerval") +
  facet_wrap(~ species_names[species], scales = "free_y", ncol = 5) +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.spacing.x = unit(5, "mm"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.text = element_text(size = 5),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 6),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 5))
ggsave("Trends.png", width = 17.3, height = 20, units = "cm", dpi = 600, bg = "white")

melt(fit_constrained[,which(substr(colnames(fit_constrained), 1, 15) == "detection_beta[")]) %>%
  cbind(expand.grid(draw = 1:nrow(fit_constrained),
                    cov = 1:datalist$N_detection_covs,
                    Species = species_names)[,-1]) %>%
  mutate(Species = factor(Species, levels = rev(tree$tip.label)),
         cov = factor(c("Detection\nintercept", "Singleton list\neffect", "Short list\neffect")[cov],
                      levels = c("Detection\nintercept", "Singleton list\neffect", "Short list\neffect"))) %>%
  filter(cov == "Detection\nintercept") %>%
  mutate(value = expit(value)) %>%
  group_by(Species) %>%
  summarize(value = median(value)) %>%
  arrange(value) %>%
  View()

fit_constrained %>%
  spread_draws(occupancy_phylo_fraction) %>%
  summarize(Median = median(occupancy_phylo_fraction),
            Lquan = quantile(occupancy_phylo_fraction, 0.025),
            Uquan = quantile(occupancy_phylo_fraction, 0.975))

fit_constrained %>%
  spread_draws(detection_phylo_fraction) %>%
  summarize(Median = median(detection_phylo_fraction),
            Lquan = quantile(detection_phylo_fraction, 0.025),
            Uquan = quantile(detection_phylo_fraction, 0.975))

fit_constrained %>%
  spread_draws(spatial_structured_fraction[dim]) %>%
  ungroup() %>%
  summarize(Median = median(spatial_structured_fraction),
            Lquan = quantile(spatial_structured_fraction, 0.025),
            Uquan = quantile(spatial_structured_fraction, 0.975))

rbind(visitmetadata %>%
  group_by(year) %>%
  summarize(n = sum(N_visits)) %>%
  mutate(Outcome = "Number of visits"),
  visitdata %>%
    select(year, scrambled_user_id) %>%
    distinct() %>%
    group_by(year) %>%
    tally() %>%
    mutate(Outcome = "Number of observers")) %>%
  ggplot() +
  geom_bar(aes(x = year, y = n), stat = "identity", color = "#274a59", fill = "#274a5950", linewidth = 0.2) +
  scale_x_discrete("Year", expand = c(0,0)) +
  scale_y_continuous("Total number of visits", limits = c(0,NA), expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ Outcome, scales = "free_y", switch = "y") +
  theme(plot.title = element_text(face = "bold", size = 6, hjust = 0.08),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.y = element_line(color = "grey93"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        strip.placement = "outside",
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0.5, "mm"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5))
ggsave("Figure_S9.png", width = 17.3, height = 7, units = "cm", dpi = 600, bg = "white")
