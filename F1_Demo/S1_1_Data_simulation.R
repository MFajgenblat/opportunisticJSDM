#-------------------------------------------------------------------------------
# Setting working directory
#-------------------------------------------------------------------------------

# Choose the top-level directory
setwd("~/opportunisticJSDM")

#-------------------------------------------------------------------------------
# Loading packages
#-------------------------------------------------------------------------------

library(tidyverse)
library(ape)
library(mgcv)
library(splines)
library(cmdstanr)
library(tidybayes)
library(pROC)
library(ggtext)

#-------------------------------------------------------------------------------
# Simulation settings
#-------------------------------------------------------------------------------

N_species <- 50
N_sites <- 100
N_years <- 10
N_weeks <- 53
N_observers <- 50
N_allvisits <- 500
N_occupancy_covs <- 6
N_detection_covs <- 3
N_traits <- 3

#-------------------------------------------------------------------------------
# Generate site visits
#-------------------------------------------------------------------------------

visitdata_template <- data.frame(site = factor(sample(1:N_sites, N_allvisits, replace = T), 1:N_sites),
                                 year = factor(sample(1:N_years, N_allvisits, replace = T, prob = 1:N_years), 1:N_years),
                                 week = sample(1:N_weeks, N_allvisits, replace = T),
                                 observer = factor(sample(1:N_observers, N_allvisits, replace = T), levels = 1:N_observers)) %>%
  arrange(site, year, week, observer) %>%
  rowid_to_column(var = "visitid")

#-------------------------------------------------------------------------------
# Compute cyclic spline basis functions for weekly phenology generation
#-------------------------------------------------------------------------------

phenology_bf <- cSplineDes(1:53, seq(1, 53, length.out = 11))

#-------------------------------------------------------------------------------
# Generating an environmental data matrix
#-------------------------------------------------------------------------------

occupancy_X <- cbind(1, replicate(N_occupancy_covs - 1, rnorm(1:N_sites)))

#-------------------------------------------------------------------------------
# Generating a detection covariate matrix
#-------------------------------------------------------------------------------

detection_X <- cbind(1, replicate(N_detection_covs - 1, rnorm(1:N_allvisits)))

#-------------------------------------------------------------------------------
# Generating a trait data matrix
#-------------------------------------------------------------------------------

Tr <- cbind(1, replicate(N_traits - 1, rnorm(1:N_species)))
rownames(Tr) <- paste0("Sp.", 1:N_species)

#-------------------------------------------------------------------------------
# Generating a phylogenetic correlation matrix
#-------------------------------------------------------------------------------

tree <- rtree(N_species, tip.label = 1:N_species)
plot(tree)
C <- vcv(tree, corr = T)[order(as.numeric(tree$tip.label)), order(as.numeric(tree$tip.label))]
colnames(C) <- paste0("Sp.", 1:N_species)
rownames(C) <- paste0("Sp.", 1:N_species)

#-------------------------------------------------------------------------------
# Generate spatial data in the two-dimensional [-1,1] space
#-------------------------------------------------------------------------------

xy <- replicate(2, runif(N_sites))

#-------------------------------------------------------------------------------
# Preparing the spatial basis functions
#-------------------------------------------------------------------------------

N_spatial_bf_1D <- 10
spatial_bf <- cbind(bs(xy[,1], knots = seq(-1, 1, length.out = N_spatial_bf_1D - 2))[,-c(N_spatial_bf_1D + 1)],
                    bs(xy[,2], knots = seq(-1, 1, length.out = N_spatial_bf_1D - 2))[,-c(N_spatial_bf_1D + 1)])
spatial_bf_indices <- expand.grid(X = 1:N_spatial_bf_1D,
                                  Y = (N_spatial_bf_1D+1):(2*N_spatial_bf_1D))
spatial_bf_range <- expand.grid(x = seq(-1, 1, length.out = N_spatial_bf_1D),
                                y = seq(-1, 1, length.out = N_spatial_bf_1D))[which(colSums(spatial_bf[,spatial_bf_indices$X] * spatial_bf[,spatial_bf_indices$Y]) > 0),]
spatial_bf_indices_selected <- spatial_bf_indices[which(colSums(spatial_bf[,spatial_bf_indices$X] * spatial_bf[,spatial_bf_indices$Y]) > 0),]

#-------------------------------------------------------------------------------
# Preparing the spatiotemporal basis functions
#-------------------------------------------------------------------------------

N_spatiotemporal_bf_1D_space <- 4
N_spatiotemporal_bf_1D_time <- 6
spatiotemporal_bf_spatial <- cbind(bs(xy[,1], knots = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_space - 2))[,-c(N_spatiotemporal_bf_1D_space + 1)],
                                   bs(xy[,2], knots = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_space - 2))[,-c(N_spatiotemporal_bf_1D_space + 1)])
spatiotemporal_bf_temporal <- bs(1:N_years, knots = seq(1, N_years, length.out = N_spatiotemporal_bf_1D_time - 2))[,-c(N_spatiotemporal_bf_1D_time + 1)]
spatiotemporal_bf_indices <- expand.grid(X = (0*N_spatiotemporal_bf_1D_space+1):(1*N_spatiotemporal_bf_1D_space),
                                         Y = (1*N_spatiotemporal_bf_1D_space+1):(2*N_spatiotemporal_bf_1D_space),
                                         T = 1:N_spatiotemporal_bf_1D_time)
spatiotemporal_bf_range <- expand.grid(x = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_space),
                                       y = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_space),
                                       t = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_time))[which(colSums(spatiotemporal_bf_spatial[,spatiotemporal_bf_indices$X] * spatiotemporal_bf_spatial[,spatiotemporal_bf_indices$Y]) > 0),]
spatiotemporal_bf_indices_selected <- spatiotemporal_bf_indices[which(colSums(spatiotemporal_bf_spatial[,spatiotemporal_bf_indices$X] * spatiotemporal_bf_spatial[,spatiotemporal_bf_indices$Y]) > 0),]

#-------------------------------------------------------------------------------
# Compiling all data for the probabilistic model
#-------------------------------------------------------------------------------

datalist_simulation <- list(N_species = N_species,
                            N_sites = N_sites,
                            N_years = N_years,
                            site = as.numeric(visitdata_template$site),
                            year = as.numeric(visitdata_template$year),
                            N_allvisits = N_allvisits,
                            N_weeks = 53,
                            week = visitdata_template$week,
                            N_phenology_bf = ncol(phenology_bf),
                            phenology_bf = phenology_bf,
                            N_observers = length(levels(visitdata_template$observer)),
                            observer = as.numeric(visitdata_template$observer),
                            N_spatial_dims = 10,
                            N_spatiotemporal_dims = 5,
                            N_observer_dims = 5,
                            N_spatial_bf = nrow(spatial_bf_indices_selected),
                            N_spatial_bf_uni = ncol(spatial_bf),
                            spatial_bf_indices = as.matrix(spatial_bf_indices_selected),
                            spatial_bf = spatial_bf,
                            spatial_bf_range = as.matrix(spatial_bf_range),
                            N_spatiotemporal_bf = nrow(spatiotemporal_bf_range),
                            N_spatiotemporal_bf_uni_spatial = ncol(spatiotemporal_bf_spatial),
                            N_spatiotemporal_bf_uni_temporal = ncol(spatiotemporal_bf_temporal),
                            spatiotemporal_bf_indices = as.matrix(spatiotemporal_bf_indices_selected),
                            spatiotemporal_bf_spatial = spatiotemporal_bf_spatial,
                            spatiotemporal_bf_temporal = spatiotemporal_bf_temporal,
                            spatiotemporal_bf_range = as.matrix(spatiotemporal_bf_range),
                            N_occupancy_covs = ncol(occupancy_X),
                            occupancy_X = occupancy_X,
                            N_detection_covs = ncol(detection_X),
                            detection_X = detection_X,
                            N_traits = ncol(Tr),
                            Tr = t(Tr),
                            C = C,
                            year_range = seq(-1, 1, length.out = N_years))

#-------------------------------------------------------------------------------
# Sampling simulated data from the prior distribution
#-------------------------------------------------------------------------------

simulator <- cmdstan_model("F1_Demo/Data/opportunisticJSDM_simulator.stan")
stan_realization <- simulator$sample(data = datalist_simulation, chains = 1, iter_sampling = 1, fixed_param = T)
realization <- read_cmdstan_csv(stan_realization$output_files(), format = "draws_matrix")$post_warmup_draws

#-------------------------------------------------------------------------------
# Derive true occupancies and detections
#-------------------------------------------------------------------------------

Z <- data.frame(expand.grid(site = 1:N_sites,
                            year = 1:N_years,
                            species = 1:N_species),
                Z = c(realization[1,substr(colnames(realization), 1, 2) == "Z["])) %>%
  pivot_wider(names_from = species, values_from = Z)
image(t(as.matrix(select(Z, -site, -year))))
colMeans(Z)

Y <- data.frame(expand.grid(visitid = 1:N_allvisits,
                            species = 1:N_species),
                Y = c(realization[1,substr(colnames(realization), 1, 2) == "Y["])) %>%
  pivot_wider(names_from = species, values_from = Y) %>%
  select(-visitid) %>%
  as.matrix()
colnames(Y) <- paste0("Sp.", 1:N_species)
image(t(Y))
colMeans(Y)

#-------------------------------------------------------------------------------
# Derive visit and community data from the sightings 
#-------------------------------------------------------------------------------

visitdata <- data.frame(expand.grid(visitid = 1:N_allvisits,
                                    species = 1:N_species),
                        Y = c(realization[1,substr(colnames(realization), 1, 2) == "Y["])) %>%
  pivot_wider(names_from = species, values_from = Y) %>%
  left_join(visitdata_template, .)

#-------------------------------------------------------------------------------
# Establish confirmed presences 
#-------------------------------------------------------------------------------

any_seen <- visitdata %>%
  select(-visitid, -week, -observer) %>%
  pivot_longer(!c(site, year), names_to = "species") %>%
  group_by(site, year, species) %>%
  summarize(value = 1*(sum(value) > 0)) %>%
  ungroup() %>%
  mutate(species = factor(species, 1:N_species)) %>%
  complete(site, year, species, fill = list(value = 0)) %>%
  pivot_wider(names_from = species, names_sort = T, values_from = value)

#-------------------------------------------------------------------------------
# Compute visit metadata
#-------------------------------------------------------------------------------

# For each site and year, the first and last visit id (i.e. row number of the
# visitdata dataset) as well as the total number of visits is computed. This
# data is used for bookkeeping purposes in the model.
visitmetadata <- visitdata %>%
  group_by(site, year) %>%
  summarise(start_visitid = first(visitid),
            end_visitid = last(visitid),
            N_visits = n()) %>%
  ungroup() %>%
  complete(site, year, fill = list(start_visitid = -1, end_visitid = -1, N_visits = 0))

#-------------------------------------------------------------------------------
# Export datasets
#-------------------------------------------------------------------------------

write.table(data.frame(visitid = 1:N_allvisits, Y), "F1_Demo/Data/detection_data.csv", row.names = F, col.names = T, sep = ",")
write.table(data.frame(site = 1:N_sites, xy), "F1_Demo/Data/spatial_data.csv", row.names = F, col.names = T, sep = ",")
write.table(visitdata_template, "F1_Demo/Data/visitmetadata.csv", row.names = F, col.names = T, sep = ",")
write.table(data.frame(site = 1:N_sites, occupancy_X), "F1_Demo/Data/occupancy_predictors.csv", row.names = F, col.names = T, sep = ",")
write.table(data.frame(visitid = 1:N_allvisits, detection_X), "F1_Demo/Data/detection_predictors.csv", row.names = F, col.names = T, sep = ",")
write.table(C, "F1_Demo/Data/phylocor.csv", row.names = F, col.names = T, sep = ",")
write.table(Tr, "F1_Demo/Data/trait_data.csv", row.names = F, col.names = T, sep = ",")

