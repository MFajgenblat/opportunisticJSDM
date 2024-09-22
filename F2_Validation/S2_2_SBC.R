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

#-------------------------------------------------------------------------------
# Generating a phylogenetic correlation matrix
#-------------------------------------------------------------------------------

tree <- rtree(N_species, tip.label = 1:N_species)
plot(tree)
C <- vcv(tree, corr = T)[order(as.numeric(tree$tip.label)), order(as.numeric(tree$tip.label))]

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

#set_cmdstan_path("C:/Users/maxim/Documents/cmdstan-2.30.1")
simulator <- cmdstan_model("F2_Validation/opportunisticJSDM_simulator.stan")
stan_realization <- simulator$sample(data = datalist_simulation, chains = 1, iter_sampling = 1, fixed_param = T)
realization <- read_cmdstan_csv(stan_realization$output_files(), format = "draws_matrix")$post_warmup_draws

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
# Compiling all data for the probabilistic model
#-------------------------------------------------------------------------------

datalist <- list(N_species = ncol(select(any_seen, -site, -year)),
                 N_sites = length(levels(visitmetadata$site)),
                 N_years = length(levels(visitmetadata$year)),
                 N_visited_units = sum(visitmetadata$N_visits > 0),
                 visited_units = which(visitmetadata$N_visits > 0),
                 site = as.numeric(visitmetadata$site),
                 year = as.numeric(visitmetadata$year),
                 N_visits = visitmetadata$N_visits,
                 any_seen = as.matrix(select(any_seen, -site, -year)),
                 start_visitid = visitmetadata$start_visitid,
                 end_visitid = visitmetadata$end_visitid,
                 N_allvisits = sum(visitmetadata$N_visits),
                 Y = as.matrix(select(visitdata, -visitid, -site, -year, -week, -observer)),
                 N_weeks = 53,
                 week = visitdata$week,
                 N_phenology_bf = ncol(phenology_bf),
                 phenology_bf = phenology_bf,
                 N_observers = length(levels(visitdata$observer)),
                 observer = as.numeric(visitdata$observer),
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
                 year_range = seq(-1, 1, length.out = length(levels(visitmetadata$year))))

#-------------------------------------------------------------------------------
# MCMC sampling
#-------------------------------------------------------------------------------

model <- cmdstan_model("F2_Validation/opportunisticJSDM_analyzer.stan", cpp_options = list(stan_threads = T))
stanfit <- model$sample(data = datalist, chains = 1, iter_warmup = 100, iter_sampling = 100, refresh = 1, parallel_chains = 1, threads_per_chain = 16)
fit1 <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws

#-------------------------------------------------------------------------------
# Randomly sample or read existing focal parameters
#-------------------------------------------------------------------------------

if (exists("F2_Validation/focal_parameters.rds")) {
  focal_parameters <- readRDS("F2_Validation/focal_parameters.rds")
} else {
  focal_parameters <- unname(sapply(c("occupancy_G\\[", "occupancy_phylo_fraction", "occupancy_cor_covariates_L\\[", "occupancy_cov_covariates_scale\\[", "occupancy_beta_std\\[", "detection_G\\[", "detection_phylo_fraction", "detection_cor_covariates_L\\[", "detection_beta_std\\[", "temporal_trend_unstr_sd\\[", "temporal_trend_unstr_std\\[", "temporal_trend_str_ls\\[", "temporal_trend_str_sd\\[", "temporal_trend_str_std\\[", "phenology_overall_sd", "phenology_overall_std\\[", "phenology_sd", "phenology_std\\[", "spatial_Lambda_std\\[", "spatial_Lambda_delta\\[", "spatial_structured_fraction\\[", "spatial_Eta_str_std\\[", "spatial_Eta_str_ls\\[", "spatial_Eta_unstr\\[", "spatiotemporal_Lambda_std\\[", "spatiotemporal_Lambda_delta\\[", "spatiotemporal_Eta_str_std\\[", "spatiotemporal_Eta_str_ls_temporal\\[", "spatiotemporal_Eta_str_ls_spatial\\[", "observer_Lambda_std\\[", "observer_Lambda_delta\\[", "observer_Eta\\["),
                                    function(i) sample(colnames(fit1)[grep(i, colnames(fit1))], 1, replace = T)))
  saveRDS(focal_parameters, "F2_Validation/focal_parameters.rds")
}

#-------------------------------------------------------------------------------
# Write outcomes of simulation
#-------------------------------------------------------------------------------

write.table(data.frame(Parameter = focal_parameters,
                       Rank = unname(sapply(focal_parameters, function(i) mean(c(fit1[,i]) <= c(realization[,i])))),
                       timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")),
            file = paste0("F2_Validation/SBC_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
            sep = ",", row.names = FALSE, col.names = TRUE)

#-------------------------------------------------------------------------------
# Compile all individual simulations
#-------------------------------------------------------------------------------

SBC_data <- bind_rows(lapply(list.files(path = "F2_Validation/SBC", pattern = "^SBC_*", full.names = TRUE), read.csv))

#-------------------------------------------------------------------------------
# Plot histograms
#-------------------------------------------------------------------------------

N_bins <- 10
SBC_data %>%
  mutate(bin = cut(Rank, breaks = seq(0, 1, length.out = N_bins + 1), include.lowest = TRUE)) %>%
  group_by(Parameter, bin) %>%
  summarise(count = n(), .groups = 'drop_last') %>%
  summarise(chi_square_stat = sum((count - mean(count))^2 / mean(count)), 
            df = N_bins - 1, 
            chi_square_p_value = pchisq(chi_square_stat, df, lower.tail = FALSE)) %>%
  mutate(chi_square_p_value = sprintf("<br><i>p</i> = %.3f", chi_square_p_value)) %>%
  left_join(SBC_data, .) %>%
  ggplot() +
  geom_histogram(aes(x = Rank*100), bins = N_bins + 1, boundary = 0, closed = "right", color = "black", fill = "#a7c2c4", linewidth = 0.1) +
  geom_hline(yintercept = length(unique(SBC_data$timestamp))/(N_bins+1), linetype = "dashed", color = "#942b4a") +
  geom_hline(yintercept = 0, color = "black") +
  scale_x_continuous("Rank statistic", expand = c(0,0), breaks = seq(0, 100, by = 20)) +
  scale_y_continuous("Frequency", expand = c(0,0)) +
  facet_wrap(~ paste0("<b>", Parameter, "</b>", chi_square_p_value), ncol = 4) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.text = element_markdown(size = 6),
        strip.clip = "off",
        axis.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("SBC.png", width = 16, height = 20, units = "cm", dpi = 600)
