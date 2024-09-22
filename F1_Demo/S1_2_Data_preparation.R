#-------------------------------------------------------------------------------
# Setting working directory
#-------------------------------------------------------------------------------

# Choose the top-level directory
setwd("~/opportunisticJSDM")

#-------------------------------------------------------------------------------
# Loading packages
#-------------------------------------------------------------------------------

library(tidyverse)
library(mgcv)
library(splines)

#-------------------------------------------------------------------------------
# Reading datasets
#-------------------------------------------------------------------------------

detection_data <- read.csv("F1_Demo/Data/detection_data.csv")
spatial_data <- read.csv("F1_Demo/Data/spatial_data.csv")
visitmetadata <- read.csv("F1_Demo/Data/visitmetadata.csv")
occupancy_predictors <- read.csv("F1_Demo/Data/occupancy_predictors.csv")
detection_predictors <- read.csv("F1_Demo/Data/detection_predictors.csv")
trait_data <- read.csv("F1_Demo/Data/trait_data.csv")
phylocor <- read.csv("F1_Demo/Data/phylocor.csv")

#-------------------------------------------------------------------------------
# Derive visit and community data from the sightings 
#-------------------------------------------------------------------------------

visitdata <- detection_data %>%
  left_join(visitmetadata, .)

#-------------------------------------------------------------------------------
# Establish confirmed presences 
#-------------------------------------------------------------------------------

any_seen <- visitdata %>%
  select(-visitid, -week, -observer) %>%
  pivot_longer(!c(site, year), names_to = "species") %>%
  group_by(site, year, species) %>%
  summarize(value = 1*(sum(value) > 0)) %>%
  ungroup() %>%
  mutate(species = factor(species, paste0("Sp.", 1:N_species))) %>%
  mutate(site = factor(site, 1:nrow(spatial_data)),
         year = factor(year)) %>%
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
  mutate(site = factor(site, 1:nrow(spatial_data)),
         year = factor(year)) %>%
  complete(site, year, fill = list(start_visitid = 0, end_visitid = 0, N_visits = 0))

#-------------------------------------------------------------------------------
# Compute cyclic spline basis functions for weekly phenology estimation
#-------------------------------------------------------------------------------

phenology_bf <- cSplineDes(1:53, seq(1, 53, length.out = 11))

#-------------------------------------------------------------------------------
# Compiling the environmental data matrix
#-------------------------------------------------------------------------------

# The log-Area of waterbodies is rescaled around the mean log-area, with a range
# of 1, to feature a similar scale as the land cover data.
# The agricultural land cover class is omitted to avoid overspecification (since
# land cover data sums to 1). This specific class is chosen as reference class
# (its effect is absorbed by the intercept) since it is the dominant land cover
# class in Flanders.
occupancy_X <- occupancy_predictors %>%
  select(-site) %>%
  as.matrix()

#-------------------------------------------------------------------------------
# Compiling the detection covariate matrix
#-------------------------------------------------------------------------------

# The continuous list length is discretized into three categories: visits with
# only a single species reported, visits with 2 to 4 reported species and visits
# with 5 or more reported species. The latter is chosen as reference type (its
# effect is absorbed by the intercept), and the two other types are represented
# by means of two dummy variables.
detection_X <- detection_predictors %>%
  select(-visitid) %>%
  as.matrix()

#-------------------------------------------------------------------------------
# Compiling the trait data matrix
#-------------------------------------------------------------------------------

Tr <- trait_data %>%
  as.matrix()

#-------------------------------------------------------------------------------
# Computing the phylogenetic correlation matrix
#-------------------------------------------------------------------------------

C <- phylocor %>%
  as.matrix()

#-------------------------------------------------------------------------------
# Scale spatial data between [-1,1]
#-------------------------------------------------------------------------------

# Spatial coordinates are rescaled to fit within a square bounding box ranging
# from -1 and 1 along each dimension.
xy <- spatial_data[,c("X1","X2")] %>%
  rename(X = X1, Y = X2) %>%
  data.frame() %>%
  mutate(X = X - min(X),
         Y = Y - min(Y),
         X_rescaled = 2 * (X / max(c(X,Y))),
         Y_rescaled = 2 * (Y / max(c(X,Y))),
         X_rescaled = X_rescaled - max(X_rescaled)/2,
         Y_rescaled = Y_rescaled - max(Y_rescaled)/2) %>%
  dplyr::select(X_rescaled, Y_rescaled)

#-------------------------------------------------------------------------------
# Preparing the spatial basis functions
#-------------------------------------------------------------------------------

# The spatial random effects approach in the model relies on B-splines Projected
# Gaussian Process. For this, two-dimensional B-spline basis functions are set
# up and provided as data to the model. A memory-efficient representation is
# chosen, with the tensor product only being computed in a parallelized fashion
# within the model (likelihood).
N_spatial_bf_1D <- 15
spatial_bf <- cbind(bs(xy[,1], knots = seq(-1, 1, length.out = N_spatial_bf_1D - 2))[,-c(N_spatial_bf_1D + 1)],
                    bs(xy[,2], knots = seq(-1, 1, length.out = N_spatial_bf_1D - 2))[,-c(N_spatial_bf_1D + 1)])
spatial_bf_indices <- expand.grid(X = 1:N_spatial_bf_1D,
                                  Y = (N_spatial_bf_1D+1):(2*N_spatial_bf_1D))
spatial_bf_range <- expand.grid(x = seq(-1, 1, length.out = N_spatial_bf_1D),
                                y = seq(-1, 1, length.out = N_spatial_bf_1D))[which(colSums(spatial_bf[,spatial_bf_indices$X] * spatial_bf[,spatial_bf_indices$Y]) > 0),]
spatial_bf_indices_selected <- spatial_bf_indices[which(colSums(spatial_bf[,spatial_bf_indices$X] * spatial_bf[,spatial_bf_indices$Y]) > 0),]

# To get a taste of the maximal spatial resolution that can be achieved, the
# following ggplot command can be used to generate randomly drawn spatial
# patterns at the highest available resolution.
#ggplot(data.frame(xy, value = (spatial_bf[,spatial_bf_indices_selected$X] * spatial_bf[,spatial_bf_indices_selected$Y]) %*% rnorm(nrow(spatial_bf_indices_selected)))) + geom_point(aes(x = X_rescaled, y = Y_rescaled, color = value)) + scale_color_viridis_c() + coord_equal()

#-------------------------------------------------------------------------------
# Preparing the spatiotemporal basis functions
#-------------------------------------------------------------------------------

# In an analogous fashion, spatiotemporal basis functions are computed.
N_spatiotemporal_bf_1D_space <- 9
N_spatiotemporal_bf_1D_time <- 6
spatiotemporal_bf_spatial <- cbind(bs(xy[,1], knots = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_space - 2))[,-c(N_spatiotemporal_bf_1D_space + 1)],
                                   bs(xy[,2], knots = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_space - 2))[,-c(N_spatiotemporal_bf_1D_space + 1)])
spatiotemporal_bf_temporal <- bs(1:length(levels(visitmetadata$year)), knots = seq(1, length(levels(visitmetadata$year)), length.out = N_spatiotemporal_bf_1D_time - 2))[,-c(N_spatiotemporal_bf_1D_time + 1)]
spatiotemporal_bf_indices <- expand.grid(X = (0*N_spatiotemporal_bf_1D_space+1):(1*N_spatiotemporal_bf_1D_space),
                                         Y = (1*N_spatiotemporal_bf_1D_space+1):(2*N_spatiotemporal_bf_1D_space),
                                         T = 1:N_spatiotemporal_bf_1D_time)
spatiotemporal_bf_range <- expand.grid(x = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_space),
                                       y = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_space),
                                       t = seq(-1, 1, length.out = N_spatiotemporal_bf_1D_time))[which(colSums(spatiotemporal_bf_spatial[,spatiotemporal_bf_indices$X] * spatiotemporal_bf_spatial[,spatiotemporal_bf_indices$Y]) > 0),]
spatiotemporal_bf_indices_selected <- spatiotemporal_bf_indices[which(colSums(spatiotemporal_bf_spatial[,spatiotemporal_bf_indices$X] * spatiotemporal_bf_spatial[,spatiotemporal_bf_indices$Y]) > 0),]

# Optionally, random realizations of the spatiotemporal splines can be generated
# For computational ease, each year is computed separately
#temp <- rnorm(sum(spatiotemporal_bf_indices[,3] == 1))
#ggplot(data.frame(xy, value = (spatiotemporal_bf_spatial[,spatiotemporal_bf_indices_selected$X[spatiotemporal_bf_indices[,3] == 1]] * spatiotemporal_bf_spatial[,spatiotemporal_bf_indices_selected$Y[spatiotemporal_bf_indices[,3] == 1]] * spatiotemporal_bf_temporal[rep(1,nrow(spatiotemporal_bf_spatial)),spatiotemporal_bf_indices_selected$T[spatiotemporal_bf_indices[,3] == 1]]) %*% temp)) + geom_point(aes(x = X_rescaled, y = Y_rescaled, color = value)) + scale_color_viridis_c() + coord_equal()
#ggplot(data.frame(xy, value = (spatiotemporal_bf_spatial[,spatiotemporal_bf_indices_selected$X[spatiotemporal_bf_indices[,3] == 1]] * spatiotemporal_bf_spatial[,spatiotemporal_bf_indices_selected$Y[spatiotemporal_bf_indices[,3] == 1]] * spatiotemporal_bf_temporal[rep(2,nrow(spatiotemporal_bf_spatial)),spatiotemporal_bf_indices_selected$T[spatiotemporal_bf_indices[,3] == 1]]) %*% temp)) + geom_point(aes(x = X_rescaled, y = Y_rescaled, color = value)) + scale_color_viridis_c() + coord_equal()
#ggplot(data.frame(xy, value = (spatiotemporal_bf_spatial[,spatiotemporal_bf_indices_selected$X[spatiotemporal_bf_indices[,3] == 1]] * spatiotemporal_bf_spatial[,spatiotemporal_bf_indices_selected$Y[spatiotemporal_bf_indices[,3] == 1]] * spatiotemporal_bf_temporal[rep(3,nrow(spatiotemporal_bf_spatial)),spatiotemporal_bf_indices_selected$T[spatiotemporal_bf_indices[,3] == 1]]) %*% temp)) + geom_point(aes(x = X_rescaled, y = Y_rescaled, color = value)) + scale_color_viridis_c() + coord_equal()
# ...

#-------------------------------------------------------------------------------
# Compiling all data for the probabilistic model
#-------------------------------------------------------------------------------

datalist <- list(N_species = ncol(select(any_seen, starts_with("Sp."))),
                 N_sites = length(unique(visitmetadata$site)),
                 N_years = length(unique(visitmetadata$year)),
                 N_visited_units = sum(visitmetadata$N_visits > 0),
                 visited_units = which(visitmetadata$N_visits > 0),
                 site = as.numeric(visitmetadata$site),
                 year = as.numeric(visitmetadata$year),
                 N_visits = visitmetadata$N_visits,
                 any_seen = as.matrix(select(any_seen, starts_with("Sp."))),
                 start_visitid = visitmetadata$start_visitid,
                 end_visitid = visitmetadata$end_visitid,
                 N_allvisits = sum(visitmetadata$N_visits),
                 Y = as.matrix(select(visitdata, starts_with("Sp."))),
                 N_weeks = 53,
                 week = visitdata$week,
                 N_phenology_bf = ncol(phenology_bf),
                 phenology_bf = phenology_bf,
                 N_observers = length(unique(visitdata$observer)),
                 observer = as.numeric(factor(visitdata$observer)),
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
                 year_range = seq(-1, 1, length.out = length(unique(any_seen$year))))
saveRDS(datalist, "F1_Demo/Data/datalist.rds")
