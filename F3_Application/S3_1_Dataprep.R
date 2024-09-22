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

#-------------------------------------------------------------------------------
# Reading datasets
#-------------------------------------------------------------------------------

# Waterbody area and coordinates (Lambert 72 projection, epsg:31370)
waterbodies_metadata <- read.csv("F3_Application/Data/waterbodies_metadata.csv")
# Land cover in a 200m buffer surrounding each waterbody
waterbodies_landcover <- read.csv("F3_Application/Data/waterbodies_landcover.csv")
# Taxonomy and trait data for all species
species_metadata <- read.csv("F3_Application/Data/species_metadata.csv", sep=";")
# Reading a dataset with sightings for each waterbody and temporal metadata
# The public repository only contains a desensitized one
if (file.exists("F3_Application/Data/linked_sightings.csv")) {
  sightings <- read.csv("F3_Application/Data/linked_sightings.csv")
} else {
  sightings <- read.csv("F3_Application/Data/linked_sightings_desensitized.csv")
}

#-------------------------------------------------------------------------------
# Derive visit and community data from the sightings 
#-------------------------------------------------------------------------------

# By pivoting the sightings data to a wide format, a new dataset is obtained
# containing whether or not each species has been sighted during a single
# site visit by a specific observer on a particular day.
# Only observers having performed at least 75 site visits are retained
# The list length is the number of species recorded during each site visit
visitdata <- sightings %>%
  pivot_wider(id_cols = c(waterbody, year, week, scrambled_user_id),
              names_from = Species, names_sort = T, values_from = Species,
              values_fn = list(Species = function(i) as.numeric(length(i) > 0)), values_ = 0) %>%
  arrange(waterbody, year, week, scrambled_user_id) %>%
  filter(scrambled_user_id %in% as.numeric(names(which(table(pull(., scrambled_user_id)) >= 75)))) %>%
  filter(waterbody %in% as.numeric(names(which(table(pull(., waterbody)) > 2)))) %>%
  rowid_to_column(var = "visitid") %>%
  mutate(scrambled_user_id = factor(scrambled_user_id),
         listlength = rowSums(select(., starts_with("Sp"))), .after = scrambled_user_id)

# Count the number of retained observers and waterbodies
length(unique(visitdata$scrambled_user_id))
length(unique(visitdata$waterbody))

#-------------------------------------------------------------------------------
# Establish confirmed presences 
#-------------------------------------------------------------------------------

# If the presence of a species has been confirmed for a site and year, the
# likelihood is simplified (since the scenario that the species is absent can be
# ruled out). Therefore, a dataset is created containing information on the
# confirmed presence of species for each site and year.
any_seen <- sightings %>%
  filter(scrambled_user_id %in% visitdata$scrambled_user_id,
         paste0(waterbody, "_",year) %in% paste0(visitdata$waterbody, "_", visitdata$year)) %>%
  pivot_wider(id_cols = c(waterbody, year),
              names_from = Species, names_sort = T, names_expand = T, values_from = Species,
              values_fn = list(Species = function(i) as.numeric(length(i) > 0)), values_ = 0) %>%
  arrange(waterbody, year) %>%
  mutate(waterbody = factor(waterbody, levels = waterbodies_metadata$waterbody),
         year = factor(year))

#-------------------------------------------------------------------------------
# Compute visit metadata
#-------------------------------------------------------------------------------

# For each site and year, the first and last visit id (i.e. row number of the
# visitdata dataset) as well as the total number of visits is computed. This
# data is used for bookkeeping purposes in the model.
visitmetadata <- visitdata %>%
  group_by(waterbody, year) %>%
  summarise(start_visitid = first(visitid),
            end_visitdid = last(visitid),
            N_visits = n()) %>%
  ungroup() %>%
  mutate(waterbody = factor(waterbody, waterbodies_metadata$waterbody),
         year = factor(year))

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
occupancy_X <- data.frame(Intercept = 1,
                          Area = (log(waterbodies_metadata$Area) - mean(log(waterbodies_metadata$Area))) / (max(log(waterbodies_metadata$Area)) - min(log(waterbodies_metadata$Area))),
                          select(waterbodies_landcover, -waterbody, -Agricultural.land)) %>%
  as.matrix()

#-------------------------------------------------------------------------------
# Compiling the detection covariate matrix
#-------------------------------------------------------------------------------

# The continuous list length is discretized into three categories: visits with
# only a single species reported, visits with 2 to 4 reported species and visits
# with 5 or more reported species. The latter is chosen as reference type (its
# effect is absorbed by the intercept), and the two other types are represented
# by means of two dummy variables.
detection_X <- data.frame(Intercept = 1,
                          Singleton_list = as.numeric(visitdata$listlength == 1),
                          Short_list = as.numeric(visitdata$listlength %in% 2:4)) %>%
  as.matrix()

#-------------------------------------------------------------------------------
# Compiling the trait data matrix
#-------------------------------------------------------------------------------

# To avoid confounding between the suborder (dragonflies vs damselflies) and
# body size, body size is subtracted with the mean suborder body size.
species_metadata <- species_metadata %>%
  group_by(Suborder) %>%
  mutate(Body_size_scaled = Body_size - mean(Body_size))

Tr <- data.frame(Intercept = 1,
                 scale(species_metadata[,c("Body_size_scaled", "Migratory", "Permanent_ponds", "Nutrient_level", "Antipredation_index", "STI")])) %>%
  as.matrix()

#-------------------------------------------------------------------------------
# Computing the phylogenetic correlation matrix
#-------------------------------------------------------------------------------

# An approximate phylogenetic correlation matrix is constructed using the
# taxonomic classification.
tree <- species_metadata %>%
  mutate(Species = factor(Species),
         Genus = factor(Genus),
         Family = factor(Family),
         Superfamily = factor(Superfamily),
         Suborder = factor(Suborder),
         Order = factor(Order)) %>%
  as.phylo(~ Order/Suborder/Superfamily/Family/Genus/Species, data = ., collapse = FALSE)
tree$edge.length = rep(1, length(tree$edge))
C <- vcv(tree, corr = T)[species_metadata$Species, species_metadata$Species]
#plot(tree)
#image(C)

#-------------------------------------------------------------------------------
# Scale spatial data between [-1,1]
#-------------------------------------------------------------------------------

# Spatial coordinates are rescaled to fit within a square bounding box ranging
# from -1 and 1 along each dimension.
xy <- waterbodies_metadata[,c("X","Y")] %>%
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
N_spatial_bf_1D <- 40
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
                 N_sites = length(levels(droplevels(visitmetadata$waterbody))),
                 N_years = length(levels(visitmetadata$year)),
                 N_visited_units = nrow(visitmetadata),
                 visited_units = 1:nrow(visitmetadata),
                 site = as.numeric(droplevels(visitmetadata$waterbody)),
                 year = as.numeric(visitmetadata$year),
                 N_visits = visitmetadata$N_visits,
                 any_seen = as.matrix(select(any_seen, starts_with("Sp."))),
                 start_visitid = visitmetadata$start_visitid,
                 end_visitid = visitmetadata$end_visitdid,
                 N_allvisits = sum(visitmetadata$N_visits),
                 Y = as.matrix(select(visitdata, starts_with("Sp."))),
                 N_weeks = 53,
                 week = visitdata$week,
                 N_phenology_bf = ncol(phenology_bf),
                 phenology_bf = phenology_bf,
                 N_observers = length(levels(visitdata$scrambled_user_id)),
                 observer = as.numeric(visitdata$scrambled_user_id),
                 N_spatial_dims = 10,
                 N_spatiotemporal_dims = 5,
                 N_observer_dims = 5,
                 N_spatial_bf = nrow(spatial_bf_indices_selected),
                 N_spatial_bf_uni = ncol(spatial_bf),
                 spatial_bf_indices = as.matrix(spatial_bf_indices_selected),
                 spatial_bf = spatial_bf[as.numeric(levels(droplevels(visitmetadata$waterbody))),],
                 spatial_bf_range = as.matrix(spatial_bf_range),
                 N_spatiotemporal_bf = nrow(spatiotemporal_bf_range),
                 N_spatiotemporal_bf_uni_spatial = ncol(spatiotemporal_bf_spatial),
                 N_spatiotemporal_bf_uni_temporal = ncol(spatiotemporal_bf_temporal),
                 spatiotemporal_bf_indices = as.matrix(spatiotemporal_bf_indices_selected),
                 spatiotemporal_bf_spatial = spatiotemporal_bf_spatial[as.numeric(levels(droplevels(visitmetadata$waterbody))),],
                 spatiotemporal_bf_temporal = spatiotemporal_bf_temporal,
                 spatiotemporal_bf_range = as.matrix(spatiotemporal_bf_range),
                 N_occupancy_covs = ncol(occupancy_X),
                 occupancy_X = occupancy_X[as.numeric(levels(droplevels(visitmetadata$waterbody))),],
                 N_detection_covs = ncol(detection_X),
                 detection_X = detection_X,
                 N_traits = ncol(Tr),
                 Tr = t(Tr),
                 C = C,
                 year_range = seq(-1, 1, length.out = length(levels(visitmetadata$year))))
saveRDS(datalist, "F3_Application/Data/datalist.rds")
