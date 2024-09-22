#-------------------------------------------------------------------------------
# Setting working directory
#-------------------------------------------------------------------------------

# Choose the top-level directory
setwd("~/opportunisticJSDM")

#-------------------------------------------------------------------------------
# Loading packages
#-------------------------------------------------------------------------------

library(cmdstanr)

#-------------------------------------------------------------------------------
# Reading data list
#-------------------------------------------------------------------------------

datalist <- readRDS("F1_Demo/Data/datalist.rds")

#-------------------------------------------------------------------------------
# MCMC sampling
#-------------------------------------------------------------------------------

model <- cmdstan_model("F1_Demo/opportunisticJSDM.stan", cpp_options = list(stan_threads = T))
stanfit <- model$sample(data = datalist, chains = 12, iter_warmup = 100, iter_sampling = 100, refresh = 1, parallel_chains = 12, adapt_delta = 0.8, threads_per_chain = 16)

#-------------------------------------------------------------------------------
# Extracting and saving MCMC draws
#-------------------------------------------------------------------------------

fit <- read_cmdstan_csv(stanfit$output_files(), format = "draws_matrix")$post_warmup_draws
saveRDS(fit, "F1_Demo/fit.rds")
