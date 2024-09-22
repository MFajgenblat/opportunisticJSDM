functions {
  
  vector GP_1D(real[] x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9))*eta);
  }
  
  vector GP_nD_isotropic(vector[] x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9))*eta);
  }
  
  vector GP_nD_anisotropic(vector[] x, real[] rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9))*eta);
  }

}

data {
  
  int<lower=0> N_species;
  int<lower=0> N_sites;
  int<lower=0> N_years;
  int<lower=0> N_allvisits;
  int site[N_allvisits];
  int year[N_allvisits];
  
  int N_weeks;
  int week[N_allvisits];
  int N_phenology_bf;
  matrix[N_weeks,N_phenology_bf] phenology_bf;
  
  int N_observers;
  int observer[N_allvisits];
  
  int N_observer_dims;
  int N_spatial_dims;
  int N_spatiotemporal_dims;
  
  int N_spatial_bf;
  int N_spatial_bf_uni;
  int spatial_bf_indices[N_spatial_bf,2];
  matrix[N_sites,N_spatial_bf_uni] spatial_bf;
  vector<lower=-1,upper=1>[2] spatial_bf_range[N_spatial_bf];
                 
  int N_spatiotemporal_bf;
  int N_spatiotemporal_bf_uni_spatial;
  int N_spatiotemporal_bf_uni_temporal;
  int spatiotemporal_bf_indices[N_spatiotemporal_bf,3];
  matrix[N_sites,N_spatiotemporal_bf_uni_spatial] spatiotemporal_bf_spatial;
  matrix[N_years,N_spatiotemporal_bf_uni_temporal] spatiotemporal_bf_temporal;
  vector<lower=-1,upper=1>[3] spatiotemporal_bf_range[N_spatiotemporal_bf];
  
  int<lower=0> N_occupancy_covs;
  matrix[N_sites,N_occupancy_covs] occupancy_X;
  int<lower=0> N_detection_covs;
  matrix[N_allvisits,N_detection_covs] detection_X;
  
  int<lower=0> N_traits;
  matrix[N_traits,N_species] Tr;
  
  corr_matrix[N_species] C;
  
  real year_range[N_years];

}

generated quantities {
  
  // Parameter definition
  
  matrix[N_occupancy_covs,N_traits] occupancy_G;
  real<lower=0,upper=1> occupancy_phylo_fraction;
  cholesky_factor_corr[N_occupancy_covs] occupancy_cor_covariates_L;
  vector<lower=0>[N_occupancy_covs] occupancy_cov_covariates_scale;
  matrix[N_occupancy_covs,N_species] occupancy_beta_std;
  
  matrix[N_detection_covs,N_traits] detection_G;
  real<lower=0,upper=1> detection_phylo_fraction;
  cholesky_factor_corr[N_detection_covs] detection_cor_covariates_L;
  vector<lower=0>[N_detection_covs] detection_cov_covariates_scale;
  matrix[N_detection_covs,N_species] detection_beta_std;
  
  real<lower=0> temporal_trend_unstr_sd[N_species];
  vector[N_years] temporal_trend_unstr_std[N_species];
  real<lower=0> temporal_trend_str_ls[N_species];
  real<lower=0> temporal_trend_str_sd[N_species];
  vector[N_years] temporal_trend_str_std[N_species];
  
  real<lower=0> phenology_overall_sd;
  vector[N_phenology_bf] phenology_overall_std;
  real<lower=0> phenology_sd[N_species];
  vector[N_phenology_bf] phenology_std[N_species];
  
  matrix[N_spatial_dims,N_species] spatial_Lambda_std;
  vector<lower=0>[N_spatial_dims] spatial_Lambda_delta;
  row_vector<lower=0,upper=1>[N_spatial_dims] spatial_structured_fraction;
  matrix[N_spatial_bf,N_spatial_dims] spatial_Eta_str_std;
  real<lower=0> spatial_Eta_str_ls[N_spatial_dims];
  matrix[N_sites,N_spatial_dims] spatial_Eta_unstr;
  
  matrix[N_spatiotemporal_dims,N_species] spatiotemporal_Lambda_std;
  vector<lower=0>[N_spatiotemporal_dims] spatiotemporal_Lambda_delta;
  matrix[N_spatiotemporal_bf,N_spatiotemporal_dims] spatiotemporal_Eta_str_std;
  real<lower=0> spatiotemporal_Eta_str_ls_temporal[N_spatiotemporal_dims];
  real<lower=0> spatiotemporal_Eta_str_ls_spatial[N_spatiotemporal_dims];
  
  matrix[N_observer_dims,N_species] observer_Lambda_std;
  vector<lower=0>[N_observer_dims] observer_Lambda_delta;
  matrix[N_observers,N_observer_dims] observer_Eta;
  
  // Transformed parameter definition
  
  matrix[N_occupancy_covs,N_species] occupancy_beta;
  cholesky_factor_cov[N_occupancy_covs] occupancy_cov_covariates_L;
  cholesky_factor_corr[N_species] occupancy_cor_species_L;
  
  matrix[N_detection_covs,N_species] detection_beta;
  cholesky_factor_cov[N_detection_covs] detection_cov_covariates_L;
  cholesky_factor_corr[N_species] detection_cor_species_L;
  
  matrix[N_years,N_species] temporal_trend;
  
  vector[N_weeks] phenology_overall;
  matrix[N_weeks,N_species] phenology;
  
  vector[N_spatial_dims] spatial_Lambda_sd;
  matrix[N_spatial_dims,N_species] spatial_Lambda;
  
  vector[N_spatiotemporal_dims] spatiotemporal_Lambda_sd;
  matrix[N_spatiotemporal_dims,N_species] spatiotemporal_Lambda;
  
  vector[N_observer_dims] observer_Lambda_sd;
  matrix[N_observer_dims,N_species] observer_Lambda;
  
  matrix[N_spatial_bf,N_spatial_dims] spatial_Eta_str;
  matrix[N_spatiotemporal_bf,N_spatiotemporal_dims] spatiotemporal_Eta_str;
  
  // Key outcome quantities definition
  
  int<lower=0,upper=1> Z[N_sites,N_years,N_species];
  int<lower=0,upper=1> Y[N_allvisits,N_species];
  
  // Sampling from the prior distributions
  
  for (i in 1:N_occupancy_covs) {
    for (j in 1:N_traits) {
      if (i == 1 && j == 1) {
        occupancy_G[i,j] = normal_rng(-15, 1);
      } else {
        occupancy_G[i,j] = normal_rng(0, 3);
      }
    }
    occupancy_cov_covariates_scale[i] = abs(normal_rng(0, 3));
    for (j in 1:N_species) {
      occupancy_beta_std[i,j] = std_normal_rng();
    }
  }
  occupancy_phylo_fraction = uniform_rng(0, 1);
  occupancy_cor_covariates_L = lkj_corr_cholesky_rng(N_occupancy_covs, 2.0);
  
  for (i in 1:N_detection_covs) {
    for (j in 1:N_traits) {
      if (i == 1 && j == 1) {
        detection_G[i,j] = normal_rng(-15, 1);
      } else {
        detection_G[i,j] = normal_rng(0, 3);
      }
    }
    detection_cov_covariates_scale[i] = abs(normal_rng(0, 3));
    for (j in 1:N_species) {
      detection_beta_std[i,j] = std_normal_rng();
    }
  }
  detection_phylo_fraction = uniform_rng(0, 1);
  detection_cor_covariates_L = lkj_corr_cholesky_rng(N_detection_covs, 2.0);
  
  for (i in 1:N_species) {
    temporal_trend_unstr_sd[i]= abs(normal_rng(0, 3));
    for (j in 1:N_years) {
      temporal_trend_unstr_std[i,j] = std_normal_rng();
      temporal_trend_str_std[i,j] = std_normal_rng();
    }
    temporal_trend_str_ls[i]= inv_gamma_rng(5, 5);
    temporal_trend_str_sd[i]= abs(normal_rng(0, 3));
  }
  
  phenology_overall_sd = abs(normal_rng(0, 3));
  for (i in 1:N_phenology_bf) {
    phenology_overall_std[i] = std_normal_rng();
  }
  
  for (i in 1:N_species) {
    phenology_sd[i] = abs(normal_rng(0, 3));
    for (j in 1:N_phenology_bf) {
      phenology_std[i,j] = std_normal_rng();
    }
  }
  
  for (i in 1:N_spatial_dims) {
    for (j in 1:N_species) {
      spatial_Lambda_std[i,j] = std_normal_rng();
    }
    if (i == 1) {spatial_Lambda_delta[i] = gamma_rng(2, 1);}
    if (i > 1) {spatial_Lambda_delta[i] = gamma_rng(6, 1);}
    spatial_structured_fraction[i] = uniform_rng(0, 1);
    for (j in 1:N_spatial_bf) {
      spatial_Eta_str_std[j,i] = std_normal_rng();
    }
    spatial_Eta_str_ls[i] = inv_gamma_rng(5, 5);
    for (j in 1:N_sites) {
      spatial_Eta_unstr[j,i] = std_normal_rng();
    }
  }
  
  for (i in 1:N_spatiotemporal_dims) {
    for (j in 1:N_species) {
      spatiotemporal_Lambda_std[i,j] = std_normal_rng();
    }
    if (i == 1) {spatiotemporal_Lambda_delta[i] = gamma_rng(2, 1);}
    if (i > 1) {spatiotemporal_Lambda_delta[i] = gamma_rng(6, 1);}
    for (j in 1:N_spatiotemporal_bf) {
      spatiotemporal_Eta_str_std[j,i] = std_normal_rng();
    }
    spatiotemporal_Eta_str_ls_spatial[i] = inv_gamma_rng(5, 5);
    spatiotemporal_Eta_str_ls_temporal[i] = inv_gamma_rng(5, 5);
  }
  
  for (i in 1:N_observer_dims) {
    for (j in 1:N_species) {
      observer_Lambda_std[i,j] = std_normal_rng();
    }
    if (i == 1) {observer_Lambda_delta[i] = gamma_rng(2, 1);}
    if (i > 1) {observer_Lambda_delta[i] = gamma_rng(6, 1);}
    for (j in 1:N_observers) {
      observer_Eta[j,i] = std_normal_rng();
    }
  }
  
  // Transformed parameter computation
  
  occupancy_cov_covariates_L = quad_form_diag(occupancy_cor_covariates_L, occupancy_cov_covariates_scale);
  occupancy_cor_species_L = cholesky_decompose(occupancy_phylo_fraction*C + (1-occupancy_phylo_fraction)*identity_matrix(N_species));
  occupancy_beta = occupancy_G*Tr + occupancy_cov_covariates_L*occupancy_beta_std*occupancy_cor_species_L';
  
  detection_cov_covariates_L = quad_form_diag(detection_cor_covariates_L, detection_cov_covariates_scale);
  detection_cor_species_L = cholesky_decompose(detection_phylo_fraction*C + (1-detection_phylo_fraction)*identity_matrix(N_species));
  detection_beta = detection_G*Tr + detection_cov_covariates_L*detection_beta_std*detection_cor_species_L';
  
  //for (i in 1:N_species) {temporal_trend[,i] = to_vector(year_range).*rep_vector(trend_slopes[i]/2,N_years) + GP_1D(year_range, temporal_trend_str_ls[i], temporal_trend_str_sd[i], temporal_trend_str_std[i]) + temporal_trend_unstr_sd[i]*temporal_trend_unstr_std[i];}
  //for (i in 1:N_species) {temporal_trend[,i] = to_vector(year_range).*rep_vector(trend_slopes[i]/2,N_years) + temporal_trend_unstr_sd[i]*temporal_trend_unstr_std[i];}
  for (i in 1:N_species) {temporal_trend[,i] = GP_1D(year_range, temporal_trend_str_ls[i], temporal_trend_str_sd[i], temporal_trend_str_std[i]) + temporal_trend_unstr_sd[i]*temporal_trend_unstr_std[i];}
  
  phenology_overall = phenology_bf*phenology_overall_std*phenology_overall_sd;
  for (i in 1:N_species) {phenology[,i] = phenology_overall + phenology_bf*phenology_std[i]*phenology_sd[i];}
  
  for (i in 1:N_spatial_dims) {spatial_Eta_str[,i] = GP_nD_isotropic(spatial_bf_range, spatial_Eta_str_ls[i], 1.0, spatial_Eta_str_std[,i]);}
  for (i in 1:N_spatiotemporal_dims) {spatiotemporal_Eta_str[,i] = GP_nD_anisotropic(spatiotemporal_bf_range, {spatiotemporal_Eta_str_ls_spatial[i], spatiotemporal_Eta_str_ls_spatial[i], spatiotemporal_Eta_str_ls_temporal[i]}, 1.0, spatiotemporal_Eta_str_std[,i]);}
  
  spatial_Lambda_sd = 1/sqrt(exp(cumulative_sum(log(spatial_Lambda_delta))));
  spatial_Lambda = spatial_Lambda_std .* rep_matrix(spatial_Lambda_sd, N_species);
  
  spatiotemporal_Lambda_sd = 1/sqrt(exp(cumulative_sum(log(spatiotemporal_Lambda_delta))));
  spatiotemporal_Lambda = spatiotemporal_Lambda_std .* rep_matrix(spatiotemporal_Lambda_sd, N_species);
  
  observer_Lambda_sd = 1/sqrt(exp(cumulative_sum(log(observer_Lambda_delta))));
  observer_Lambda = observer_Lambda_std .* rep_matrix(observer_Lambda_sd, N_species);
  
  // Key outcome computation
  
  for (i in 1:N_sites) {
    for (j in 1:N_years) {
      Z[i,j,] = bernoulli_logit_rng(temporal_trend[j,] + occupancy_X[i,]*occupancy_beta + (((1 - spatial_structured_fraction) .* spatial_Eta_unstr[i,]) + (spatial_structured_fraction .* ((spatial_bf[i,spatial_bf_indices[,1]] .* spatial_bf[i,spatial_bf_indices[,2]]) * spatial_Eta_str))) * spatial_Lambda + ((spatiotemporal_bf_spatial[i,spatiotemporal_bf_indices[,1]] .* spatiotemporal_bf_spatial[i,spatiotemporal_bf_indices[,2]] .* spatiotemporal_bf_temporal[j,spatiotemporal_bf_indices[,3]]) * spatiotemporal_Eta_str) * spatiotemporal_Lambda);
    }
  }
  
  for (i in 1:N_allvisits) {
    for (j in 1:N_species) {
      if (Z[site[i],year[i],j] == 0) {
        Y[i,j] = 0;
      } else {
        Y[i,j] = bernoulli_logit_rng(detection_X[i,]*detection_beta[,j] + phenology[week[i],j] + observer_Eta[observer[i],]*observer_Lambda[,j]);
      }
    }
  }
  
}
