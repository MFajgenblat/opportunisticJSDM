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
  
  real partial_sum(int[] visited_units_slice,
    int start, int end,
    int N_species,
    int[] N_visits,
    int[,] any_seen,
    int[,] Y,
    int[] start_visitid,
    int[] end_visitid,
    int[] year,
    int[] site,
    matrix occupancy_X,
    matrix occupancy_beta,
    matrix temporal_trend,
    matrix spatial_Lambda,
    matrix spatial_Eta_str,
    matrix spatial_bf,
    int[,] spatial_bf_indices,
    matrix spatial_Eta_unstr,
    row_vector spatial_structured_fraction,
    matrix spatiotemporal_Lambda,
    matrix spatiotemporal_Eta_str,
    matrix spatiotemporal_bf_spatial,
    matrix spatiotemporal_bf_temporal,
    int[,] spatiotemporal_bf_indices) {
      real ll = 0;
      for (i in visited_units_slice) {
        ll += bernoulli_logit_lpmf(any_seen[i,] | occupancy_X[site[i],]*occupancy_beta + temporal_trend[year[i],] + (((1 - spatial_structured_fraction) .* spatial_Eta_unstr[site[i],]) + (spatial_structured_fraction .* ((spatial_bf[site[i],spatial_bf_indices[,1]] .* spatial_bf[site[i],spatial_bf_indices[,2]]) * spatial_Eta_str))) * spatial_Lambda + ((spatiotemporal_bf_spatial[site[i],spatiotemporal_bf_indices[,1]] .* spatiotemporal_bf_spatial[site[i],spatiotemporal_bf_indices[,2]] .* spatiotemporal_bf_temporal[year[i],spatiotemporal_bf_indices[,3]]) * spatiotemporal_Eta_str) * spatiotemporal_Lambda);
      }
    return ll;
  } 
  
}

data {
  
  int<lower=0> N_species;
  int<lower=0> N_sites;
  int<lower=0> N_years;
  int<lower=0> N_visited_units;
  int visited_units[N_visited_units];
  int site[N_sites*N_years];
  int year[N_sites*N_years];
  int N_visits[N_sites*N_years];
  int any_seen[N_sites*N_years,N_species];
  int start_visitid[N_sites*N_years];
  int end_visitid[N_sites*N_years];
  int N_allvisits;
  int Y[N_allvisits,N_species];
  
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

transformed data {
  
  int grainsize = 1;
  
}

parameters {
  
  matrix[N_occupancy_covs,N_traits] occupancy_G;
  real<lower=0,upper=1> occupancy_phylo_fraction;
  cholesky_factor_corr[N_occupancy_covs] occupancy_cor_covariates_L;
  vector<lower=0>[N_occupancy_covs] occupancy_cov_covariates_scale;
  matrix[N_occupancy_covs,N_species] occupancy_beta_std;
  
  real<lower=0> temporal_trend_unstr_sd[N_species];
  vector[N_years] temporal_trend_unstr_std[N_species];
  real<lower=0> temporal_trend_str_ls[N_species];
  real<lower=0> temporal_trend_str_sd[N_species];
  vector[N_years] temporal_trend_str_std[N_species];
  
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
  
}

transformed parameters {
  
  matrix[N_occupancy_covs,N_species] occupancy_beta;
  cholesky_factor_cov[N_occupancy_covs] occupancy_cov_covariates_L;
  cholesky_factor_corr[N_species] occupancy_cor_species_L;
  
  matrix[N_years,N_species] temporal_trend;
  
  vector[N_spatial_dims] spatial_Lambda_sd = 1/sqrt(exp(cumulative_sum(log(spatial_Lambda_delta))));
  matrix[N_spatial_dims,N_species] spatial_Lambda = spatial_Lambda_std .* rep_matrix(spatial_Lambda_sd, N_species);
  
  vector[N_spatiotemporal_dims] spatiotemporal_Lambda_sd = 1/sqrt(exp(cumulative_sum(log(spatiotemporal_Lambda_delta))));
  matrix[N_spatiotemporal_dims,N_species] spatiotemporal_Lambda = spatiotemporal_Lambda_std .* rep_matrix(spatiotemporal_Lambda_sd, N_species);
  
  matrix[N_spatial_bf,N_spatial_dims] spatial_Eta_str;
  matrix[N_spatiotemporal_bf,N_spatiotemporal_dims] spatiotemporal_Eta_str;
  
  occupancy_cov_covariates_L = quad_form_diag(occupancy_cor_covariates_L, occupancy_cov_covariates_scale);
  occupancy_cor_species_L = cholesky_decompose(occupancy_phylo_fraction*C + (1-occupancy_phylo_fraction)*identity_matrix(N_species));
  occupancy_beta = occupancy_G*Tr + occupancy_cov_covariates_L*occupancy_beta_std*occupancy_cor_species_L';
  
  for (i in 1:N_species) {temporal_trend[,i] = GP_1D(year_range, temporal_trend_str_ls[i], temporal_trend_str_sd[i], temporal_trend_str_std[i]) + temporal_trend_unstr_sd[i]*temporal_trend_unstr_std[i];}
  
  for (i in 1:N_spatial_dims) {spatial_Eta_str[,i] = GP_nD_isotropic(spatial_bf_range, spatial_Eta_str_ls[i], 1.0, spatial_Eta_str_std[,i]);}
  for (i in 1:N_spatiotemporal_dims) {spatiotemporal_Eta_str[,i] = GP_nD_anisotropic(spatiotemporal_bf_range, {spatiotemporal_Eta_str_ls_spatial[i], spatiotemporal_Eta_str_ls_spatial[i], spatiotemporal_Eta_str_ls_temporal[i]}, 1.0, spatiotemporal_Eta_str_std[,i]);}
  
  
}

model {

  target += normal_lpdf(to_vector(occupancy_G) | 0, 3);
  target += uniform_lpdf(occupancy_phylo_fraction| 0, 1);
  target += lkj_corr_cholesky_lpdf(occupancy_cor_covariates_L | 2.0);
  target += normal_lpdf(occupancy_cov_covariates_scale | 0, 3);
  target += std_normal_lpdf(to_vector(occupancy_beta_std));

  for (i in 1:N_species) {
    target += normal_lpdf(temporal_trend_str_sd[i] | 0, 3);
    target += std_normal_lpdf(temporal_trend_str_std[i]);
    target += inv_gamma_lpdf(temporal_trend_str_ls[i] | 5, 5);
    target += normal_lpdf(temporal_trend_unstr_sd[i] | 0, 3);
    target += std_normal_lpdf(temporal_trend_unstr_std[i]);
  }
  
  for (i in 1:N_spatial_dims) {
    target += std_normal_lpdf(spatial_Lambda_std[i,]);
    if (i == 1) {target += gamma_lpdf(spatial_Lambda_delta[i] | 2, 1);}
    if (i > 1) {target += gamma_lpdf(spatial_Lambda_delta[i] | 6, 1);}
    target += uniform_lpdf(spatial_structured_fraction[i] | 0, 1);
    target += std_normal_lpdf(spatial_Eta_str_std[,i]);
    target += inv_gamma_lpdf(spatial_Eta_str_ls[i] | 5, 5);
    target += std_normal_lpdf(spatial_Eta_unstr[,i]);
  }
  
  for (i in 1:N_spatiotemporal_dims) {
    target += std_normal_lpdf(spatiotemporal_Lambda_std[i,]);
    if (i == 1) {target += gamma_lpdf(spatiotemporal_Lambda_delta[i] | 2, 1);}
    if (i > 1) {target += gamma_lpdf(spatiotemporal_Lambda_delta[i] | 6, 1);}
    target += std_normal_lpdf(spatiotemporal_Eta_str_std[,i]);
    target += inv_gamma_lpdf(spatiotemporal_Eta_str_ls_spatial[i] | 5, 5);
    target += inv_gamma_lpdf(spatiotemporal_Eta_str_ls_temporal[i] | 5, 5);
  }
  
  target += reduce_sum(partial_sum, visited_units, grainsize,
    N_species, N_visits,
    any_seen, Y,
    start_visitid, end_visitid,
    year, site,
    occupancy_X, occupancy_beta,
    temporal_trend,
    spatial_Lambda, spatial_Eta_str, spatial_bf, spatial_bf_indices, spatial_Eta_unstr, spatial_structured_fraction,
    spatiotemporal_Lambda, spatiotemporal_Eta_str, spatiotemporal_bf_spatial, spatiotemporal_bf_temporal, spatiotemporal_bf_indices);
  
}

generated quantities {
  
  row_vector<lower=0,upper=1>[N_species] Z[N_sites,N_years];
  
  for (i in 1:N_sites) {
    for (j in 1:N_years) {
      Z[i,j] = inv_logit(occupancy_X[i,]*occupancy_beta + temporal_trend[j,] + (((1 - spatial_structured_fraction) .* spatial_Eta_unstr[i,]) + (spatial_structured_fraction .* ((spatial_bf[i,spatial_bf_indices[,1]] .* spatial_bf[i,spatial_bf_indices[,2]]) * spatial_Eta_str))) * spatial_Lambda + ((spatiotemporal_bf_spatial[i,spatiotemporal_bf_indices[,1]] .* spatiotemporal_bf_spatial[i,spatiotemporal_bf_indices[,2]] .* spatiotemporal_bf_temporal[j,spatiotemporal_bf_indices[,3]]) * spatiotemporal_Eta_str) * spatiotemporal_Lambda);
    }
  }
  
}

