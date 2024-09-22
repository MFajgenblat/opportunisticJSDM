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
    int[] week,
    int[] observer,
    matrix occupancy_X,
    matrix occupancy_beta,
    matrix detection_X,
    matrix detection_beta,
    matrix temporal_trend,
    matrix phenology,
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
    int[,] spatiotemporal_bf_indices,
    matrix observer_Lambda,
    matrix observer_Eta) {
      real ll = 0;
      for (i in visited_units_slice) {
        row_vector[N_species] occupancy_logit = occupancy_X[site[i],]*occupancy_beta + temporal_trend[year[i],] + (((1 - spatial_structured_fraction) .* spatial_Eta_unstr[site[i],]) + (spatial_structured_fraction .* ((spatial_bf[site[i],spatial_bf_indices[,1]] .* spatial_bf[site[i],spatial_bf_indices[,2]]) * spatial_Eta_str))) * spatial_Lambda + ((spatiotemporal_bf_spatial[site[i],spatiotemporal_bf_indices[,1]] .* spatiotemporal_bf_spatial[site[i],spatiotemporal_bf_indices[,2]] .* spatiotemporal_bf_temporal[year[i],spatiotemporal_bf_indices[,3]]) * spatiotemporal_Eta_str) * spatiotemporal_Lambda;
        matrix[N_visits[i],N_species] detection_logit = detection_X[start_visitid[i]:end_visitid[i],]*detection_beta + phenology[week[start_visitid[i]:end_visitid[i]],] + observer_Eta[observer[start_visitid[i]:end_visitid[i]],]*observer_Lambda;
        for (j in 1:N_species) {
          if (any_seen[i,j] > 0) {
            ll += log_inv_logit(occupancy_logit[j]) + bernoulli_logit_lpmf(Y[start_visitid[i]:end_visitid[i],j] | detection_logit[,j]);
          } else {
            ll += log_sum_exp(log_inv_logit(occupancy_logit[j]) + bernoulli_logit_lpmf(Y[start_visitid[i]:end_visitid[i],j] | detection_logit[,j]), log1m_inv_logit(occupancy_logit[j]));
          }
        }
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
  int site[N_visited_units];
  int year[N_visited_units];
  int N_visits[N_visited_units];
  int any_seen[N_visited_units,N_species];
  int start_visitid[N_visited_units];
  int end_visitid[N_visited_units];
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
  
  int<lower=0> N_gradientpred;
  matrix[N_gradientpred,N_occupancy_covs] occupancy_X_gradient;

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
  
}

transformed parameters {
  
  matrix[N_occupancy_covs,N_species] occupancy_beta;
  cholesky_factor_cov[N_occupancy_covs] occupancy_cov_covariates_L;
  cholesky_factor_corr[N_species] occupancy_cor_species_L;
  
  occupancy_cov_covariates_L = quad_form_diag(occupancy_cor_covariates_L, occupancy_cov_covariates_scale);
  occupancy_cor_species_L = cholesky_decompose(occupancy_phylo_fraction*C + (1-occupancy_phylo_fraction)*identity_matrix(N_species));
  occupancy_beta = occupancy_G*Tr + occupancy_cov_covariates_L*occupancy_beta_std*occupancy_cor_species_L';
  
}


generated quantities {
  
  real Z[N_gradientpred];
  
  for (i in 1:N_gradientpred) {
    Z[i] = sum(inv_logit(occupancy_X_gradient[i,]*occupancy_beta));
  }

}
