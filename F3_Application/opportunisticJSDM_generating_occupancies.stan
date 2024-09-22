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
  
  int<lower=0> N_units_pred;
  int<lower=0> N_sites_pred;
  matrix[N_sites_pred,N_occupancy_covs] occupancy_X_pred;
  matrix[N_sites_pred,N_spatial_bf_uni] spatial_bf_pred;
  matrix[N_sites_pred,N_spatiotemporal_bf_uni_spatial] spatiotemporal_bf_spatial_pred;
  int site_pred[N_units_pred];
  int year_pred[N_units_pred];
  int<lower=0,upper=1> included[N_sites_pred];
  int<lower=0,upper=N_sites_pred> site_crossid[N_sites_pred];
  int<lower=0,upper=N_species> species_pred;

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
  
}

transformed parameters {
  
  matrix[N_occupancy_covs,N_species] occupancy_beta;
  cholesky_factor_cov[N_occupancy_covs] occupancy_cov_covariates_L;
  cholesky_factor_corr[N_species] occupancy_cor_species_L;
  
  matrix[N_detection_covs,N_species] detection_beta;
  cholesky_factor_cov[N_detection_covs] detection_cov_covariates_L;
  cholesky_factor_corr[N_species] detection_cor_species_L;
  
  matrix[N_years,N_species] temporal_trend;
  
  vector[N_weeks] phenology_overall;
  matrix[N_weeks,N_species] phenology;
  
  vector[N_spatial_dims] spatial_Lambda_sd = 1/sqrt(exp(cumulative_sum(log(spatial_Lambda_delta))));
  matrix[N_spatial_dims,N_species] spatial_Lambda = spatial_Lambda_std .* rep_matrix(spatial_Lambda_sd, N_species);
  
  vector[N_spatiotemporal_dims] spatiotemporal_Lambda_sd = 1/sqrt(exp(cumulative_sum(log(spatiotemporal_Lambda_delta))));
  matrix[N_spatiotemporal_dims,N_species] spatiotemporal_Lambda = spatiotemporal_Lambda_std .* rep_matrix(spatiotemporal_Lambda_sd, N_species);
  
  vector[N_observer_dims] observer_Lambda_sd = 1/sqrt(exp(cumulative_sum(log(observer_Lambda_delta))));
  matrix[N_observer_dims,N_species] observer_Lambda = observer_Lambda_std .* rep_matrix(observer_Lambda_sd, N_species);
  
  matrix[N_spatial_bf,N_spatial_dims] spatial_Eta_str;
  matrix[N_spatiotemporal_bf,N_spatiotemporal_dims] spatiotemporal_Eta_str;
  
  occupancy_cov_covariates_L = quad_form_diag(occupancy_cor_covariates_L, occupancy_cov_covariates_scale);
  occupancy_cor_species_L = cholesky_decompose(occupancy_phylo_fraction*C + (1-occupancy_phylo_fraction)*identity_matrix(N_species));
  occupancy_beta = occupancy_G*Tr + occupancy_cov_covariates_L*occupancy_beta_std*occupancy_cor_species_L';
  
  detection_cov_covariates_L = quad_form_diag(detection_cor_covariates_L, detection_cov_covariates_scale);
  detection_cor_species_L = cholesky_decompose(detection_phylo_fraction*C + (1-detection_phylo_fraction)*identity_matrix(N_species));
  detection_beta = detection_G*Tr + detection_cov_covariates_L*detection_beta_std*detection_cor_species_L';
  
  for (i in 1:N_species) {temporal_trend[,i] = GP_1D(year_range, temporal_trend_str_ls[i], temporal_trend_str_sd[i], temporal_trend_str_std[i]) + temporal_trend_unstr_sd[i]*temporal_trend_unstr_std[i];}
  
  phenology_overall = phenology_bf*phenology_overall_std*phenology_overall_sd;
  for (i in 1:N_species) {phenology[,i] = phenology_overall + phenology_bf*phenology_std[i]*phenology_sd[i];}
  
  for (i in 1:N_spatial_dims) {spatial_Eta_str[,i] = GP_nD_isotropic(spatial_bf_range, spatial_Eta_str_ls[i], 1.0, spatial_Eta_str_std[,i]);}
  for (i in 1:N_spatiotemporal_dims) {spatiotemporal_Eta_str[,i] = GP_nD_anisotropic(spatiotemporal_bf_range, {spatiotemporal_Eta_str_ls_spatial[i], spatiotemporal_Eta_str_ls_spatial[i], spatiotemporal_Eta_str_ls_temporal[i]}, 1.0, spatiotemporal_Eta_str_std[,i]);}
  
}


generated quantities {
  
  real Z[N_units_pred];
  
  for (i in 1:N_units_pred) {
    if (species_pred > 0) {
      if (included[site_pred[i]] == 1) {
        Z[i] = inv_logit(occupancy_X_pred[site_pred[i],]*occupancy_beta + temporal_trend[year_pred[i],] + (((1 - spatial_structured_fraction) .* spatial_Eta_unstr[site_crossid[site_pred[i]],]) + (spatial_structured_fraction .* ((spatial_bf_pred[site_pred[i],spatial_bf_indices[,1]] .* spatial_bf_pred[site_pred[i],spatial_bf_indices[,2]]) * spatial_Eta_str))) * spatial_Lambda + ((spatiotemporal_bf_spatial_pred[site_pred[i],spatiotemporal_bf_indices[,1]] .* spatiotemporal_bf_spatial_pred[site_pred[i],spatiotemporal_bf_indices[,2]] .* spatiotemporal_bf_temporal[year_pred[i],spatiotemporal_bf_indices[,3]]) * spatiotemporal_Eta_str) * spatiotemporal_Lambda)[species_pred];
      } else {
        Z[i] = inv_logit(occupancy_X_pred[site_pred[i],]*occupancy_beta + temporal_trend[year_pred[i],] + ((spatial_structured_fraction .* ((spatial_bf_pred[site_pred[i],spatial_bf_indices[,1]] .* spatial_bf_pred[site_pred[i],spatial_bf_indices[,2]]) * spatial_Eta_str))) * spatial_Lambda + ((spatiotemporal_bf_spatial_pred[site_pred[i],spatiotemporal_bf_indices[,1]] .* spatiotemporal_bf_spatial_pred[site_pred[i],spatiotemporal_bf_indices[,2]] .* spatiotemporal_bf_temporal[year_pred[i],spatiotemporal_bf_indices[,3]]) * spatiotemporal_Eta_str) * spatiotemporal_Lambda)[species_pred];
      }
    } else if (species_pred == 0) {
      if (included[site_pred[i]] == 1) {
        Z[i] = sum(inv_logit(occupancy_X_pred[site_pred[i],]*occupancy_beta + temporal_trend[year_pred[i],] + (((1 - spatial_structured_fraction) .* spatial_Eta_unstr[site_crossid[site_pred[i]],]) + (spatial_structured_fraction .* ((spatial_bf_pred[site_pred[i],spatial_bf_indices[,1]] .* spatial_bf_pred[site_pred[i],spatial_bf_indices[,2]]) * spatial_Eta_str))) * spatial_Lambda + ((spatiotemporal_bf_spatial_pred[site_pred[i],spatiotemporal_bf_indices[,1]] .* spatiotemporal_bf_spatial_pred[site_pred[i],spatiotemporal_bf_indices[,2]] .* spatiotemporal_bf_temporal[year_pred[i],spatiotemporal_bf_indices[,3]]) * spatiotemporal_Eta_str) * spatiotemporal_Lambda));
      } else {
        Z[i] = sum(inv_logit(occupancy_X_pred[site_pred[i],]*occupancy_beta + temporal_trend[year_pred[i],] + ((spatial_structured_fraction .* ((spatial_bf_pred[site_pred[i],spatial_bf_indices[,1]] .* spatial_bf_pred[site_pred[i],spatial_bf_indices[,2]]) * spatial_Eta_str))) * spatial_Lambda + ((spatiotemporal_bf_spatial_pred[site_pred[i],spatiotemporal_bf_indices[,1]] .* spatiotemporal_bf_spatial_pred[site_pred[i],spatiotemporal_bf_indices[,2]] .* spatiotemporal_bf_temporal[year_pred[i],spatiotemporal_bf_indices[,3]]) * spatiotemporal_Eta_str) * spatiotemporal_Lambda));
      }
    }
  }

}
