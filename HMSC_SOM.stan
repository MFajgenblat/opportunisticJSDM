functions {
  
  vector GP_1D(real[] x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9))*eta);
  }
  
  vector GP_nD_isotropic(vector[] x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9))*eta);
  }
  
  real partial_sum(int[] visited_slice,
  int start, int end,
  int N_species,
  int[,] any_seen,
  int[,] Y,
  int[] N_visits,
  int[] start_visitid,
  int[] end_visitid,
  int[] year,
  int[] site,
  int[] day,
  int[] observer,
  matrix psi_X,
  matrix psi_beta,
  matrix psi_trend,
  matrix psi_Lambda_space,
  matrix psi_Eta_space,
  matrix psi_Lambda_spacetime,
  matrix psi_Eta_spacetime,
  real[] year_scaled,
  matrix p_X,
  matrix p_beta,
  matrix p_Eta_observer,
  matrix p_Lambda_observer,
  matrix p_phenology) {
    real ll = 0;
    for (i in visited_slice) {
      row_vector[N_species] psi_logit = psi_X[site[i],]*psi_beta + psi_trend[year[i],] + psi_Eta_space[site[i],]*psi_Lambda_space + psi_Eta_spacetime[site[i],]*psi_Lambda_spacetime.*rep_row_vector(year_scaled[year[i]], N_species);
      matrix[N_visits[i],N_species] p_logit = p_X[start_visitid[i]:end_visitid[i],]*p_beta + p_Eta_observer[observer[start_visitid[i]:end_visitid[i]],]*p_Lambda_observer + p_phenology[day[start_visitid[i]:end_visitid[i]],];
      for (j in 1:N_species) {
        if (any_seen[i,j] > 0) {
          ll += log_inv_logit(psi_logit[j]) + bernoulli_logit_lpmf(Y[start_visitid[i]:end_visitid[i],j] | p_logit[,j]);
        } else {
          ll += log_sum_exp(log_inv_logit(psi_logit[j]) + bernoulli_logit_lpmf(Y[start_visitid[i]:end_visitid[i],j] | p_logit[,j]), log1m_inv_logit(psi_logit[j]));
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
  int<lower=0> N_units;
  int<lower=0> N_visited;
  int visited[N_visited];
  int site[N_sites*N_years];
  int year[N_sites*N_years];
  int N_visits[N_sites*N_years];
  int any_seen[N_sites*N_years,N_species];
  int start_visitid[N_sites*N_years];
  int end_visitid[N_sites*N_years];
  int N_allvisits;
  int Y[N_allvisits,N_species];
  
  int N_days;
  int day[N_allvisits];
  int N_bf_phenology;
  matrix[N_days,N_bf_phenology] bf_phenology;
  
  int N_observers;
  int observer[N_allvisits];
  
  int N_p_dims;
  int N_psi_dims_space;
  int N_psi_dims_spacetime;
  
  int N_space_bf;
  matrix[N_sites,N_space_bf] space_bf;
  vector<lower=-1,upper=1>[2] space_bf_range[N_space_bf];
  
  int<lower=0> N_psi_X;
  matrix[N_sites,N_psi_X] psi_X;
  int<lower=0> N_p_X;
  matrix[N_allvisits,N_p_X] p_X;
  
  int<lower=0> N_Tr;
  matrix[N_Tr,N_species] Tr;
  
  corr_matrix[N_species] C;
  
  real year_range[N_years];
  real year_scaled[N_years];

}

transformed data {
  
  int grainsize = 1;
  
}

parameters {
  
  matrix[N_psi_X,N_Tr] psi_G;
  real<lower=0,upper=1> psi_a_phylo;
  cholesky_factor_corr[N_psi_X] psi_corr_covariates_L;
  vector<lower=0>[N_psi_X] psi_cov_covariates_scale;
  matrix[N_psi_X,N_species] psi_beta_std;
  
  matrix[N_p_X,N_Tr] p_G;
  real<lower=0,upper=1> p_a_phylo;
  cholesky_factor_corr[N_p_X] p_corr_covariates_L;
  vector<lower=0>[N_p_X] p_cov_covariates_scale;
  matrix[N_p_X,N_species] p_beta_std;
  
  real<lower=0> trend_species_ls[N_species];
  real<lower=0> trend_species_sd[N_species];
  vector[N_years] trend_species_std[N_species];
  
  matrix[N_psi_dims_space,N_species] psi_Lambda_space_std;
  vector<lower=0>[N_psi_dims_space] psi_Lambda_space_delta;
  real<lower=0,upper=1> psi_a_space[N_psi_dims_space];
  matrix[N_space_bf,N_psi_dims_space] psi_Eta_space_str_std;
  real<lower=0> psi_Eta_space_str_ls[N_psi_dims_space];
  matrix[N_sites,N_psi_dims_space] psi_Eta_space_unstr_std;
  
  matrix[N_psi_dims_spacetime,N_species] psi_Lambda_spacetime_std;
  vector<lower=0>[N_psi_dims_spacetime] psi_Lambda_spacetime_delta;
  real<lower=0,upper=1> psi_a_spacetime[N_psi_dims_spacetime];
  matrix[N_space_bf,N_psi_dims_spacetime] psi_Eta_spacetime_str_std;
  real<lower=0> psi_Eta_spacetime_str_ls[N_psi_dims_spacetime];
  matrix[N_sites,N_psi_dims_spacetime] psi_Eta_spacetime_unstr_std;
  
  matrix[N_p_dims,N_species] p_Lambda_observer_std;
  vector<lower=0>[N_p_dims] p_Lambda_observer_delta;
  matrix[N_observers,N_p_dims] p_Eta_observer_std;
  
  real<lower=0> p_phenology_overall_sd;
  vector[N_bf_phenology] p_phenology_overall_std;
  real<lower=0> p_phenology_sd[N_species];
  vector[N_bf_phenology] p_phenology_std[N_species];
  
}

transformed parameters {
  
  matrix[N_psi_dims_space,N_species] psi_Lambda_space;
  matrix[N_sites,N_psi_dims_space] psi_Eta_space;
  matrix[N_psi_dims_spacetime,N_species] psi_Lambda_spacetime;
  matrix[N_sites,N_psi_dims_spacetime] psi_Eta_spacetime;
  
  matrix[N_p_dims,N_species] p_Lambda_observer;
  matrix[N_observers,N_p_dims] p_Eta_observer;
  
  matrix[N_psi_X,N_species] psi_beta;
  cholesky_factor_cov[N_psi_X] psi_cov_covariates_L;
  cholesky_factor_corr[N_species] psi_cor_species_L;
  
  matrix[N_p_X,N_species] p_beta;
  cholesky_factor_cov[N_p_X] p_cov_covariates_L;
  cholesky_factor_corr[N_species] p_cor_species_L;
  
  matrix[N_years,N_species] psi_trend;
  
  vector[N_days] p_phenology_overall;
  matrix[N_days,N_species] p_phenology;
  
  psi_cov_covariates_L = quad_form_diag(psi_corr_covariates_L, psi_cov_covariates_scale);
  psi_cor_species_L = cholesky_decompose(psi_a_phylo*C + (1-psi_a_phylo)*identity_matrix(N_species));
  psi_beta = psi_G*Tr + psi_cov_covariates_L*psi_beta_std*psi_cor_species_L';
  
  p_cov_covariates_L = quad_form_diag(p_corr_covariates_L, p_cov_covariates_scale);
  p_cor_species_L = cholesky_decompose(p_a_phylo*C + (1-p_a_phylo)*identity_matrix(N_species));
  p_beta = p_G*Tr + p_cov_covariates_L*p_beta_std*p_cor_species_L';
  
  for (i in 1:N_species) {psi_trend[,i] = GP_1D(year_range, trend_species_ls[i], trend_species_sd[i], trend_species_std[i]);}
  
  psi_Lambda_space = psi_Lambda_space_std ./ rep_matrix(sqrt(exp(cumulative_sum(log(psi_Lambda_space_delta)))), N_species);
  for (i in 1:N_psi_dims_space) {psi_Eta_space[,i] = psi_a_space[i]*(space_bf*GP_nD_isotropic(space_bf_range, psi_Eta_space_str_ls[i], 1.0, psi_Eta_space_str_std[,i])) + (1 - psi_a_space[i])*psi_Eta_space_unstr_std[,i];}
  
  psi_Lambda_spacetime = psi_Lambda_spacetime_std ./ rep_matrix(sqrt(exp(cumulative_sum(log(psi_Lambda_spacetime_delta)))), N_species);
  for (i in 1:N_psi_dims_spacetime) {psi_Eta_spacetime[,i] = psi_a_spacetime[i]*space_bf*GP_nD_isotropic(space_bf_range, psi_Eta_spacetime_str_ls[i], 1.0, psi_Eta_spacetime_str_std[,i]) + (1 - psi_a_spacetime[i])*psi_Eta_spacetime_unstr_std[,i];}
  
  p_Lambda_observer = p_Lambda_observer_std ./ rep_matrix(sqrt(exp(cumulative_sum(log(p_Lambda_observer_delta)))), N_species);
  for (i in 1:N_p_dims) {p_Eta_observer[,i] = p_Eta_observer_std[,i];}
  
  p_phenology_overall = bf_phenology*p_phenology_overall_std*p_phenology_overall_sd;
  for (i in 1:N_species) {p_phenology[,i] = p_phenology_overall + bf_phenology*p_phenology_std[i]*p_phenology_sd[i];}
  
}

model {

  target += normal_lpdf(to_vector(psi_G) | 0, 3);
  target += uniform_lpdf(psi_a_phylo| 0, 1);
  target += lkj_corr_cholesky_lpdf(psi_corr_covariates_L | 2.0);
  target += normal_lpdf(psi_cov_covariates_scale | 0, 3);
  target += std_normal_lpdf(to_vector(psi_beta_std));

  target += normal_lpdf(to_vector(p_G) | 0, 3);
  target += uniform_lpdf(p_a_phylo | 0, 1);
  target += lkj_corr_cholesky_lpdf(p_corr_covariates_L | 2.0);
  target += normal_lpdf(p_cov_covariates_scale | 0, 3);
  target += std_normal_lpdf(to_vector(p_beta_std));
  
  for (i in 1:N_species) {
    target += std_normal_lpdf(trend_species_std[i]);
    target += normal_lpdf(trend_species_sd[i] | 0, 3);
    target += inv_gamma_lpdf(trend_species_ls[i] | 5, 5);
  }
  
  for (i in 1:N_psi_dims_space) {
    target += std_normal_lpdf(psi_Lambda_space_std[i,]);
    if (i == 1) {target += gamma_lpdf(psi_Lambda_space_delta[i] | 2, 1);}
    if (i > 1) {target += gamma_lpdf(psi_Lambda_space_delta[i] | 6, 1);}
    target += uniform_lpdf(psi_a_space[i] | 0, 1);
    target += std_normal_lpdf(psi_Eta_space_str_std[,i]);
    target += inv_gamma_lpdf(psi_Eta_space_str_ls[i] | 5, 5);
    target += std_normal_lpdf(psi_Eta_space_unstr_std[,i]);
  }
  
  for (i in 1:N_psi_dims_spacetime) {
    target += std_normal_lpdf(psi_Lambda_spacetime_std[i,]);
    if (i == 1) {target += gamma_lpdf(psi_Lambda_spacetime_delta[i] | 2, 1);}
    if (i > 1) {target += gamma_lpdf(psi_Lambda_spacetime_delta[i] | 6, 1);}
    target += uniform_lpdf(psi_a_spacetime[i] | 0, 1);
    target += std_normal_lpdf(psi_Eta_spacetime_str_std[,i]);
    target += inv_gamma_lpdf(psi_Eta_spacetime_str_ls[i] | 5, 5);
    target += std_normal_lpdf(psi_Eta_spacetime_unstr_std[,i]);
  }
  
  for (i in 1:N_p_dims) {
    target += std_normal_lpdf(p_Lambda_observer_std[i,]);
    if (i == 1) {target += gamma_lpdf(p_Lambda_observer_delta[i] | 2, 1);}
    if (i > 1) {target += gamma_lpdf(p_Lambda_observer_delta[i] | 6, 1);}
    target += std_normal_lpdf(p_Eta_observer_std[,i]);
  }
  
  target += normal_lpdf(p_phenology_overall_sd | 0, 3);
  target += std_normal_lpdf(p_phenology_overall_sd);
  
  for (i in 1:N_species) {
    target += normal_lpdf(p_phenology_sd[i] | 0, 3);
    target += std_normal_lpdf(p_phenology_std[i]);
  }
  
  target += reduce_sum(partial_sum, visited, grainsize, N_species, any_seen, Y, N_visits, start_visitid, end_visitid, year, site, day, observer, psi_X, psi_beta, psi_trend, psi_Lambda_space, psi_Eta_space, psi_Lambda_spacetime, psi_Eta_spacetime, year_scaled, p_X, p_beta, p_Eta_observer, p_Lambda_observer, p_phenology);
  
}

generated quantities {
  
//  real Z[N_units,N_species];

  matrix[N_psi_X,N_psi_X] psi_corr_covariates = multiply_lower_tri_self_transpose(psi_corr_covariates_L);
  matrix[N_p_X,N_p_X] p_corr_covariates = multiply_lower_tri_self_transpose(p_corr_covariates_L);
  matrix[N_species,N_species] psi_Omega_space = crossprod(psi_Lambda_space);
  matrix[N_species,N_species] psi_Omega_spacetime = crossprod(psi_Lambda_spacetime);
  matrix[N_species,N_species] p_Omega_observer = crossprod(p_Lambda_observer);
  
//  for (i in 1:N_units) {
//    row_vector[N_species] p_occupancy = inv_logit(psi_X[site[i],]*psi_beta + psi_trend[year[i],] + psi_Eta_space[site[i],]*psi_Lambda_space + psi_Eta_spacetime[site[i],]*psi_Lambda_spacetime.*rep_row_vector(year_scaled[year[i]], N_species));
//    if (N_visits[i] > 0) {
//      for (j in 1:N_species) {
//        if (any_seen[i,j] == 1) {
//          Z[i,j] = 1;
//        } else {
//          real p_nodetections_given_presence = exp(bernoulli_logit_lpmf(Y[start_visitid[i]:end_visitid[i],j] | p_X[start_visitid[i]:end_visitid[i],]*p_beta[,j] + p_Eta_observer[observer[start_visitid[i]:end_visitid[i]],]*p_Lambda_observer[,j] + p_phenology[day[start_visitid[i]:end_visitid[i]],j]));
//          Z[i,j] = (p_occupancy[j]*p_nodetections_given_presence)/((1 - p_occupancy[j]) + (p_occupancy[j]*p_nodetections_given_presence));
//        }
//      }
//    } else {
//      Z[i,] = to_array_1d(p_occupancy);
//    }
//  }

}
