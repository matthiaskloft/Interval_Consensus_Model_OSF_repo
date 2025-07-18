//////////////////////////////////////////////////////////////////////////////// 
// Interval Truth Model
//////////////////////////////////////////////////////////////////////////////// 
functions {
  // simplex to bivariate normal, array to matrix
  matrix simplex_to_bvn(array[] vector Y_splx) {
    
    int N = size(Y_splx);
    matrix[N, 2] Y;
    
    vector[N] log_ratio_1 = log(to_vector(Y_splx[, 1]) ./ to_vector(Y_splx[, 3]));
    vector[N] log_ratio_2 = log(to_vector(Y_splx[, 2]) ./ sqrt(to_vector(Y_splx[, 1]) .* to_vector(Y_splx[, 3])));

    Y[, 1] = sqrt(1.0 / 2) * log_ratio_1;
    Y[, 2] = sqrt(2.0 / 3) * log_ratio_2;
    
    return Y;
  }
  // simplex to bivariate normal, matrix to matrix
  matrix simplex_to_bvn(matrix Y_splx) {
    
    int N = rows(Y_splx);
    matrix[N, 2] Y;
    
    vector[N] log_ratio_1 = log(Y_splx[, 1] ./ Y_splx[, 3]);
    vector[N] log_ratio_2 = log(Y_splx[, 2] ./ sqrt(Y_splx[, 1] .* Y_splx[, 3]));

    Y[, 1] = sqrt(1.0 / 2) * log_ratio_1;
    Y[, 2] = sqrt(2.0 / 3) * log_ratio_2;
    
    return Y;
  }
  
    // bivariate normal to simplex, matrix to matrix
    matrix bvn_to_simplex(matrix Y) {
    
    int N = rows(Y);
    matrix[N, 3] Y_splx;
    
    vector[N] Sum = exp(sqrt(2) .* Y[,1]) + 
                    exp(sqrt(3.0/2) .* Y[,2] + Y[,1] ./ sqrt(2)) + 
                    1;
    Y_splx[ ,1] = exp(sqrt(2) .* Y[,1])                        ./ Sum .* 1.03 - 0.01;
    Y_splx[ ,2] = exp(sqrt(3.0/2) .* Y[,2] + Y[,1] ./ sqrt(2)) ./ Sum .* 1.03 - 0.01;
    Y_splx[ ,3] = 1                                            ./ Sum .* 1.03 - 0.01;
    
    return Y_splx;
  }
  
}
//////////////////////////////////////////////////////////////////////////////// 
data{
  int<lower=1> I; // number of persons
  int<lower=1> J; // number of items
  int<lower=1> N; // number of observed responses
  array[N] int<lower=1> ii; // person indices
  array[N] int<lower=1> jj; // item indices
  array[N] int<lower=1> nn; // response indices
  array[N] simplex[3] Y_splx ; // DRS responses as simplex
  real<lower=0,upper=1> padding; // padding constant used in zero handling; set to 0 if not applicable
}
//////////////////////////////////////////////////////////////////////////////// 
transformed data {
  // scale maximum for back transformation to simplex
  real scale_max = 1 + padding * 3;
  // arrray to array
  array[N] vector[2] Y;
  for (n in 1:N){
    Y[n,1] = sqrt(1.0/2) * log(Y_splx[n,1] / Y_splx[n,3]);
    Y[n,2] = sqrt(2.0/3) * log(Y_splx[n,2] / sqrt(Y_splx[n,1] * Y_splx[n,3]));
  }
}
//////////////////////////////////////////////////////////////////////////////// 
parameters{
// person patameters
  array[J] simplex[3] Tr_splx_model; // true intervals in simplex format
  matrix[I,5] I_raw;  
  
  // hyperpriors person
  vector[2] mu_E;
  vector[5] sigma_I;
  cholesky_factor_corr[2] L_corr_E; // correlation between proficiencies for location and width
  
// item parameters
  // raws for non-centered parameterization
  matrix[J,3] J_raw;
  cholesky_factor_corr[2] L_corr_lambda; // correlation between item discernibilities for location and width

  // hyperpriors item
  vector[2] sigma_lambda;
  vector<lower=0,upper=1> [J] omega_beta; // residual correlations
}
//////////////////////////////////////////////////////////////////////////////// 
transformed parameters{
  // person parameters
  vector[I] E_loc; // competence location
  vector[I] E_wid; // competence width
  vector[I] a_loc; // scaling bias location
  vector[I] b_loc; // shifting bias location
  vector[I] b_wid; // shifting bias width
  // item parameters
  matrix[J,3] Tr_splx_model;
  vector[J] Tr_loc; // latent truth location
  vector[J] Tr_wid; // latent truth width
  vector[J] lambda_loc; // item difficulty / discernibility location
  vector[J] lambda_wid; // item difficulty / discernibility width
  // correlation residual
  vector[J] omega;

  // person parameters
  {
  matrix[I,2] E_temp = (diag_pre_multiply(exp(sigma_I[1:2] .* .5 + log(.5)), L_corr_E) * I_raw[ ,1:2]')';
  E_loc = exp(E_temp[,1] + mu_E[1]);
  E_wid = exp(E_temp[,2] + mu_E[2]);
  }
  a_loc = exp(I_raw[,3] * exp(sigma_I[3] * .5 + log(.5))); // logN(log(0.5),0.5)
  b_loc =     I_raw[,4] * exp(sigma_I[4]      + log(.5));  // logN(log(0.5),1)
  b_wid =     I_raw[,5] * exp(sigma_I[5]      + log(.5));  // logN(log(0.5),1)

  //item parameters
  {  
  // transform marginal locations and  widths into simplex  
  Tr_splx_model[,1] = (1 - Tr_wid_beta) .* Tr_loc_beta;
  Tr_splx_model[,2] = Tr_wid_beta;
  Tr_splx_model[,3] = 1 - Tr_splx_model[,1] - Tr_splx_model[,2];
  // transform simplex to bivariate normal
  matrix[J,2] Tr_bvn = simplex_to_bvn(Tr_splx_model);
  // final bivariate locations and widths
  Tr_loc             = Tr_bvn[,1];  // latent truth
  Tr_wid             = Tr_bvn[,2];  // latent truth
  
  matrix[J,2] lambda_temp = (diag_pre_multiply(exp(sigma_lambda[1:2] * .5 + log(.5)), L_corr_lambda) * J_raw[ ,1:2]')';
  lambda_loc = exp(lambda_temp[,1]); // item difficulty / discernibility
  lambda_wid = exp(lambda_temp[,2]); // item difficulty / discernibility
  }
  // correlation residual
  omega = omega_beta * 2 - 1;

}
//////////////////////////////////////////////////////////////////////////////// 
model{
  // raws
  Tr_splx_model ~ dirichlet([1,1,1]);
  
  to_vector(I_raw) ~ std_normal();
  to_vector(J_raw) ~ std_normal();

  L_corr_E         ~ lkj_corr_cholesky(2);
  L_corr_lambda    ~ lkj_corr_cholesky(2);
  
  // hyper priors means
  mu_E             ~ std_normal();
  // hyper priors scales
  sigma_I          ~ std_normal();
  sigma_lambda     ~ std_normal();
  // prior correlation
  omega_beta       ~ beta(2,2);

// Model //
  {
    // parameters for MVN
    matrix[N,2] mu; 
    matrix[N,2] sigma;
    // Mean Vector  
    mu[ ,1]    = Tr_loc[jj] .* a_loc[ii];
    mu[ ,2]    = Tr_wid[jj] + b_wid[ii]; 
    sigma[ ,1] = exp(log(a_loc[ii]) - log(E_loc[ii]) - log(lambda_loc[jj]));
    sigma[ ,2] = exp(-log(E_wid[ii]) - log(lambda_wid[jj]));
   
  for(n in 1:N){
      // correlation between location and width
      matrix[2,2] Omega;
      // Correlation Matrix
      Omega[1,1] = 1;
      Omega[2,2] = 1;
      Omega[1,2] = omega[jj[n]];
      Omega[2,1] = Omega[1,2];
      matrix[2,2] Sigma = quad_form_diag(Omega, sigma[n,]);
    // Likelihood
    Y[n] ~ multi_normal(mu[n,], Sigma);
    } //end n
  } // end block
}

//////////////////////////////////////////////////////////////////////////////// 
 generated quantities{
   // person parameter correlations
  real rho_E = multiply_lower_tri_self_transpose(L_corr_E)[1,2];
  // item parameter correlations
  real rho_lambda = multiply_lower_tri_self_transpose(L_corr_lambda)[1,2];

  // Posterior Predicted Values
  vector[N] Y_ppc_loc;
  vector[N] Y_ppc_wid;
  matrix[N,3] Y_ppc_splx;
  vector[N] Y_ppc_loc_splx;
  vector[N] Y_ppc_wid_splx;
  
  {
    // Posterior Predicted Values
    matrix[N,2] Y_ppc;
    // parameters for MVN
    matrix[N,2] mu; 
    matrix[N,2] sigma;
    // Mean Vector  
    mu[ ,1]    = Tr_loc[jj] .* a_loc[ii];
    mu[ ,2]    = Tr_wid[jj] + b_wid[ii]; 
    sigma[ ,1] = exp(log(a_loc[ii]) - log(E_loc[ii]) - log(lambda_loc[jj]));
    sigma[ ,2] = exp(-log(E_wid[ii]) - log(lambda_wid[jj]));
   
    for(n in 1:N){
      // correlation between location and width
      matrix[2,2] Omega;
      // Correlation Matrix
      Omega[1,1] = 1;
      Omega[2,2] = 1;
      Omega[1,2] = omega[jj[n]];
      Omega[2,1] = Omega[1,2];
      matrix[2,2] Sigma = quad_form_diag(Omega, sigma[n,]);
      // Predicted Responses
      Y_ppc[n] = multi_normal_rng(mu[n,], Sigma)';
    } //end n
  
    // Posterior Predicted Values for unbounded location and width  
    Y_ppc_loc = Y_ppc[,1];
    Y_ppc_wid = Y_ppc[,2];
    Y_ppc_splx = bvn_to_simplex(Y_ppc);
    Y_ppc_loc_splx = Y_ppc_splx[,1] + 0.5 .* Y_ppc_splx[,2];
    Y_ppc_wid_splx = Y_ppc_splx[,2];
  } // end block
  
  // Latent consensus simplex with reversed padding
  matrix[J,3] Tr_splx;
  Tr_splx[,1] = Tr_splx_model[,1] .* scale_max - padding;
  Tr_splx[,2] = Tr_splx_model[,2] .* scale_max - padding;
  Tr_splx[,3] = Tr_splx_model[,3] .* scale_max - padding;
  vector[J] Tr_loc_splx = Tr_splx[,1] + 0.5 .* Tr_splx[,2];
  vector[J] Tr_wid_splx = Tr_splx[,2];
}

