//////////////////////////////////////////////////////////////////////////////// 
// Interval Truth Model
//////////////////////////////////////////////////////////////////////////////// 
functions {
  // simplex to bivariate normal, array to matrix
  matrix simplex_to_bvn(array[] vector Y_splx) {
    
    int N = size(Y_splx);
    matrix[N, 2] Y;
    
    vector[N] log_ratio_1 = 
      log(to_vector(Y_splx[, 1]) ./ to_vector(Y_splx[, 3]));
    vector[N] log_ratio_2 = 
      log(to_vector(Y_splx[, 2]) ./ 
      sqrt(to_vector(Y_splx[, 1]) .* to_vector(Y_splx[, 3])));

    Y[, 1] = sqrt(1.0 / 2) .* log_ratio_1;
    Y[, 2] = sqrt(2.0 / 3) .* log_ratio_2;
    
    return Y;
  }
  // simplex to bivariate normal, matrix to matrix
  matrix simplex_to_bvn(matrix Y_splx) {
    
    int N = rows(Y_splx);
    matrix[N, 2] Y;
    
    vector[N] log_ratio_1 = log(Y_splx[, 1] ./ Y_splx[, 3]);
    vector[N] log_ratio_2 = log(Y_splx[, 2] ./ sqrt(Y_splx[, 1] .* Y_splx[, 3]));

    Y[, 1] = sqrt(1.0 / 2) .* log_ratio_1;
    Y[, 2] = sqrt(2.0 / 3) .* log_ratio_2;
    
    return Y;
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
  int<lower=1,upper=2> link; // link function: 1 = ILR, 2 = SB
}
//////////////////////////////////////////////////////////////////////////////// 
transformed data {
  // simplex to bvn, array to array
  array[N] vector[2] Y;
  
  for (n in 1:N){
    if(link == 1){
      // ILR
      Y[n,1] = sqrt(1.0/2) * log(Y_splx[n,1] / Y_splx[n,3]);
      Y[n,2] = sqrt(2.0/3) * log(Y_splx[n,2] / sqrt(Y_splx[n,1] * Y_splx[n,3]));
    }else{
      // SB
      Y[n,1] = log(Y_splx[n,1] / Y_splx[n,3]);
      Y[n,2] = log(Y_splx[n,2] / (Y_splx[n,1] + Y_splx[n,3]));
    } // end if
  } // end n
}
//////////////////////////////////////////////////////////////////////////////// 
parameters{
// person patameters
  vector<lower=0,upper=1> [J] Tr_loc_beta; // marginal locations on bounded scale
  vector<lower=0,upper=1> [J] Tr_wid_beta; // marginal widths on bounded scale
  matrix[I,5] I_raw;
  // hyperpriors person
  vector[2] mu_I;
  vector[5] sigma_I;
// item parameters
  // raws for non-centered parameterization
  matrix[J,2] J_raw;
  // hyperpriors item
  vector[2] sigma_J;
  vector <lower=0,upper=1>[J] omega_beta;
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

  //// person parameters ////
  E_loc = exp(I_raw[,1] * exp(sigma_I[1] * .5 + log(.5)) + mu_I[1]); // N(mu_I,logN(log(0.5),0.5))
  E_wid = exp(I_raw[,2] * exp(sigma_I[2] * .5 + log(.5)) + mu_I[2]); // N(mu_I,logN(log(0.5),0.5))
  a_loc = exp(I_raw[,3] * exp(sigma_I[3] * .5 + log(.5))); // N(0,logN(log(0.5),0.5))
  b_loc =     I_raw[,4] * exp(sigma_I[4]      + log(.5)); // N(0,logN(log(0.5),1))
  b_wid =     I_raw[,5] * exp(sigma_I[5]      + log(.5)); // N(0,logN(log(0.5),1))
  
  //// item parameters ////
  // transform marginal locations and  widths into simplex  
  Tr_splx_model[,1] = (1 - Tr_wid_beta) .* Tr_loc_beta;
  Tr_splx_model[,2] = Tr_wid_beta;
  Tr_splx_model[,3] = 1 - Tr_splx_model[,1] - Tr_splx_model[,2];
  {
  // transform simplex to bivariate normal
  matrix[J,2] Tr_bvn = simplex_to_bvn(Tr_splx_model);
  Tr_loc             = Tr_bvn[,1];  // latent true location
  Tr_wid             = Tr_bvn[,2];  // latent true width
  }
  lambda_loc =  exp(J_raw[,1] * exp(sigma_J[1] * .5 + log(.5)));
  lambda_wid =  exp(J_raw[,2] * exp(sigma_J[2] * .5 + log(.5)));
  // residual correlation
  omega      =  omega_beta .* 2 - 1;
}
//////////////////////////////////////////////////////////////////////////////// 
model{
  // Priors //
  // marginal true locations and widths
  Tr_loc_beta    ~ beta(1,1); // marginal locations on bounded scale
  Tr_wid_beta    ~ beta(1.01,3); // marginal widths on bounded scale
  // raws
  to_vector(I_raw) ~ std_normal();
  to_vector(J_raw) ~ std_normal();
  // hyper priors means
  mu_I    ~ std_normal();
  // hyper priors scales
  sigma_I ~ std_normal();
  sigma_J ~ std_normal();
  omega_beta ~ beta(2,2);
  
  //// Model ////
  {
    // parameters for MVN
    matrix[N,2] mu; 
    matrix[N,2] sigma;
    // Mean Vector  
    mu[ ,1]    = a_loc[ii] .* Tr_loc[jj] + b_loc[ii];
    mu[ ,2]    = Tr_wid[jj] + b_wid[ii]; 
    sigma[ ,1] = exp(log(a_loc[ii]) - log(E_loc[ii]) - log(lambda_loc[jj]));
    sigma[ ,2] = exp(               - log(E_wid[ii]) - log(lambda_wid[jj]));
   
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
  // Latent truth simplex
  matrix[J,3] Tr_splx;
  Tr_splx[,1] = Tr_splx_model[,1] .* 1.03 - 0.01;
  Tr_splx[,2] = Tr_splx_model[,2] .* 1.03 - 0.01;
  Tr_splx[,3] = Tr_splx_model[,3] .* 1.03 - 0.01;
  vector[J] Tr_loc_splx = Tr_splx[,1] + 0.5 .* Tr_splx[,2];
  vector[J] Tr_wid_splx = Tr_splx[,2];
}
