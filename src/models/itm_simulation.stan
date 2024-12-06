//////////////////////////////////////////////////////////////////////////////// 
// Interval Truth Model
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
  matrix[I,5] I_raw;
    // hyperpriors person
  vector[2] mu_I;
  vector[5] sigma_I;
// item parameters
  // raws for non-centered parameterization
  matrix[J,5] J_raw;
  // hyperpriors item
  vector[3] mu_J;
  vector[5] sigma_J;
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
  vector[J] Tr_loc; // latent truth location
  vector[J] Tr_wid; // latent truth width
  vector[J] lambda_loc; // item difficulty / discernibility location
  vector[J] lambda_wid; // item difficulty / discernibility width
  // correlation residual
  vector[J] omega;

  // person parameters
  E_loc = exp(I_raw[,1] * exp(sigma_I[1]) + mu_I[1]);
  E_wid = exp(I_raw[,2] * exp(sigma_I[2]) + mu_I[2]);
  a_loc = exp(I_raw[,3] * exp(sigma_I[3]));
  b_loc = I_raw[,4] * exp(sigma_I[4]);
  b_wid = I_raw[,5] * exp(sigma_I[5]);
  //item parameters
  Tr_loc     =      J_raw[,1] * exp(sigma_J[1]) + mu_J[1];  // latent truth
  Tr_wid     =      J_raw[,2] * exp(sigma_J[2]) + mu_J[2];  // latent truth
  lambda_loc =  exp(J_raw[,3] * exp(sigma_J[3])); // item difficulty / discernibility
  lambda_wid =  exp(J_raw[,4] * exp(sigma_J[4])); // item difficulty / discernibility
  // correlation residual
  omega      = tanh(J_raw[,5] * exp(sigma_J[5]) + mu_J[3]);

}
//////////////////////////////////////////////////////////////////////////////// 
model{
  // raws
  to_vector(I_raw) ~ std_normal();
  to_vector(J_raw) ~ std_normal();
  // hyper priors means
  mu_I    ~ student_t(3,0,2);
  mu_J    ~ student_t(3,0,2);
  // hyper priors scales
  sigma_I ~ student_t(3,0,2);
  sigma_J ~ student_t(3,0,2);
// Model //
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
      matrix[2,2] L = cholesky_decompose(Omega);
    // Likelihood
    Y[n] ~ multi_normal_cholesky(mu[n,], diag_pre_multiply(sigma[n,], L));
    } //end n
  } // end block
}
////////////////////////////////////////////////////////////////////////////////
 generated quantities{
  // Latent truth simplex
  matrix[J,3] Tr_splx;
  {
  if(link == 1){  
    // ILR
    vector[J] Sum = exp(sqrt(2) .* Tr_loc) + 
                    exp(sqrt(3.0/2) .* Tr_wid + Tr_loc ./ sqrt(2)) + 
                    1;
    Tr_splx[ ,1] = exp(sqrt(2) .* Tr_loc)                         ./ Sum;
    Tr_splx[ ,2] = exp(sqrt(3.0/2) .* Tr_wid + Tr_loc ./ sqrt(2)) ./ Sum;
    Tr_splx[ ,3] = 1                                              ./ Sum;
    }else{
    // SB
    Tr_splx[ ,1] = exp(Tr_loc) ./ ((exp(Tr_loc) + 1) .* (exp(Tr_wid) + 1));
    Tr_splx[ ,2] = exp(Tr_wid) ./  (exp(Tr_wid) + 1);
    Tr_splx[ ,3] = 1 ./           ((exp(Tr_loc) + 1) .* (exp(Tr_wid) + 1));
    } // end if
  } // end block
}
