#######################################################################################################
#######################################################################################################
##############################################Stan models##############################################
#######################################################################################################
#######################################################################################################

library(rstan)

############################################################################
###################################Models###################################
############################################################################

approach_2="
  data {
    int<lower=0> n_patients;         // number of observations (patients)
    int<lower=0,upper=1> y[n_patients]; // binary toxicity response
    vector[n_patients] logR;         // Log tranformation of the PD response 
    real prior_beta0[2];    // normal prior
    real prior_beta1[2];    // gamma prior
  }
  parameters {
    real beta0;
    real<lower=0> beta1;
  }
  model {
    y ~ bernoulli_logit(beta0+beta1*logR);

    beta0 ~ normal(prior_beta0[1],prior_beta0[2]);
    beta1 ~ gamma(prior_beta1[1],prior_beta1[2]);

  }"

approach_3="
  data {
    int<lower=0> n_patients;         // number of patients
    int<lower=0> N;         // number of patients x admin
    int<lower=0,upper=n_patients> id[N];// patients' id
    int <lower=0,upper=1> y[N];// binary toxicity response
    real logR[N];
    real<lower=0> sigma_mu;
    real<lower=0> sigma_tau;
  }

  transformed data {
    vector[n_patients] lower_bound;
    vector[n_patients] upper_bound;

    for(i in 1:n_patients){
      lower_bound[i]=negative_infinity();
      upper_bound[i]=positive_infinity();
    }

    for(i in 1:N){
      if (y[i] == 0) {
        if(logR[i]>lower_bound[id[i]]){
          lower_bound[id[i]]=logR[i];
        }
      } else {
        if(logR[i]<upper_bound[id[i]]){
          upper_bound[id[i]]=logR[i];
        }
      }
    }
  }

  parameters {
    real mu;
    real<lower=0> tau;
    vector[n_patients] z_raw;
  }

  transformed parameters {
    vector[n_patients] z;
    vector[n_patients] log_jacobian;
  
    for(i in 1:n_patients) {
      if(is_inf(lower_bound[i])) {
        if(is_inf(upper_bound[i])) {
          //No bounds, direct copy
          z[i] = z_raw[i];
          log_jacobian[i] = 0;
        } else {
          //upper bound only
          z[i] = upper_bound[i] - exp(z_raw[i]);
          log_jacobian[i] = z_raw[i];
        }
      } else {
        if(is_inf(upper_bound[i])) {
          //Lower bound only
          z[i] = lower_bound[i] + exp(z_raw[i]);
          log_jacobian[i] = z_raw[i];
        } else {
          //Both bounds
          real inv_logit_z = inv_logit(z_raw[i]);
          z[i] = lower_bound[i] + (upper_bound[i] - lower_bound[i]) * inv_logit_z;
          log_jacobian[i] = log(upper_bound[i] - lower_bound[i])  + log(inv_logit_z) + log(1 - inv_logit_z);
        }
      }
    }
  }

  model {
    z ~ normal(mu,tau);
    target += sum(log_jacobian); //Add the log Jacobian to target density
    mu ~ normal(0,sigma_mu);
    tau ~ cauchy(0,sigma_tau);
  }"



############################################################################
##################################Sampling##################################
############################################################################

sm_approach_2=stan_model(model_code = approach_2,auto_write = T)
sm_approach_3=stan_model(model_code = approach_3,auto_write = T)



