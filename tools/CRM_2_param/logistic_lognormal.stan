data {
  int n_pat;
  int nb_doses;
  int y[n_pat];
  int dose_donnees[n_pat];
  vector[nb_doses] wm;
}
transformed data {
  vector[nb_doses] doses_trans;
  doses_trans = log(wm ./ (1-wm));
}
parameters {
  real log_beta0;
  real log_beta1;
}
transformed parameters { 
  vector[nb_doses] ptox;
  ptox = exp(log_beta0+exp(log_beta1)*doses_trans) ./ (1+exp(log_beta0+exp(log_beta1)*doses_trans));
}
model {
  log_beta0 ~ normal(0, 1);
  log_beta1 ~ normal(-0.5, 1);
  for(i in 1:n_pat){
    y[i] ~ bernoulli(ptox[dose_donnees[i]]);
  }
}

