########################################################################################################
#########################################Simulations parameters#########################################
########################################################################################################

path_simus=""

library(MASS)
library(dfcrm)

source(paste0(path_simus,"/tools/stan_functions.R"))

load(paste0(path_simus,"/scenario/scenario.Rdata"))
load(paste0(path_simus,"/scenario/pTox.Rdata"))
load(paste0(path_simus,"/scenario/Rmax/Rmax_local_population.Rdata"))
Rmax_global_population=apply(Rmax_local_population,2,max,na.rm=T)


###############################################Parameters###############################################

target=0.3
dose_ref=4
dose_50=6
wm=getprior(0.05,target,dose_ref,ncol(pTox))

# n_chains=4
# n_warmup=2500
# n_iter=5000
# n_PKPD=600

n_chains=4
n_warmup=1000
n_iter=2000
n_PKPD=200

###########################################Priors definition############################################

#############
###Approach 2
#############

####ESS = 1.6
prior_beta0_approach_2=c(
  logit(target),
  2
)
alpha_approach_2=5

beta_1_approach_2_function=function(beta0,p,logR,range=c(0,50)){
  optimize(optim_beta2,range,beta0=beta0,
           p=p,
           logR=logR)$minimum
}

beta_1_approach_2_pop=NULL

for(i in as.numeric(rownames(scenario[scenario$w_alpha==unique(scenario$w_alpha)[1],]))){

  beta_1_approach_2_pop=c(beta_1_approach_2_pop,
                          beta_1_approach_2_function(
                            beta0=prior_beta0_approach_2[1],
                            p=wm[(dose_ref-1):(dose_ref+1)],
                            logR=log(Rmax_global_population[as.numeric(scenario[i,2+1:ncol(pTox)])][(dose_ref-1):(dose_ref+1)]/
                                       Rmax_global_population[as.numeric(scenario[i,2+1:ncol(pTox)])][dose_ref]),
                            range=c(0,50)))
}


###ESS calculation
# N_ESS=1000000
# 
# beta0=rnorm(N_ESS,mean=prior_beta0_approach_2[1],sd=prior_beta0_approach_2[2])
# 
# res_ESS_logistic=NULL
# 
# for(j in as.numeric(rownames(scenario[scenario$w_alpha==unique(scenario$w_alpha)[1],]))){
#   
#   prior_beta1_approach_2_pop=c(alpha_approach_2,alpha_approach_2/beta_1_approach_2_pop[j])
#   
#   beta1=rgamma(N_ESS,shape=prior_beta1_approach_2_pop[1],rate=prior_beta1_approach_2_pop[2])
#   
#   p_ESS_2=matrix(NA,ncol=ncol(pTox),nrow=N_ESS)
#   
#   for(i in 1:ncol(pTox)){
#     p_ESS_2[,i]=Pi2(beta0,beta1,
#                     log(Rmax_global_population[as.numeric(scenario[j,2+1:ncol(pTox)])][i]/
#                           Rmax_global_population[as.numeric(scenario[j,2+1:ncol(pTox)])][dose_ref]))
#   }
#   
#   res_ESS_logistic=c(res_ESS_logistic,
#                      mean(apply(p_ESS_2,2,mean)*(1-apply(p_ESS_2,2,mean))/apply(p_ESS_2,2,var)-1))
#   
#   
# }
# mean(res_ESS_logistic)
# mean(res_ESS_logistic[1:4])
# 
# rm(p_ESS_2,beta1,prior_beta1_approach_2_pop,beta_1_approach_2_pop,beta0)


#############
###Approach 3
#############

####ESS = 1.8
sigma_mu=1
sigma_tau=1

###ESS calculation
# mu=rnorm(N_ESS,0,sigma_mu)
# tau=abs(rcauchy(N_ESS,0,sigma_tau))
# 
# res_ESS_hierarchical=NULL
# 
# for(j in as.numeric(rownames(scenario[scenario$w_alpha==unique(scenario$w_alpha)[1],]))){
#   
#   p_ESS_3=matrix(NA,ncol=ncol(pTox),nrow=N_ESS)
#   
#   for(i in 1:ncol(pTox)){
#     p_ESS_3[,i]=Pi3(mu,tau,
#                     log(Rmax_global_population[as.numeric(scenario[j,2+1:ncol(pTox)])][i]/
#                           Rmax_global_population[as.numeric(scenario[j,2+1:ncol(pTox)])][dose_50]))
#   }
#   
#   res_ESS_hierarchical=c(res_ESS_hierarchical,
#                          mean(apply(p_ESS_3,2,mean)*(1-apply(p_ESS_3,2,mean))/apply(p_ESS_3,2,var)-1))
#   
#   
# }
# mean(res_ESS_hierarchical)
# mean(res_ESS_hierarchical[1:4])
# rm(mu,tau,p_ESS_3,i,j,Rmax_local_population,Rmax_global_population)
