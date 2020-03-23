############################################################################
###############################Initial values###############################
############################################################################

initfunction_approach2=function(chain_id = 1) {
  list(
    beta0=rnorm(1,0,1),
    beta1=runif(1,0,25)
  )
} 

initfunction_approach3=function(chain_id = 1) {
  list(
    mu=rnorm(1,0,1),
    tau=rexp(1,5)
  )
}

################

logit=function(x) {
  log(x/(1-x))
}

invlogit=function (x) {
  1/(1 + exp(-x))
}

Pi2=function(beta0,beta1,logR){
  invlogit(beta0+beta1*logR)
}

Pi3=function(mu,tau,logR){
  pnorm(logR,mean=mu,sd=tau)
}

optim_beta2=function(beta1,beta0,p,logR){
  sum((p-Pi2(beta0=beta0,beta1=beta1,logR))^2)
}