source(paste0(path_simus,"/tools/pk_functions.R"))
source(paste0(path_simus,"/tools/stan_models.R"))
source(paste0(path_simus,"/tools/stan_functions.R"))

library(parallel)
library(foreach)
library(doParallel)

return_PKPD_param=function(dir,
                           names_PK,
                           names_PD,
                           fixed_param){
  
  setwd(dir = dir)
  
  pop_param=read.table(file="results/populationParameters.txt",
                       sep=",",header=T)
  
  ind_param=read.table(file="results/IndividualParameters/estimatedIndividualParameters.txt",
                       sep=",",header=T)
  
  #############Population parameters#############
  
  ####Fixed effects vector####
  mu_PK=numeric(length(names_PK))
  names(mu_PK)=names_PK
  
  mu_PK_index=unlist(sapply(
    1:length(names_PK),
    function(x) grep(paste0(names_PK[x],"_pop"),pop_param$parameter)
  ))
  mu_PK_index2=unlist(sapply(
    1:length(as.character(pop_param$parameter[mu_PK_index])),
    function(x) grep(as.character(pop_param$parameter)[mu_PK_index[x]],paste0(names(mu_PK),"_pop"))
  ))
  
  mu_PK[mu_PK_index2]=pop_param$value[mu_PK_index]
  
  
  mu_PD=numeric(length(names_PD))
  names(mu_PD)=names_PD
  
  mu_PD_index=unlist(sapply(
    1:length(names_PD),
    function(x) grep(paste0(names_PD[x],"_pop"),pop_param$parameter)
  ))
  mu_PD_index2=unlist(sapply(
    1:length(as.character(pop_param$parameter[mu_PD_index])),
    function(x) grep(as.character(pop_param$parameter)[mu_PD_index[x]],paste0(names(mu_PD),"_pop"))
  ))
  
  mu_PD[mu_PD_index2]=pop_param$value[mu_PD_index]
  
  
  mu_PK_index_fixed=unlist(sapply(
    1:length(names_PK),
    function(x) grep(names_PK[x],names(fixed_param))
  ))
  
  if(length(mu_PK_index_fixed)>0){
    mu_PK_index_fixed2=unlist(sapply(
      1:length(fixed_param[mu_PK_index_fixed]),
      function(x) grep(names(fixed_param[mu_PK_index_fixed[x]]),names(mu_PK))
    ))
    
    mu_PK[mu_PK_index_fixed2]=fixed_param[mu_PK_index_fixed]
  }
  
  
  mu_PD_index_fixed=unlist(sapply(
    1:length(names_PD),
    function(x) grep(names_PD[x],names(fixed_param))
  ))
  
  if(length(mu_PD_index_fixed)>0){
    mu_PD_index_fixed2=unlist(sapply(
      1:length(fixed_param[mu_PD_index_fixed]),
      function(x) grep(names(fixed_param[mu_PD_index_fixed[x]]),names(mu_PD))
    ))
    
    mu_PD[mu_PD_index_fixed2]=fixed_param[mu_PD_index_fixed]
  }
  
  
  
  ####Variance covariance matrix####
  w_PK=matrix(0,length(mu_PK),length(mu_PK))
  rownames(w_PK)=colnames(w_PK)=names_PK
  
  w_PK_index=unlist(sapply(
    1:length(names_PK),
    function(x) grep(paste0("omega_",names_PK[x]),pop_param$parameter)
  ))
  w_PK_index2=unlist(sapply(
    1:length(as.character(pop_param$parameter[w_PK_index])),
    function(x) grep(as.character(pop_param$parameter)[w_PK_index[x]],paste0("omega_",names(mu_PK)))
  ))
  
  diag(w_PK)[w_PK_index2]=pop_param$value[w_PK_index]^2
  
  
  
  w_PD=matrix(0,length(mu_PD),length(mu_PD))
  rownames(w_PD)=colnames(w_PD)=names_PD
  
  w_PD_index=unlist(sapply(
    1:length(names_PD),
    function(x) grep(paste0("omega_",names_PD[x]),pop_param$parameter)
  ))
  w_PD_index2=unlist(sapply(
    1:length(as.character(pop_param$parameter[w_PD_index])),
    function(x) grep(as.character(pop_param$parameter)[w_PD_index[x]],paste0("omega_",names(mu_PD)))
  ))
  
  diag(w_PD)[w_PD_index2]=pop_param$value[w_PD_index]^2
  
  
  PK_param_pop=list(mu=mu_PK,
                    w=w_PK)
  PD_param_pop=list(mu=mu_PD,
                    w=w_PD)
  
  
  #############Individual parameters#############
  PK_param_ind=matrix(0,nrow=nrow(ind_param),ncol=length(names_PK)+1)
  
  colnames(PK_param_ind)=c("id",names_PK)
  PK_param_ind[,1]=ind_param$id
  
  ind_PK_index=unlist(sapply(
    1:length(colnames(PK_param_ind)),
    function(x) grep(paste0(colnames(PK_param_ind)[x],"_SAEM"),colnames(ind_param))
  ))
  
  ind_PK_index2=unlist(sapply(
    1:length(colnames(ind_param)[ind_PK_index]),
    function(x) grep(colnames(ind_param)[ind_PK_index[x]],paste0(colnames(PK_param_ind),"_SAEM"))
  ))
  
  PK_param_ind[,ind_PK_index2]=as.matrix(ind_param[,ind_PK_index])
  
  
  
  PD_param_ind=matrix(0,nrow=nrow(ind_param),ncol=length(names_PD)+1)
  
  colnames(PD_param_ind)=c("id",names_PD)
  PD_param_ind[,1]=ind_param$id
  
  ind_PD_index=unlist(sapply(
    1:length(colnames(PD_param_ind)),
    function(x) grep(paste0(colnames(PD_param_ind)[x],"_SAEM"),colnames(ind_param))
  ))
  
  ind_PD_index2=unlist(sapply(
    1:length(colnames(ind_param)[ind_PD_index]),
    function(x) grep(colnames(ind_param)[ind_PD_index[x]],paste0(colnames(PD_param_ind),"_SAEM"))
  ))
  
  PD_param_ind[,ind_PD_index2]=as.matrix(ind_param[,ind_PD_index])
  
  
  if(length(mu_PK_index_fixed)>0){
    
    PK_param_ind[,1+mu_PK_index_fixed2]=matrix(rep(fixed_param[mu_PK_index_fixed],each=nrow(PK_param_ind)),nrow=nrow(PK_param_ind))
    
  }
  
  if(length(mu_PD_index_fixed)>0){
    
    PD_param_ind[,1+mu_PD_index_fixed2]=matrix(rep(fixed_param[mu_PD_index_fixed],each=nrow(PD_param_ind)),nrow=nrow(PD_param_ind))
    
  }
  
  return(list(
    PK_param_pop=PK_param_pop,
    PD_param_pop=PD_param_pop,
    PK_param_ind=PK_param_ind,
    PD_param_ind=PD_param_ind
  ))
  
  
  
}


return_data=function(dir,
                     t,
                     DL,
                     t_doses,
                     t_inf,
                     PK_model,
                     PK_param_ind,
                     PD_model,
                     PD_param_ind,
                     PD_initial_values,
                     method_ode,
                     predict_tox=F){
  
  data=read.table(file=paste0(dir,"/data.txt"),header=T)
  
  max_admin=sapply(unique(data$id),function(x) sum(data$id==x))
  
  data$admin=as.numeric(unlist(sapply(max_admin,function(x) 1:x)))
  
  data$R=0
  
  if(predict_tox){
    data_predict=data.frame(
      id=rep(unique(data$id),each=length(t_doses)),
      tox=NA,
      dose=rep(data[!duplicated(data[,c(1,3)]),"dose"],each=length(t_doses)),
      admin=rep(1:length(t_doses),length(unique(data$id)))
    )
    
    data_merge=merge(data_predict[,c("id","admin")],data[c("id","tox","admin")],
                     by=c("id","admin"),all.x=T)
    
    data_merge$tox[is.na(data_merge$tox)]=1
    
    data_predict$tox=data_merge$tox
  }
  
  for(i in 1:nrow(PK_param_ind)){
    
    if(predict_tox){
      Rmax_ind=PK_PD_simu_Rmax_local_ind(t=t,
                                         doses=DL[data[data$id==i,][1,"dose"],],
                                         t_doses=t_doses,
                                         t_inf=t_inf[data[data$id==i,][1,"dose"],],
                                         PK_model=PK_model,
                                         PK_param_ind=as.list(PK_param_ind[i,2:ncol(PK_param_ind)]),
                                         PD_model=PD_model,
                                         PD_param_ind=as.list(PD_param_ind[i,2:ncol(PD_param_ind)]),
                                         PD_initial_values=PD_initial_values,
                                         method_ode=method_ode)
      
      data_predict$R[data_predict$id==i]=as.numeric(Rmax_ind)
      
    } else{
      Rmax_ind=PK_PD_simu_Rmax_local_ind(t=t,
                                         doses=DL[data[data$id==i,][1,"dose"],1:sum(data$id==i)],
                                         t_doses=t_doses[1:sum(data$id==i)],
                                         t_inf=t_inf[data[data$id==i,][1,"dose"],1:sum(data$id==i)],
                                         PK_model=PK_model,
                                         PK_param_ind=as.list(PK_param_ind[i,2:ncol(PK_param_ind)]),
                                         PD_model=PD_model,
                                         PD_param_ind=as.list(PD_param_ind[i,2:ncol(PD_param_ind)]),
                                         PD_initial_values=PD_initial_values,
                                         method_ode=method_ode)
      
      data$R[data$id==i]=as.numeric(Rmax_ind)
    }
    
  }
  
  if(predict_tox){
    return(data=data_predict)
  } else{
    return(data=data)
  }
  
}




removeError=function(data){
  
  data_no_error=data
  
  data_index=NULL
  
  nb_error=0
  
  data$R_error_init=data$R
  
  
  for(index in unique(data$id)[tapply(data$tox, data$id, sum)==1]){
    data_index=data[data$id==index,]
    
    rownames_pbm=rownames(data_index[data_index$R>data_index$R[data_index$tox==1],])
    
    nb_error=nb_error+length(rownames_pbm)
    
    if(length(rownames_pbm)>0){
      data_no_error=data_no_error[!rownames(data_no_error) %in% rownames_pbm,]
      
      data_init=data_no_error[data_no_error$id==index & data_no_error$tox==0,]
      
      if(nrow(data_init)>0){
        data$R_error_init[rownames(data) %in% rownames_pbm]=max(data_init$R)
      } else{
        data$R_error_init[rownames(data) %in% rownames_pbm]=0.1
      }
      
    }
    
    
  }
  
  return(list(data_no_error=data_no_error,prop_error=nb_error/nrow(data),R_error_init=data$R_error_init))
}




estim_one_trial=function(itrial,scenario,
                         dir,names_PK,
                         names_PD,
                         fixed_param,
                         t,
                         DL,
                         t_doses,
                         t_inf,
                         PK_model,
                         PD_model,
                         PD_initial_values,
                         method_ode,
                         design,
                         target,
                         dose_ref,
                         dose_50,
                         seed,
                         n_chains,
                         n_warmup,
                         n_iter,
                         n_PKPD,
                         prior_beta0_approach_2,
                         beta_1_approach_2_function,
                         alpha_approach_2,
                         wm,
                         sigma_mu,
                         sigma_tau,
                         predict,
                         DL_predict,
                         t_inf_predict){
  
  ######Results definition
  seed_trial=seed+itrial
  n_doses=nrow(DL)
  
  method=c("logistic","logistic_predict","hierarchical")
  
  if(design=="3p3"){
    method=c(method,"3p3")
  } 
  if(design=="crm"){
    method=c(method,"crm_logistic")
  }
  
  n_predict=ifelse(predict,length(DL_predict)/length(t_doses),0)
  
  res=data.frame(matrix(NA,
                        nrow=length(method),
                        ncol=2*n_doses+9+n_predict))
  if(predict){
    
    colnames(res)=c("trial","Scenario","n","design","Method","Seed",
                    paste0("pTox_DL",1:n_doses),
                    paste0("pTox_predict_DL",1:n_predict),
                    paste0("n_DL",1:n_doses),
                    "Rhat_issue","n_eff_issue","data_error")
    
  } else{

    colnames(res)=c("trial","Scenario","n","design","Method","Seed",
                    paste0("pTox_DL",1:n_doses),
                    paste0("n_DL",1:n_doses),
                    "Rhat_issue","n_eff_issue","data_error")
  }
  
  
  
  ######Data definition
  PKPD_param=return_PKPD_param(dir,names_PK,names_PD,fixed_param)
  
  PK_param_pop=PKPD_param$PK_param_pop
  PD_param_pop=PKPD_param$PD_param_pop
  PK_param_ind=PKPD_param$PK_param_ind
  PD_param_ind=PKPD_param$PD_param_ind
  
  
  data=return_data(dir,t,DL,t_doses,t_inf,
                   PK_model,PK_param_ind,
                   PD_model,PD_param_ind,
                   PD_initial_values,
                   method_ode,predict_tox=F)
  
  
  
  data_predict=return_data(dir,t,DL,t_doses,t_inf,
                           PK_model,PK_param_ind,
                           PD_model,PD_param_ind,
                           PD_initial_values,
                           method_ode,predict_tox=T)
  
  error_removal=removeError(data)
  data_no_error=error_removal$data_no_error
  prop_error=error_removal$prop_error
  
  
  Rmax_pop=PK_PD_simu_Rmax(t,doses=DL,t_doses,t_inf,
                           PK_model,
                           PK_param_pop=list(mu=PK_param_pop$mu,
                                             w=matrix(0,length(PK_param_pop$mu),length(PK_param_pop$mu))),
                           PD_model,
                           PD_param_pop=list(mu=PD_param_pop$mu,
                                             w=matrix(0,length(PD_param_pop$mu),length(PD_param_pop$mu))),
                           PD_initial_values,
                           method_ode=method_ode)
  
  data_aggreg=aggregate(data,by=list(data$id),max)[,-1]
  
  data_aggreg_predict=aggregate(data_predict,by=list(data_predict$id),max)[,-1]
  
  
  beta_1_approach_2=beta_1_approach_2_function(
    beta0=prior_beta0_approach_2[1],
    p=wm[(dose_ref-1):(dose_ref+1)],
    logR=log(Rmax_pop[(dose_ref-1):(dose_ref+1)]/Rmax_pop[dose_ref]),
    range=c(0,50))
  
  prior_beta1_approach_2=c(alpha_approach_2,alpha_approach_2/beta_1_approach_2)
  
  
  data_stan_logistic=list(n_patients=nrow(data_aggreg),
                          y=data_aggreg$tox,
                          logR=log(data_aggreg$R/Rmax_pop[dose_ref]),
                          prior_beta0=prior_beta0_approach_2,
                          prior_beta1=prior_beta1_approach_2)
  
  data_stan_logistic_predict=list(n_patients=nrow(data_aggreg_predict),
                                  y=data_aggreg_predict$tox,
                                  logR=log(data_aggreg_predict$R/Rmax_pop[dose_ref]),
                                  prior_beta0=prior_beta0_approach_2,
                                  prior_beta1=prior_beta1_approach_2)
  
  data_stan_hierarchical=list(n_patients=length(unique(data_no_error$id)),N=nrow(data_no_error),id=data_no_error$id,
                              y=data_no_error$tox,logR=log(data_no_error$R/Rmax_pop[dose_50]),
                              sigma_mu=sigma_mu,
                              sigma_tau=sigma_tau)
  

  set.seed(seed_trial)
  init_approach2=lapply(1:n_chains, function(id) initfunction_approach2(chain_id = id))
  init_approach3=lapply(1:n_chains, function(id) initfunction_approach3(chain_id = id))
  
  sampl_logistic=sampling(sm_approach_2, data = data_stan_logistic, 
                          iter=n_iter, warmup=n_warmup,
                          chains = n_chains,control = list(adapt_delta = 0.99),
                          cores=1,seed=seed_trial,init=init_approach2)
  
  sampl_logistic_predict=sampling(sm_approach_2, data = data_stan_logistic_predict, 
                                  iter=n_iter, warmup=n_warmup,
                                  chains = n_chains,control = list(adapt_delta = 0.99),
                                  cores=1,seed=seed_trial,init=init_approach2)
  
  sampl_hierarchical=sampling(sm_approach_3, data = data_stan_hierarchical, iter=n_iter, warmup=n_warmup,
                              chains = n_chains,control = list(adapt_delta = 0.99),
                              cores=1,seed=seed_trial,init=init_approach3)
  
  ##########
  
  if(predict){
    DL=rbind(DL,DL_predict)
    t_inf=rbind(t_inf,t_inf_predict)
  }
  
  
  Rmax_function=function(i){
    cbind(
      rep(i,length(t_doses)),
      1:length(t_doses),
      PK_PD_simu_Rmax_local(t=t,
                            doses=DL,
                            t_doses=t_doses,
                            t_inf=t_inf,
                            PK_model=PK_model,
                            PK_param_pop=PK_param_pop,
                            PD_model=PD_model,
                            PD_param_pop=PD_param_pop,
                            PD_initial_values=PD_initial_values))
  }
  
  set.seed(seed_trial)
  Rmax_MC_local=do.call(rbind,lapply(1:n_PKPD,Rmax_function))
  colnames(Rmax_MC_local)=c("id","admin",paste0("DL",1:nrow(DL)))
  
  
  Rmax_MC_global=aggregate(Rmax_MC_local[,grep("DL",colnames(Rmax_MC_local))],
                           by=list(Rmax_MC_local[,"id"]),
                           max,na.rm=T)[,-1]
  
  
  logR_MC_approach2=as.matrix(log(Rmax_MC_global/Rmax_pop[dose_ref]))
  logR_MC_approach3=as.matrix(log(Rmax_MC_global/Rmax_pop[dose_50])) 
  
  beta0_logistic=extract(sampl_logistic)$beta0
  beta1_logistic=extract(sampl_logistic)$beta1
  
  beta0_logistic_predict=extract(sampl_logistic_predict)$beta0
  beta1_logistic_predict=extract(sampl_logistic_predict)$beta1
  
  mu_hierarchical=extract(sampl_hierarchical)$mu
  tau_hierarchical=extract(sampl_hierarchical)$tau
  
  logR_MC_logistic_big=do.call(rbind,
                               lapply(1:length(beta0_logistic),function(i) logR_MC_approach2))
  
  logR_MC_hierarchical_big=do.call(rbind,
                                   lapply(1:length(mu_hierarchical),function(i) logR_MC_approach3))
  
  
  beta0_logistic_big=rep(beta0_logistic,each=n_PKPD)
  beta1_logistic_big=rep(beta1_logistic,each=n_PKPD)
  
  beta0_logistic_predict_big=rep(beta0_logistic_predict,each=n_PKPD)
  beta1_logistic_predict_big=rep(beta1_logistic_predict,each=n_PKPD)
  
  mu_hierarchical_big=rep(mu_hierarchical,each=n_PKPD)
  tau_hierarchical_big=rep(tau_hierarchical,each=n_PKPD)
  
  p_estim_logistic_big=matrix(sapply(1:length(beta0_logistic_big),function(x){
    Pi2(beta0=beta0_logistic_big[x],beta1=beta1_logistic_big[x],
        logR=logR_MC_logistic_big[x,])
  }),nrow=length(beta0_logistic_big),byrow = T)
  
  p_estim_logistic_predict_big=matrix(sapply(1:length(beta0_logistic_predict_big),function(x){
    Pi2(beta0=beta0_logistic_predict_big[x],beta1=beta1_logistic_predict_big[x],
        logR=logR_MC_logistic_big[x,])
  }),nrow=length(beta0_logistic_predict_big),byrow = T)
  
  p_estim_hierarchical_big=matrix(sapply(1:length(mu_hierarchical_big),function(x){
    Pi3(mu=mu_hierarchical_big[x],tau=tau_hierarchical_big[x],
        logR=logR_MC_hierarchical_big[x,])
  }),nrow=length(mu_hierarchical_big),byrow = T)
  
  ###############Results
  ##Logistic model
  res[1,]=c(itrial,scenario,length(unique(data$id)),design,method[1],seed_trial,
            colMeans(p_estim_logistic_big),
            as.numeric(table(factor(data_aggreg$dose,levels=1:n_doses))),
            as.numeric(sum(summary(sampl_logistic)$summary[,"Rhat"]>=1.1)>0),
            as.numeric(sum(summary(sampl_logistic)$summary[,"n_eff"]/(n_iter-n_warmup)<0.001)>0),
            NA)
  
  
  ##Logistic model (predict)
  res[2,]=c(itrial,scenario,length(unique(data$id)),design,method[2],seed_trial,
            colMeans(p_estim_logistic_predict_big),
            as.numeric(table(factor(data_aggreg_predict$dose,levels=1:n_doses))),
            as.numeric(sum(summary(sampl_logistic_predict)$summary[,"Rhat"]>=1.1)>0),
            as.numeric(sum(summary(sampl_logistic_predict)$summary[,"n_eff"]/(n_iter-n_warmup)<0.001)>0),
            NA)  
  
  ##Hierarchical model
  res[3,]=c(itrial,scenario,length(unique(data_no_error$id)),design,method[3],seed_trial,
            colMeans(p_estim_hierarchical_big),
            as.numeric(table(factor(data_aggreg$dose,levels=1:n_doses))),
            as.numeric(sum(summary(sampl_hierarchical)$summary[,"Rhat"]>=1.1)>0),
            as.numeric(sum(summary(sampl_hierarchical)$summary[,"n_eff"]/(n_iter-n_warmup)<0.001)>0),
            prop_error)
  
  ##3+3 
  if(design=="3p3"){
    res_algo=numeric(n_doses)
    res_algo[as.numeric(unlist(read.table(file=paste0(dir,"/MTD.txt"))))]=target
    res[4,]=c(itrial,scenario,length(unique(data$id)),design,method[5],seed_trial,
              res_algo,
              rep(NA,n_predict),
              as.numeric(table(factor(data_aggreg$dose,levels=1:n_doses))),
              NA,NA,NA)
  }
  
  ##CRM
  if(design=="crm"){
    
    res[4,]=c(itrial,scenario,length(unique(data$id)),design,method[5],seed_trial,
              as.numeric(unlist(read.table(file=paste0(dir,"/pTox.txt")))),
              rep(NA,n_predict),
              as.numeric(table(factor(data_aggreg$dose,levels=1:n_doses))),
              NA,NA,NA)
  }
  
  return(res)
}


estim_trials=function(n_trials,scenario,
                      dir_trials,names_PK,
                      names_PD,
                      fixed_param,
                      t,
                      DL,
                      t_doses,
                      t_inf,
                      PK_model,
                      PD_model,
                      PD_initial_values,
                      method_ode,
                      design,
                      target,
                      dose_ref,
                      dose_50,
                      seed,
                      n_chains,
                      n_warmup,
                      n_iter,
                      n_PKPD,
                      prior_beta0_approach_2,
                      beta_1_approach_2_function,
                      alpha_approach_2,
                      wm,
                      sigma_mu,
                      sigma_tau,
                      predict,
                      DL_predict,
                      t_inf_predict,
                      N_sim,
                      no_cores,
                      save_output,
                      start=0){
  
  
  setwd(dir_trials)
  
  vect_packages=c("deSolve","rstan","MASS","reshape2","dfcrm")
  
  N_pack=n_trials/N_sim
  
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  res_trials<- foreach(ip=1:N_pack,.combine=rbind,.inorder=FALSE,.packages = vect_packages,.export = ls(globalenv())) %dopar% { 
   
    res=NULL
    
    for(i in 1:N_sim){
      
      cat(ip,"/",N_pack,"-",i,"/",N_sim,"\n")
      
      res=rbind(
        res,
        estim_one_trial(itrial=(ip-1)*N_sim+i+start,
                        scenario=scenario,
                        dir=paste0(dir_trials,"/trial",(ip-1)*N_sim+i+start),
                        names_PK=names_PK,
                        names_PD=names_PD,
                        fixed_param=fixed_param,
                        t=t,
                        DL=DL,
                        t_doses=t_doses,
                        t_inf=t_inf,
                        PK_model=PK_model,
                        PD_model=PD_model,
                        PD_initial_values=PD_initial_values,
                        method_ode=method_ode,
                        design=design,
                        target=target,
                        dose_ref=dose_ref,
                        dose_50=dose_50,
                        seed=seed,
                        n_chains=n_chains,
                        n_warmup=n_warmup,
                        n_iter=n_iter,
                        n_PKPD=n_PKPD,
                        prior_beta0_approach_2=prior_beta0_approach_2,
                        beta_1_approach_2_function=beta_1_approach_2_function,
                        alpha_approach_2=alpha_approach_2,
                        wm=wm,
                        sigma_mu=sigma_mu,
                        sigma_tau=sigma_tau,
                        predict=predict,
                        DL_predict=DL_predict,
                        t_inf_predict=t_inf_predict)
      )
      
      
      
    } ## End loop i (simulations)
    
    setwd(dir_trials)
    
    res
    
  }
  
  stopCluster(cl)
  
  if(is.null(save_output)){
    return(res_trials)
  } else{
    
    setwd(dir_trials)
    
    if (!file.exists("results")){
      dir.create("results",showWarnings = FALSE)
    }
    
    save(res_trials,file=paste0("results/",save_output,".RData"))
    return("OK")
  }
  
}
