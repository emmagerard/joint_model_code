#######################################################################################################
#######################################################################################################
######################################Functions to simulate Rmax#######################################
#######################################################################################################
#######################################################################################################

library(parallel)
library(foreach)
library(doParallel)


simu_Rmax_population=function(path,res_save=TRUE){
  source(paste0(path,"/PK_PD_param.R"))
  
  t=seq(0,650,by=0.1)
  
  PK_param_fixed=list(mu=mu_PK,w=matrix(0,length(mu_PK),length(mu_PK)))
  PD_param_fixed=list(mu=mu_PD,w=matrix(0,length(mu_PD),length(mu_PD)))
  
  pop_profile=PK_PD_simu(t=t,doses=DL,t_doses=t_doses,t_inf=t_inf,PK_model=PK_model,PK_param_pop=PK_param_fixed,
                         PD_model=PD_model,PD_param_pop=PD_param_fixed,
                         PD_initial_values=PD_initial_values,method_ode=method_ode)
  
  Rmax_local_population=apply(matrix(pop_profile$R,ncol=nrow(DL)),2,extract_max_local,t,doses=rep(1,length(t_doses)),t_doses)
  
  Rmax_local_population[t(DL==0)]=NA
  
  if(res_save){
    save(Rmax_local_population, file = paste0(path,"/scenario/Rmax/Rmax_local_population.Rdata"))
  }
  
  return(Rmax_local_population)
}


simu_Rmax=function(path,n,Nsim,no_cores,parallel=TRUE,res_save=TRUE){
  
  setwd(path)
  
  source(paste0(path,"/PK_PD_param.R"))
  
  t=seq(0,650,by=0.1)
  
  vect_packages=c("deSolve","MASS")
  
  N_pack=n/Nsim
  
  if(parallel){
    #####Parallel
    
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    
    Rmax <- foreach(ip=1:N_pack,.combine=rbind,.packages = vect_packages,.export = ls(globalenv())) %dopar% {
      res=matrix(NA,nrow=Nsim*length(t_doses),ncol=nrow(DL)+2)
      colnames(res)=c("id","admin",paste0("DL",1:nrow(DL)))
      
      for(i in 1:Nsim){
        
        set.seed(seed+(ip-1)*Nsim+i)
        res[(i-1)*length(t_doses)+1:length(t_doses),]=cbind(rep((ip-1)*Nsim+i,length(t_doses)),1:length(t_doses),
                                                            PK_PD_simu_Rmax_local(t=t,
                                                                                  doses=DL,
                                                                                  t_doses=t_doses,
                                                                                  t_inf=t_inf,
                                                                                  PK_model=PK_model,
                                                                                  PK_param_pop=PK_param_pop,
                                                                                  PD_model=PD_model,
                                                                                  PD_param_pop=PD_param_pop,
                                                                                  PD_initial_values=PD_initial_values))
        
      } ## End loop i (simulations)
      
      res
      
    } ## Foreach loop (parallel computing)
    
    ### Release cluster
    stopCluster(cl)
  } else{
    
    Rmax <- foreach(ip=1:N_pack,.combine=rbind,.inorder=FALSE,.packages = vect_packages) %do% {  
      res=matrix(NA,nrow=Nsim*length(t_doses),ncol=nrow(DL)+2)
      colnames(res)=c("id","admin",paste0("DL",1:nrow(DL)))
      
      for(i in 1:Nsim){
        
        cat(ip,"/",N_pack,"-",i,"/",Nsim,"\n")
        
        set.seed(seed+(ip-1)*Nsim+i)
        res[(i-1)*length(t_doses)+1:length(t_doses),]=cbind(rep((ip-1)*Nsim+i,length(t_doses)),1:length(t_doses),
                                                            PK_PD_simu_Rmax_local(t=t,
                                                                                  doses=DL,
                                                                                  t_doses=t_doses,
                                                                                  t_inf=t_inf,
                                                                                  PK_model=PK_model,
                                                                                  PK_param_pop=PK_param_pop,
                                                                                  PD_model=PD_model,
                                                                                  PD_param_pop=PD_param_pop,
                                                                                  PD_initial_values=PD_initial_values))
        
      } ## End loop i (simulations)
      
      res
      
    } ## Foreach loop (parallel computing)
  }
  
  if(res_save){
    save(Rmax, file = paste0(path,"/scenario/Rmax/Rmax_local.Rdata"))
    return("OK")
  } else{
    return(Rmax)
  }
  
}
