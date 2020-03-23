#######################################################################################################
#######################################################################################################
###############################Functions to obtain simulation scenarios################################
#######################################################################################################
#######################################################################################################


pTox_calculation=function(scenario,Rmax){
  res=NULL
  threshold=scenario[1]
  w_alpha=scenario[2]
  ####
  if(length(scenario)>2){
    id_doses=as.numeric(scenario[3:length(scenario)])
  } else{
    id_doses=1:ncol(Rmax)
  }
  
  Rmax=Rmax[,id_doses]
  names(Rmax)=paste0("DL",1:ncol(Rmax))
  
  if(w_alpha==0){
    res=colMeans(Rmax>=threshold,na.rm=T)
  } else{
    res=colMeans(1-pnorm((log(threshold)-log(as.matrix(Rmax)))/w_alpha),na.rm=T)
  }
  return(res)
}


optim_function=function(threshold,w_alpha,Rmax_global,target,index_target){
  if(w_alpha==0){
    res=mean(Rmax_global[,index_target]>=threshold,na.rm=T)
  } else{
    res=mean(1-pnorm((log(threshold)-log(Rmax_global[,index_target]))/w_alpha),na.rm=T)
  }
  
  return(abs(res-target))
}

scenario_function=function(w_alpha,
                           target,
                           index_target,
                           Rmax_global,
                           interval){
  
  res=data.frame(
    threshold=numeric(length(w_alpha)*length(target)),
    w_alpha=rep(w_alpha,each=length(target))
  )
  
  target_tot=rep(target,length(w_alpha))
  index_target_tot=rep(index_target,length(w_alpha))
  
  for(i in 1:nrow(res)){
    res$threshold[i]=optimize(optim_function,interval,w_alpha=res$w_alpha[i],
                              Rmax_global=Rmax_global,target=target_tot[i],index_target=index_target_tot[i])$minimum
  }
  
  return(res)
  
}
