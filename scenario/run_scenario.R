########################################################################################################
##########################################Simulation scenarios##########################################
########################################################################################################
path_simus=""

source(paste0(path_simus,"/tools/scenario_functions.R"))

path=paste0(path_simus,"/scenario")

load(paste0(path,"/Rmax/Rmax_local.Rdata"))
load(paste0(path,"/Rmax/Rmax_local_population.Rdata"))

Rmax_global=aggregate(Rmax[,grep("DL",colnames(Rmax))],by=list(Rmax[,"id"]), max,na.rm=T)[,-1]

w_alpha=c(0.25,0.5)

scenario_doses=scenario_function(w_alpha=w_alpha,
                                 target=c(rep(0.3,3),0.32,0.3),
                                 index_target=c(5,5,9,8,7),
                                 Rmax_global=Rmax_global,
                                 interval=c(0,5000))


index_doses=data.frame(matrix(0,nrow=nrow(scenario_doses),ncol=6))
names(index_doses)=paste0("s",1:6)

index_doses[(1:length(w_alpha))*5-5+1,]=rep(c(1:3,5:7),each=length(w_alpha))
index_doses[(1:length(w_alpha))*5-5+2,]=rep(c(3,5:8,10),each=length(w_alpha))
index_doses[(1:length(w_alpha))*5-5+3,]=rep(c(1:2,4:5,7,9),each=length(w_alpha))
index_doses[(1:length(w_alpha))*5-5+4,]=rep(c(2,4:5,7:9),each=length(w_alpha))
index_doses[(1:length(w_alpha))*5-5+5,]=rep(c(2,4:5,8:10),each=length(w_alpha))

scenario=cbind(scenario_doses,index_doses)
pTox=t(apply(scenario,1,pTox_calculation,Rmax_global))
colnames(pTox)=paste0("s",1:ncol(pTox))

save(scenario, file = paste0(path,"/scenario.Rdata"))
save(pTox, file = paste0(path,"/pTox.Rdata"))

