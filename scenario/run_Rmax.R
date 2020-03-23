#########################################################################################################
#############################################Rmax simulation#############################################
#########################################################################################################

path_simus=""

source(paste0(path_simus,"/tools/simu_Rmax.R"))

####Population profiles
simu_Rmax_population(path=path_simus,res_save=TRUE)

####Individual profiles
simu_Rmax(path=path_simus,n=100,Nsim=2,no_cores=5,parallel=TRUE,res_save=TRUE)


