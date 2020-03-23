######
path_simus=""

design="crm"
tot_index_scenario=1
#tot_index_scenario=1:5

n_trials=5
no_cores=5
N_sim=1

######
source(paste0(path_simus,"/tools/simu_trials_PKPD.R"))
source(paste0(path_simus,"/PK_PD_param.R"))

dir=path_simus

# n_chains=3
# n_smooth=200
# n_explo=500

n_chains=1
n_smooth=50
n_explo=50

target=0.3
wm=getprior(0.05,target,4,6) #According to sc1
coh_size=3
n_tot=30

seed=25011995

create_mlxtran=create_mlxtran_Chen
create_model=create_model_Chen

load(paste0(dir,"/scenario/scenario.Rdata"))


for(index_scenario in tot_index_scenario){
  
  if((index_scenario %% 5)==0){
    predict=T
    DL_predict=DL[7,]
    t_inf_predict=t_inf[7,]
    
  } else{
    predict=F
    DL_predict=NULL
    t_inf_predict=NULL
  }
  
  dir.create(paste0(dir,"/",design,"_trials_sc",index_scenario))
  
  save(design,index_scenario,n_trials,target,seed,
       wm, coh_size,n_tot,
       predict,DL_predict,t_inf_predict,
       file = paste0(dir,"/",design,"_trials_sc",index_scenario,"/trial_param.RData"))
  
  simu_trials_PKPD_parallel(n_trials,
                            DL=DL[as.numeric(scenario[index_scenario,3:ncol(scenario)]),],
                            t,t_doses,
                            t_inf=t_inf[as.numeric(scenario[index_scenario,3:ncol(scenario)]),],
                            PK_model,PK_param_pop,
                            PD_model,PD_param_pop,
                            fixed_param,param_var,init,
                            b_PK,b_PD,PD_initial_values,method_ode,
                            n_chains,n_smooth,n_explo,
                            create_model,
                            create_mlxtran,
                            threshold=scenario$threshold[index_scenario],
                            w_alpha=scenario$w_alpha[index_scenario],
                            design,
                            wm,
                            target,
                            coh_size,n_tot,
                            dir=paste0(dir,"/",design,"_trials_sc",index_scenario),
                            seed,
                            no_cores,
                            N_sim)
}





