
path_simus=""

dir_trials=paste0(path_simus,"/crm_trials_sc1")

N_sim=1

no_cores=5

load(paste0(dir_trials,"/trial_param.RData"))
save_output=paste0(design,"_sc",index_scenario,"_test")
#####


source(paste0(path_simus,"/tools/simu_trials_stat.R"))
source(paste0(path_simus,"/PK_PD_param.R"))
source(paste0(path_simus,"/simu_param.R"))


estim_trials(n_trials,scenario=index_scenario,
             dir_trials=dir_trials,
             names_PK=names(mu_PK),
             names_PD=names(mu_PD),
             fixed_param=fixed_param,
             t=t,
             DL=DL[as.numeric(scenario[index_scenario,3:ncol(scenario)]),],
             t_doses=t_doses,
             t_inf=t_inf[as.numeric(scenario[index_scenario,3:ncol(scenario)]),],
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
             t_inf_predict=t_inf_predict,
             N_sim=N_sim,
             no_cores=no_cores,
             save_output=save_output,
             start=0)

