# Implementation of the method to determine the maximal tolerated dose-regimen using PK/PD

This is the code for the paper "Bayesian dose-regimen assessment in early phase oncology incorporating pharmacokinetics and pharmacodynamics" authored by Emma Gerard, Sarah Zohar, Hoai-Thu Thai, Christelle Lorenzato, Marie-Karelle Riviere and Moreno Ursino. 

## Prerequisite

PK/PD estimation is performed via `R` using the software Monolix (http://lixoft.com/products/monolix/) through the package `lixoftConnectors`.

Bayesian analysis is performed with Stan via the package `rstan` (https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).


## Implementation

### Tools folder

The tools folder contains all generic functions used for the simulations:
  - `pk_functions.R` contains all PK/PD functions and functions to generate PK/PD data  
  - `scenario_functions.R` contains the functions to define the toxicity scenarios
  - `simu_Rmax.R` contains the functions to simulate the PD response for a large number of patients to define the toxicity scenarios
  - `simu_trials_PKPD.R` contains the functions to simulate the data under either the 3+3 or the CRM design and estimate the PKPD models with Monolix
  - `simu_trials_stat.R` contains the functions to extract PK/PD estimates and fit the statistical models on each dataset to estimate the probability of toxicity at each dose-regimen
  - `stan_functions.R` contains the functions to initialize parameters for Stan and other useful generic functions
  
### Scenario folder

The scenario folder contains programs to define the toxicity scenarios:
  1. `run_Rmax.R`: Simulation of the PD response of a large number of patients (using `PK_PD_param.R`). It generates 2 files in the scenario/Rmax folder: the population Rmax `Rmax_local_population.Rdata` and the individual Rmax `Rmax_local.Rdata`
  2. `run_scenario.R`: Definition of the toxicity scenarios from the data generated with `run_Rmax.R`. It generates 2 files in the scenario folder: `scenario.Rdata` with the value of the threshold and toxicity variability and `pTox.Rdata` with the true toxicity probabilities.
  
### Run

We will develop the steps to run the simulations:

  1. Define the PK/PD parameters in `PK_PD_param.R`
  2. Define the parameters for the simulations in `simu_param.R`
  3. Run the PK/PD estimation with `run_simu_trials_PKPD.R`. For each trial, data and estimated PK/PD parameters are written in different folders.
  4. Run the statistical estimation with `run_simu_trials_stat.R`. Results are written in a `.Rdata` file in the results folder of the design x scenario chosen.
  
## Disclaimer

Programs are provided "as is" without any warranty.
