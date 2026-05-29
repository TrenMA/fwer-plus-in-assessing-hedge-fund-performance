rm(list = ls())

library(qvalue)
library(lmtest)
library(sandwich)
library(dplyr)

code_folder <- "input the path to your 'Simulation part' folder here"

folder_simulation_results <- paste0(code_folder,'/Results')

setwd(code_folder)

source('FWER Procedure.R')
source('functions_for_simulations.R')

# return data

df_ret <- readRDS('Unsmoothed Return for Simulation.RDS')

if (!dir.exists(folder_simulation_results)) {
  dir.create(folder_simulation_results)
} 

################################################################################
# Performance comparisons among fwer+, stepM and StepSPA with varying alpha
################################################################################

b_list <- list(c(2,1,1),c(2.4,1.5,1.5),c(2.8,2,2))
# cases of b_list to run
B <- c(1,2,3)[2]

# number of funds in population:
N_funds <- c(500,1000,2000,3000)[2]

# cases of number of observations per fund to run: 
n_obs <- c(24,36,48)[2]

# alpha cases to run
fund_alpha <- c(0.5,1,1.5)

# Dependence types of covariates to run: 
Z_type_set <- c('independent','dependent','noniid')[1]

# if there is no noise in covariate consider, set include_noise <- NULL; 
# if yes, set include_noise <- '_Noise_Covariates'
include_noise <- NULL
#include_noise <- '_Noise_Covariates'

# If running StepM and StepSPA, change FALSE below to TRUE 
include_stepM <- include_spa <- TRUE;  

q_stepM <- q_spa <- 1/6; 

# and choose number of cores for parallel (mc.cores) to speed up the process: 
# if you are running in windows operation system, you need to set mc.cores <- 1
# otherwise you can set mc.cores <- parallel::detectCores() - 2 
mc.cores <- 32

source('Simulation Core.R')  # FYI: for mc.cores = 32, it takes around 30 hours for this line to complete.

source('Plot - Performance comparison among fwer plus StepM and StepSPA under varying alpha.R')

setwd(code_folder)

################################################################################
# Performance of fwer+ under varying signal settings.
################################################################################
# covariate coefficient cases:
b_list <- list(c(2,1,1),c(2.4,1.5,1.5),c(2.8,2,2))
# cases of b_list to run
B <- c(1,2,3)

# number of funds in population:
N_funds <- c(500,1000,2000,3000)[2]

# cases of number of observations per fund to run: 
n_obs <- c(24,36,48)[2]

# alpha cases to run
fund_alpha <- c(0.5,1,1.5)

# Dependence types of covariates to run: 
Z_type_set <- c('independent','dependent','noniid')[1:2]

# if there is no noise in covariate consider, set include_noise <- NULL; 
# if yes, set include_noise <- '_Noise_Covariates'
include_noise <- NULL
#include_noise <- '_Noise_Covariates'

# If running StepM and StepSPA, change FALSE below to TRUE 
include_stepM <- include_spa <- FALSE; 

source('Simulation Core.R') 

source('Plot - Performance under varying signal settings.R')

setwd(code_folder)

################################################################################
# Performance of fwer+ under use of insufficient information
################################################################################
# covariate coefficient cases:
b_list <- list(c(2,1,1),c(2.4,1.5,1.5),c(2.8,2,2))
# cases of b_list to run
B <- c(1,2,3)

# number of funds in population:
N_funds <- c(500,1000,2000,3000)[2]

# cases of number of observations per fund to run: 
n_obs <- c(24,36,48)[2]

# alpha cases to run
fund_alpha <- c(0.5,1,1.5)

# Dependence types of covariates to run: 
Z_type_set <- c('independent','dependent','noniid')[1]

# if there is no noise in covariate consider, set include_noise <- NULL; 
# if yes, set include_noise <- '_Noise_Covariates'
#include_noise <- NULL
include_noise <- '_Noise_Covariates'

# If running StepM and StepSPA, change FALSE below to TRUE 
include_stepM <- include_spa <- FALSE; 

source('Simulation Core.R') 

source('Plot - Performance under use of insufficient information.R')

setwd(code_folder)

################################################################################
# Performance of fwer+ under Varying sample size and number of observations.
################################################################################
# covariate coefficient cases:
b_list <- list(c(2,1,1),c(2.4,1.5,1.5),c(2.8,2,2))
# cases of b_list to run
B <- c(1,2,3)[2]

# number of funds in population:
N_funds <- c(500,1000,2000,3000)

# cases of number of observations per fund to run: 
n_obs <- c(24,36,48)

# alpha cases to run
fund_alpha <- c(0.5,1,1.5)[2]

# Dependence types of covariates to run: 
Z_type_set <- c('independent','dependent','noniid')[1]

# if there is no noise in covariate consider, set include_noise <- NULL; 
# if yes, set include_noise <- '_Noise_Covariates'
include_noise <- NULL
#include_noise <- '_Noise_Covariates'

# If running StepM and StepSPA, change FALSE below to TRUE 
include_stepM <- include_spa <- FALSE; 

source('Simulation Core.R') 

source('Plot - Performance under varying sample size and number of observations.R')

setwd(code_folder)

################################################################################
# Performance of fwer+ under Cross Sectional Dependence in Covariates
################################################################################
# covariate coefficient cases:
b_list <- list(c(2,1,1),c(2.4,1.5,1.5),c(2.8,2,2))
# cases of b_list to run
B <- c(1,2,3)[2]

# number of funds in population:
N_funds <- c(500,1000,2000,3000)[2]

# cases of number of observations per fund to run: 
n_obs <- c(24,36,48)[2]

# alpha cases to run
fund_alpha <- c(0.5,1,1.5)[2]

# Dependence types of covariates to run: 
Z_type_set <- c('independent','dependent','noniid')[3]

# if there is no noise in covariate consider, set include_noise <- NULL; 
# if yes, set include_noise <- '_Noise_Covariates'
include_noise <- NULL
#include_noise <- '_Noise_Covariates'

# If running StepM and StepSPA, change FALSE below to TRUE 
include_stepM <- include_spa <- FALSE; 

for (k in c(10,100,200)[3]) {
  source('Simulation Core.R') 
}

source('Plot - Performance under Cross Sectional Dependence in Covariates.R')

setwd(code_folder)

