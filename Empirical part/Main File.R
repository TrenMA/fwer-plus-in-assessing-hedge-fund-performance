rm(list = ls())

library(lmtest)
library(sandwich)
library(FinTS)
library(dplyr)
library(qvalue)
library(xtable)

Empirical_folder <- "input the path to your 'Empirical part' folder here"

# create a folder in Empirical_folder named "Results" for result ouput

results_folder <- paste(Empirical_folder,'/', 'Results',sep = '')

if (!dir.exists(results_folder)){dir.create(results_folder)}

# create a folder inside "Results" named "pvalue and covariates" to store 
# p-value and covariates for all factor models.

pvalue_covariates_folder <- paste(results_folder,'/', 'pvalue and covariates',sep = '')

if (!dir.exists(pvalue_covariates_folder)){dir.create(pvalue_covariates_folder)}

# create a folder inside "Results" named "tables and figures" to store 
# all table and figure results.

tables_figures_folder <- paste(results_folder,'/', 'tables and figures',sep = '')

if (!dir.exists(tables_figures_folder)){dir.create(tables_figures_folder)}

################################################################################
# 1. Produce results for 7-factor model
#    portfolios selected fwer plus using all covariates 
#    and re-balance every year (Holding_horizon = 1).
#    Results include:
#      - Summary OOS performance of FWERtau portolios using individual covariates  
#      - OOS performance of (equally weighted and past alpha weighted) FWERtau portfolios using all covariates
#      - Performance of (equal weight, equal weight plus, and past alpha weighted) portfolios of eligible funds
#      - Alpha and wealth evolution comparisons
#      - Sub-sample performance reports for FWERtau portfolios
#      - Outperformance and retention rates table for FWERtau portfolios
# 
n_factors <- 7

# first we calculate p-value and covariates for 7-factor model
# results for this step are stored in the folder pvalue_covariates_folder

setwd(Empirical_folder)

source('PvalueCalculator.R')

source('CovariatesCalculator.R')

# now we produce empirical results
# all empirical results will be returned in the folder tables_figures_folder

Holding_horizon <- 1 #years

setwd(Empirical_folder)

source('Core code for fwer plus using individual covariates.R')

source('Core code for fwer plus using all covariates.R')

################################################################################
# 2. Produce results for re-balancing every 2, 3 and 4 years.
#    Results include: 
#      - OOS performance of (equally weighted and past alpha weighted) 
#        FWERtau portfolios rebalancing every 2, 3 and 4 years

n_factors <- 7

for (Holding_horizon in c(2,3,4)) {
  
  setwd(Empirical_folder)
  
  source('Core code for fwer plus using all covariates.R')
}

################################################################################
# 3. Produce results for 4, 6, 9 and Chen et al (2025)'s 9-factor models
#    portfolios selected fwer plus using all covariates 
#    and re-balance every year (Holding_horizon = 1).
#    Results include: 
#      - OOS performance of (equally weighted and past alpha weighted) FWERtau portfolios

Holding_horizon <- 1 #years

for (n_factors in c(4,6,9,'Chen_9Factor')) {
  
  setwd(Empirical_folder)
  # calculate p-value and covariate
  source('PvalueCalculator.R')

  source('CovariatesCalculator.R')  

  setwd(Empirical_folder)
  source('Core code for fwer plus using all covariates.R')
}
################################################################################
# 4. Plot communality heatmap
#    Open and run file 'Communality Plot.R'
