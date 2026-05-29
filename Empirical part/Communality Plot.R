rm(list = ls())

# Note: This file utilizes covariates and p-value of all factor models
# To produce them, you need to run two files "Covariates.R" and "PvalueCalculator.R"

library(lmtest)
library(sandwich)
library(FinTS)
library(dplyr)
library(qvalue)
library(xtable)

n_factors <- 7

Empirical_folder <- "C:/Users/User/Dropbox/Hedge fund project/Data and Code for Gihub - Empirical - clean version"

results_folder <- paste(Empirical_folder,'/', 'Results',sep = '')

if (!dir.exists(results_folder)){dir.create(results_folder)}

pvalue_covariates_folder <- paste(results_folder,'/', 'pvalue and covariates',sep = '')

if (!dir.exists(pvalue_covariates_folder)){dir.create(pvalue_covariates_folder)}

tables_figures_folder <- paste(results_folder,'/', 'tables and figures',sep = '')

if (!dir.exists(tables_figures_folder)){dir.create(tables_figures_folder)}


BS <- seq(199701,202101,100) # start date of in-samples
ES <- seq(199912,202312,100) # end date of in-samples

source('FWER Procedure.R')

source('Unsmooth return function.R')

oos_bg_date <- 200001; oos_end_date <- 202312

if (n_factors == 'Chen_9Factor') {
  
  factors_list <- c("MKT", "AGR", "BAB", "LRSK", "ROA", "TSMOM", "X10Y", "CRDT", "TRM")
  
  factors <- readRDS('Chen Nine factors.RDS')%>% 
    filter(Date >= oos_bg_date & Date <= oos_end_date) %>%  
    select(Date,any_of(factors_list),RF) %>%
    as_tibble()
  
}else{
  
  if (n_factors == 7) {
    factors_list <- c('EQ','Size','CS','DGS10C','PTFSBD','PTFSFX','PTFSCOM')
  }else if (n_factors == 4) {
    factors_list <- c('MktRF','SMB','MOM','HML')
  }else if (n_factors == 6) {
    factors_list <- c('MktRF','SMB','MOM','HML','CS','DGS10C')
  }else if (n_factors == 9) {
    factors_list <- c('MktRF','SMB','MOM','HML','CS','DGS10C','PTFSBD','PTFSFX','PTFSCOM')
  }
  
  factors <- readRDS('Fake Factors.RDS')  %>% 
    filter(Date >= oos_bg_date & Date <= oos_end_date) %>%  
    select(Date,any_of(factors_list),RF) %>%
    as_tibble()
}

df_ret <- readRDS('Fake Data for Running Replication Code.RDS')

df_ret <- df_ret %>% filter(!is.na(Ret))

# make sure that Date - one of the key variables for merging - has the same type:
df_ret$Date <- df_ret$Date %>% as.double()
factors$Date <- factors$Date %>% as.double()


# import pvalue and covariates:

setwd(pvalue_covariates_folder)

pvalue_file_name <- paste0('Pvalues and R2 for ', n_factors, ' factor-model.RDS')

pvalue <- readRDS(pvalue_file_name) %>% bind_rows() 

covariates_file_name <- paste0('Covariates ',n_factors,' factor-model.RDS')

covariates <- readRDS(covariates_file_name)

pvalue$ISEnd <- as.numeric(pvalue$ISEnd)

pc <- pvalue %>% left_join(covariates,by = c('Name','ISEnd')) 

df_rets <- df_ret %>% arrange(Date) %>% group_split(Name)

characteristics <- c('aum',
                     'IncentiveFee',
                     'LockUpPeriod',
                     'HighWaterMark',
                     'Age',
                     'ManagementFee',
                     'RedemptionNoticePeriod',
                     'MinimumInvestment'
)

pc <- pc %>% left_join(df_ret[,c('Date','Name','Ret',characteristics)],
                       by = c('Name','ISEnd'='Date'))

pc <- pc %>% 
  mutate(FundSize = log(aum)) %>% 
  as.data.frame()

# in case some data values are not defined

for (k in 2:ncol(pc)) {
  pc[is.infinite(as.vector(pc[,k])),k] <- NA
}

pc <- pc %>% as_tibble()

target_set <- c(0.001,0.005,0.01,0.05)

characteristics_covs <- c('FundSize',
                          'IncentiveFee',
                          'ManagementFee',
                          'Age',
                          'LockUpPeriod',
                          'HighWaterMark',
                          'RedemptionNoticePeriod',
                          'MinimumInvestment')

macro_covs <- c('EPU',
                'Inflation',
                'MacroRisk',
                'InvSenChgBeta',
                'DefaultSpread',
                'EMU',
                'TermSpread',
                'Liquidity'
)

return_based_covs <- c('std12',
                       'skewness12',
                       'kurtosis12',
                       'ShR12',
                       'cumret12')


managerial_skill <- c('IVol',
                      'IVolBaliWeigert',
                      'model_rsq')

cov_set <-   c(characteristics_covs,
               macro_covs,
               return_based_covs,
               managerial_skill)

################################################################################
# construct communality data

comm <- list()

for (k in seq(from = 1,to = length(ES)-1, by = 1)) {
  
  df <- pc %>% filter(ISEnd == ES[k]) 
  
  df <- df %>% select(Name,pvalue,alpha,aum,any_of(cov_set)) %>% 
    na.omit() %>%
    filter(aum > 0)
  
  pi0.var <- df %>% select(any_of(cov_set)) %>% as.data.frame()
  
  # standardize covariates (to align with the code for all covariates, FA procedure below will also standardize its inputs)
  
  for (v in 1:ncol(pi0.var)) {
    V <- pi0.var[,v]
    pi0.var[,v] <-  (V-mean(V))/sd(V)
  }
  
  library(psych)
  
  # 1. Compute the correlation matrix
  cor_matrix <- cor(pi0.var, use = "pairwise.complete.obs")
  
  # 2. Get eigenvalues
  eigen_result <- eigen(cor_matrix)
  eigenvalues <- eigen_result$values
  
  # 3. Apply Kaiser Criterion: keep components with eigenvalue > 1
  kaiser_factors <- sum(eigenvalues > 1)
  
  # Assuming your_data is a dataframe or matrix of numeric variables
  fa_result <- fa(pi0.var, nfactors = kaiser_factors, rotate = "varimax", fm = "ml",co)
  
  # View loadings
  # print(fa_result$loadings, cutoff = 0.3)
  
  comm[[k]] <- data.frame(EndIS = ES[k], t(data.frame(fa_result$communality)))%>% as_tibble()
  
} # done portfolios of funds

##############################################################################
# communality plot

library(tidyr)
library(ggplot2)

# reshape
df_long <- comm %>% bind_rows() %>%
  pivot_longer(-EndIS, names_to = "Variable", values_to = "Value") %>%
  mutate(EndIS = as.character(EndIS),
         EndIS = substring(EndIS, 1, 4))

# Create a named vector for mapping
abbr_map <- c(
  "cumret12" = "Cumulative Return",
  "std12" = "Return Volatility",
  "skewness12" = "Skewness",
  "kurtosis12" = "Kurtosis",
  "ShR12" = "Sharpe Ratio",
  "FundSize" = "Fund Size",
  "Age" = "Fund Age",
  "IncentiveFee" = "Incentive Fee",
  "ManagementFee" = "Management Fee",
  "MinimumInvestment" = "Minimum Investment",
  "LockUpPeriod" = "Lockup Period",
  "HighWaterMark" = "High-Water Mark",
  "RedemptionNoticePeriod" = "Redemption Notice Period",
  "model_rsq" = "Rsquare",
  "IVolBaliWeigert" = "Six-factor Model Idiosyncratic Risk",
  "IVol" = "Idiosyncratic Risk",
  "InvSenChgBeta" = "Investor Sentiment Change",
  "MacroRisk" = "Macro Risk",
  "Liquidity" = "Liqudity Timing",
  "Inflation" = "Inflation",
  "TermSpread" = "Term Spread",
  "DefaultSpread" = "Default Yield Spread",
  "EMU" = "Equity Market Uncertainty",
  "EPU" = "Economic Policy Uncertainty"
)

# Replace variable names with abbreviations
df_long$Variable <- abbr_map[df_long$Variable]

# compute average value for each variable (after renaming!)
avg_vals <- df_long %>%
  group_by(Variable) %>%
  summarise(Average = mean(Value, na.rm = TRUE))

# add an "Average" row into the data
df_with_avg <- df_long %>%
  bind_rows(
    avg_vals %>%
      mutate(EndIS = "Average") %>%
      rename(Value = Average)
  )

# reorder variables by average value
df_with_avg <- df_with_avg %>%
  left_join(avg_vals, by = "Variable") %>%
  mutate(Variable = reorder(Variable, Average))

df_with_avg_plot <- df_with_avg %>%
  filter(EndIS != "Average")

setwd(tables_figures_folder)

# Ensure Windows font is registered
windowsFonts(A = windowsFont("Times New Roman"))

file_name <- paste0('Communality heatmap for ',n_factors,'_factors','.png')

png(filename = file_name,
    width = 2720, height = 1300, res = 250)

ggplot(df_with_avg_plot, aes(x = EndIS, y = Variable, fill = Value)) +
  geom_tile(color = "grey90") +
  # Dashed line for Average (adjust if needed)
  scale_fill_gradientn(
    colours = c("#deebf7", "#9ecae1", "#3182bd")
  ) +
  labs(x = "Year End", y = "Covariates", fill = "Value") +
  theme_minimal(base_family = "A") +   #Apply Times New Roman
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

dev.off()

setwd(Empirical_folder)