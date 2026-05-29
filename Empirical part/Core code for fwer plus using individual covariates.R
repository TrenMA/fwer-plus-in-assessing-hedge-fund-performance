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

factors$Date <- as.double(factors$Date)

BS <- seq(199701,202101,100)
ES <- seq(199912,202312,100)

df_ret <- readRDS('Fake Data for Running Replication Code.RDS')

setwd(pvalue_covariates_folder)

pvalue_file_name <- paste0('Pvalues and R2 for ', n_factors, ' factor-model.RDS')

pvalue <- readRDS(pvalue_file_name) %>% bind_rows() 

covariatesfile_name <- paste0('Covariates ',n_factors,' factor-model.RDS')
  
covariates <- readRDS(covariatesfile_name)

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

covs <- c(
  lapply(characteristics_covs, function(x) c(x)),
  lapply(macro_covs, function(x) c(x)),
  lapply(return_based_covs, function(x) c(x)),
  lapply(managerial_skill, function(x) c(x))
)

pc %>% names()

samples_size <- list(); l <- 1

res_unsmoothed_recorded <- list()

for (cov_set in covs) {

  cov_name <- cov_set  
  portfolios <- list()
  sample_size <- c()
  
  for (k in seq(from = 1,to = length(ES)-1, by = Holding_horizon)) {
    
    df <- pc %>% filter(ISEnd == ES[k]) 
    
    df <- df %>% select(Name,pvalue,alpha,aum,any_of(cov_set)) %>% 
      na.omit() %>%
      filter(aum > 0)
    
    sample_size[k] <- nrow(df)
    
    pi0.var <- df %>% select(any_of(cov_set)) %>% as.data.frame()
    
    # standardize each covariate:
    
    for (v in 1:ncol(pi0.var)) {
      V <- pi0.var[,v]
      pi0.var[,v] <-  (V-mean(V))/sd(V)
    }
    
    fwer <- fwer_twosided(pvals = df$pvalue,pi0.var = pi0.var,alpha = target_set)
    
    funds_selected <- list()
    
    for (i in 1:length(target_set)) {
      df_selected <- df[fwer$rejection[[i]],] %>% filter(alpha > 0)
      funds_selected[[i]] <- df_selected$Name 
    }
    
    portfolios[[k]] <- list(funds = funds_selected,ISEnd = ES[k])
  }
  
  samples_size[[l]] <- sample_size; 
  
  ##############################################
  # portfolio of selected funds
  
  portfolios_return <- list()
  
  for (j in 1:length(target_set)) {
    rets <- list()
    for (i in 2:length(portfolios)) {
      
      portfolio_j <- portfolios[[i]]$funds[[j]]
      
      portfolio_j_ret <- df_ret %>% filter(Name %in% portfolio_j) %>% 
        filter(Date > ES[i] & Date <= ES[i+Holding_horizon]) %>%
        group_by(Date) %>% 
        summarise(ret = mean(Ret,na.rm = TRUE),
                  size = n())
      rets[[i]] <-  portfolio_j_ret   
    }
    portfolios_return[[j]] <- rets %>% bind_rows()
  }
  
  factors <- factors %>% 
    select(Date,any_of(factors_list),RF) %>%
    filter(Date >= '200001' & Date <= 202312)
  
  reg_metrics <- function(dfi){
    model_reg  <- lm(ER ~., data = dfi)
    
    # lag for HAC check, chosen based on Newey West recommendation
    
    n_obs <- nrow(dfi)
    
    L <- floor(4*(n_obs/100)^(2/9))
    
    NW <- coeftest(model_reg,vcov=NeweyWest(model_reg,
                                            prewhite = FALSE,
                                            adjust = TRUE,
                                            lag = L))
    
    data.frame(Alpha = NW[1,1]*12,
               tstat = NW[1,3],
               pvalue = NW[1,4],
               ER = mean(dfi$ER)*12)
  }
  
  df_unsmoothed <- res_unsmoothed <- list()
  
  sizes <- list()
  
  for (j in 1:length(target_set)) {
    df <- portfolios_return[[j]]
    
    df_temp <- factors %>% 
      left_join(df,by = c('Date')) 
    
    # count number of months where portfolio is empty
    n_na <- sum(is.na(df_temp$size))
    
    if (n_na == 0) {
      n_na <- NA
    }else{
      df_temp$size[is.na(df_temp$size)] <- 0
    }
    
    sizes[[j]] <- data.frame(Average = mean(df_temp$size),Min = min(df_temp$size),
                             Max = max(df_temp$size),empty = n_na)
    
    # unsmooth portfolio return
    # empty portfolio months = invest in risk-free rate:
    df_temp$ret[is.na(df_temp$ret)] <- df_temp$RF[is.na(df_temp$ret)]
    
    df_temp_unsmoothed <- df_temp %>% 
      mutate(ret_unsmoothed = ma_unsmooth(ret,H = 3)) %>%
      mutate(ER = ret_unsmoothed - RF)
    
    dfi_unsmoothed <- df_temp_unsmoothed %>% select(ER,any_of(factors_list))
    
    df_unsmoothed[[j]] <- df_temp_unsmoothed
    res_unsmoothed[[j]] <- reg_metrics(dfi = dfi_unsmoothed)
  }
  
  n_months <- nrow(df_temp)
  
  size <- sizes %>% bind_rows() %>% round(digits = 2)  %>% mutate(empty = empty*100/n_months)
  
  res_unsmoothed <- res_unsmoothed %>% 
    bind_rows() %>% 
    round(digits = 2) 
  
  res_unsmoothed <- data.frame(target_set = target_set*100,res_unsmoothed,size) %>% select(-pvalue)
  
  res_unsmoothed_recorded[[l]] <- data.frame(covariate = rep(cov_name,nrow(res_unsmoothed)),
                                             res_unsmoothed)
  
  # index for storing results for portfolio of the next covariate in the covs list
  l <- l + 1
}

################################################################################
# Summary of portfolios using individual covariate:

df_res <- res_unsmoothed_recorded %>% 
  bind_rows() %>% as_tibble() %>%
  select(Alpha,target_set,tstat,ER,Average,Max,Min,empty)

# If not empty, indicate 0
df_res$empty[is.na(df_res$empty)] <- 0

df_res <- df_res %>%
  group_by(target_set) %>% 
  summarise(alpha_mean = mean(Alpha),alpha_min = min(Alpha),alpha_max = max(Alpha),
            e1 = NA, t_mean = mean(tstat),t_min = min(tstat),t_max = max(tstat),e2 = NA,
            ER_mean = mean(ER),ER_min = min(ER),ER_max = max(ER),e3 = NA,
            Size_mean = mean(Average),Size_min = min(Min),Size_max = max(Max),empty = mean(empty))

file_name_tex <- paste0('Summary OOS of fwer-based portfolios using individual covariates for ', n_factors, ' factor-model.html')

results_table <- xtable(df_res,digits = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1))

setwd(tables_figures_folder)

print(results_table,type = 'html',file = file_name_tex,include.rownames = FALSE)

setwd(Empirical_folder)
