
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
# construct portfolios

portfolios_ew <- portfolios_ewp <- portfolios <- list()
sample_size <- c()
#eligible_funds <- list()

for (k in seq(from = 1,to = length(ES)-1, by = Holding_horizon)) {
  
  # collect only fund with in-sample end at ES[k]:
  df <- pc %>% filter(ISEnd == ES[k]) 
  
  df <- df %>% select(Name,pvalue,alpha,aum,any_of(cov_set)) %>% 
    na.omit() %>%
    filter(aum > 0)
  
  #eligible_funds[[k]] <- data.frame(Name = df$Name, EndIS = ES[k])
  
  sample_size[k] <- nrow(df)
  
  # covariates for estimating pi0(z):
  
  pi0.var <- df %>% select(any_of(cov_set)) %>% as.data.frame()
  
  # standardize the covariates
  
  for (v in 1:ncol(pi0.var)) {
    V <- pi0.var[,v]
    pi0.var[,v] <-  (V-mean(V))/sd(V)
  }
  
  # fwer+ starts
  
  fwer <- fwer_twosided(pvals = df$pvalue,pi0.var = pi0.var,alpha = target_set)
  
  funds_selected <- list()
  
  for (i in 1:length(target_set)) {
    df_selected <- df[fwer$rejection[[i]],] %>% filter(alpha > 0)
    funds_selected[[i]] <- df_selected$Name 
  } 
  # fwer+ ended
  
  portfolios[[k]] <- list(funds = funds_selected,ISEnd = ES[k],
                          n_factor = n_factors)
  
  # eligible funds for the benchmark portfolio:
  portfolios_ew[[k]] <- list(funds = df$Name,ISEnd = ES[k])
  
  # eligible and positive IS alpha funds for the benchmark portfolio:
  portfolios_ewp[[k]] <- list(funds = df$Name[df$alpha > 0],ISEnd = ES[k])
} # done portfolios of funds

##############################################################################
# portfolio of equally weighted selected funds

portfolios_return <- list()

for (j in 1:length(target_set)) {
  rets <- list()
  for (i in 1:length(portfolios)) {
    
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


# weighted by past alpha:

portfolios_return_past_alpha_weighted <- list()

for (j in 1:length(target_set)) {
  rets <- list()
  for (i in 1:length(portfolios)) {
    
    portfolio_j <- portfolios[[i]]$funds[[j]]
    
    past_alpha <- pc %>% filter(Name %in% portfolio_j) %>% 
      filter(ISEnd == portfolios[[i]]$ISEnd) %>%
      select(Name,alpha) %>% rename('past_alpha' = 'alpha')
    
    portfolio_j_ret <- df_ret %>% filter(Name %in% portfolio_j) %>% 
      filter(Date > ES[i] & Date <= ES[i+Holding_horizon]) %>%
      left_join(past_alpha,by = 'Name') %>%
      group_by(Date) %>% 
      summarise(ret = weighted.mean(Ret,past_alpha,na.rm = TRUE),
                size = n())
    
    rets[[i]] <-  portfolio_j_ret   
  }
  portfolios_return_past_alpha_weighted[[j]] <- rets %>% bind_rows()
}

##############################################################################
# equally weighted portfolio of eligible funds

rets_past_alpha_wp <- rets_ewp <- rets_ew <- list()

for (i in 1:length(portfolios_ew)) {
  
  # ew
  
  portfolio_ew <- portfolios_ew[[i]]$funds
  
  portfolio_ew_ret <- df_ret %>% filter(Name %in% portfolio_ew) %>% 
    filter(Date > ES[i] & Date <= ES[i+Holding_horizon]) %>%
    group_by(Date) %>% 
    summarise(ret = mean(Ret,na.rm = TRUE),
              size = n())
  
  rets_ew[[i]] <-  portfolio_ew_ret   
  
  #ewp: equally weighted positive IS alpha funds
  
  portfolio_ewp <- portfolios_ewp[[i]]$funds
  
  portfolio_ewp_ret <- df_ret %>% filter(Name %in% portfolio_ewp) %>% 
    filter(Date > ES[i] & Date <= ES[i+Holding_horizon]) %>%
    group_by(Date) %>% 
    summarise(ret = mean(Ret,na.rm = TRUE),
              size = n())
  rets_ewp[[i]] <-  portfolio_ewp_ret   
  
  # past alpha weighted positive IS alpha funds 
  
  portfolio_ewp <- portfolios_ewp[[i]]$funds
  
  past_alpha <- pc %>% filter(Name %in% portfolio_ewp) %>% 
    filter(ISEnd == portfolios_ewp[[i]]$ISEnd) %>%
    select(Name,alpha) %>% rename('past_alpha' = 'alpha')
  
  portfolio_ret_past_alpha_wp <- df_ret %>% filter(Name %in% portfolio_ewp) %>% 
    filter(Date > ES[i] & Date <= ES[i+Holding_horizon]) %>%
    left_join(past_alpha,by = 'Name') %>%
    group_by(Date) %>% 
    summarise(ret =weighted.mean(Ret,past_alpha,na.rm = TRUE),
              size = n())
  
  rets_past_alpha_wp[[i]] <-  portfolio_ret_past_alpha_wp   
  
}

sample_ew_ret <- rets_ew %>% bind_rows()
sample_ewp_ret <- rets_ewp %>% bind_rows()
sample_past_alpha_wp_ret <- rets_past_alpha_wp %>% bind_rows()

################################################################################

reg_metrics <- function(dfi){
  
  dfi <- dfi %>% select(ER,any_of(factors_list))
  
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

################################################################################
# calculate performance metrics for portfolios:

# create a function to calculate performance evaluation

portfolios_evaluation <- function(portfolios_ret){
  
  # portfolios_ret: list of portfolios each with: Date, ret, size
  
  df_unsmoothed <- res <- res_unsmoothed <- list()
  
  sizes <- list()
  
  for (j in 1:length(target_set)) {
    
    df <- portfolios_ret[[j]]
    
    df_temp <- factors %>% 
      left_join(df,by = c('Date')) 
    
    n_na <- sum(is.na(df_temp$size))
    
    if (n_na == 0) {
      n_na <- NA
    }else{
      df_temp$size[is.na(df_temp$size)] <- 0
    }
    
    sizes[[j]] <- data.frame(Average = mean(df_temp$size),Min = min(df_temp$size),
                             Max = max(df_temp$size),empty = n_na)
    
    df_temp <- df_temp %>% mutate(ER = ret - RF)
    
    df_temp$ret[is.na(df_temp$ret)] <- df_temp$RF[is.na(df_temp$ret)]
    
    df_temp_unsmoothed <- df_temp %>% 
      mutate(ret_unsmoothed = ma_unsmooth(ret,H = 3)) %>%
      mutate(ER = ret_unsmoothed - RF)
    
    dfi_unsmoothed <- df_temp_unsmoothed %>% select(ER,any_of(factors_list))
    dfi_unsmoothed$ER[is.na(dfi_unsmoothed$ER)] <- 0
    df_unsmoothed[[j]] <- df_temp_unsmoothed
    res_unsmoothed[[j]] <- reg_metrics(dfi = dfi_unsmoothed)
  }
  
  
  res_unsmoothed <- res_unsmoothed %>% bind_rows() %>% round(digits = 2) 
  
  # for alpha and wealth evolution plots
  
  port_wealth <- lapply(df_unsmoothed, function(x) x  %>%
                          mutate(cumret = (cumprod(1+ret_unsmoothed/100)-1)*100) %>%
                          mutate(wealth = cumprod(1+ret_unsmoothed/100)) %>%
                          mutate(Date = as.Date(paste0(Date,'01'),format = '%Y%m%d')) 
  )
  
  n_months <- nrow(df_temp)
  
  size <- sizes %>% bind_rows() %>% round(digits = 2)  %>% mutate(empty = empty*100/n_months)
  
  res_perf <- data.frame(target_set = target_set*100,res_unsmoothed,size) %>% 
    select(target_set ,Alpha, tstat, ER, Average, Min, Max,empty) # also report empty rate   
  
  
  list(port_wealth=port_wealth,
       res_perf=res_perf, 
       ports_unsmoothed_ret = df_unsmoothed)
}

# Baseline FWER portfolio:

perf_report_baseline_ports <- portfolios_evaluation(portfolios_ret = portfolios_return)

file_name_tex <- paste0('Performance of fwer portfolios using all covariates rebalancing every ',
                        Holding_horizon,' year - ',n_factors,
                        ' factor-model.html')

latex_table <- xtable(perf_report_baseline_ports$res_perf,digits = c(0,1, 1, 1, 1, 0, 0, 0,1))

setwd(tables_figures_folder)

print(latex_table,type = 'html',file = file_name_tex,include.rownames = FALSE)

# past alpha weighted FWER portfolios 

perf_report_past_alpha_weighted_ports <- portfolios_evaluation(portfolios_ret = portfolios_return_past_alpha_weighted)

file_name_tex <- paste0('Performance past alpha weighted fwer portfolios using all covariates rebalancing every ',
                        Holding_horizon,'-year - ',n_factors,
                        ' factor-model.html')

latex_table <- xtable(perf_report_past_alpha_weighted_ports$res_perf,digits = c(0,1, 1, 1, 1, 0, 0, 0,1))

setwd(tables_figures_folder)

print(latex_table,type = 'html',file = file_name_tex,include.rownames = FALSE)

################################################################################
#
# the below results are reported for only portfolios rebalanced yearly and use 7-factor model 

if (n_factors == 7 & Holding_horizon == 1) {
  
  ##############################################################################
  # sub-period report:
  
  P1 <- c(200001,200711)
  
  P2 <- c(200801,201512)
  
  P3 <- c(201601,202312)
  
  Ps <- list(P1,P2,P3) # list of sub-periods
  
  sub_report <- function(d1,d2,df_ret_unsmoothed){
    # df_ret_unsmoothed:  unsmoothed-return and factors dataframe
    # d1 begin period; d2: end period  
    S <- df_ret_unsmoothed %>% 
      filter(Date >= d1 & Date <= d2)
    
    S$ER[is.na(S$ER)] <- 0
    
    s <- S$size
    
    if (sum(s == 0) > 0) {
      e = sum(s == 0)/length(S)*100
    }else{
      e <- NA        
    }
    
    S <- S  %>% 
      select(ER,any_of(factors_list))
    
    res <- data.frame(Period = paste(d1,'-',d2,sep = '',collapse = ''),
                      target = target_set[j]*100,
                      reg_metrics(dfi = S),
                      AverageSize = mean(s),MinSize = min(s), MaxSize = max(s), empty = e)
    res
  }
  
  # report for baseline FWER portfolios:
  
  S_report <- list()
  
  for (j in 1:length(target_set)) {
    
    sub_res <- list()
    
    p <- 0
    
    for (P in Ps) {
      p <- p + 1
      unsmoothed_ret_factors <- perf_report_baseline_ports$ports_unsmoothed_ret[[j]]
      sub_res[[p]] <- sub_report(d1 = P[1],d2 = P[2],
                                 df_ret_unsmoothed = unsmoothed_ret_factors)
    }
    S_report[[j]] <- sub_res %>% bind_rows()
  }
  
  sub_sample_res <- S_report %>% bind_rows() %>% group_by(Period) %>% group_split()
  
  setwd(tables_figures_folder)
  
  for (rs in sub_sample_res) {
    
    file_name_tex <- paste0('Performance over subperiod ',rs$Period[1],' for fwer portfolios ',
                            ' rebalancing every ',Holding_horizon,'-year.html')
    
    rs <- rs %>% select(target, Alpha, tstat, ER,AverageSize, MinSize, MaxSize, empty)
    latex_table <- xtable(rs,digits = c(0,1, 1, 1, 1, 0, 0,0,1))
    print(latex_table,type = 'html',file = file_name_tex,include.rownames = FALSE)
  }
  
  # report for past alpha weighted FWER portfolios
  
  S_report <- list()
  
  for (j in 1:length(target_set)) {
    
    sub_res <- list()
    
    p <- 0
    
    for (P in Ps) {
      p <- p + 1
      unsmoothed_ret_factors <- perf_report_past_alpha_weighted_ports$ports_unsmoothed_ret[[j]]
      sub_res[[p]] <- sub_report(d1 = P[1],d2 = P[2],
                                 df_ret_unsmoothed = unsmoothed_ret_factors)
    }
    S_report[[j]] <- sub_res %>% bind_rows()
  }
  
  sub_sample_res <- S_report %>% bind_rows() %>% group_by(Period) %>% group_split()
  
  setwd(tables_figures_folder)
  for (rs in sub_sample_res) {
    
    file_name_tex <- paste0('Performance of alpha weighted fwer portfolios rebalancing every ',
                            Holding_horizon,'-year over subperiod ',rs$Period[1],'.html')
    
    rs <- rs %>% select(target, Alpha, tstat, ER,AverageSize, MinSize, MaxSize, empty)
    latex_table <- xtable(rs,digits = c(0,1, 1, 1, 1, 0, 0,0,1))
    print(latex_table,type = 'html',file = file_name_tex,include.rownames = FALSE)
  }
  
  #################################################################################
  # ew portfolios:
  
  sample_portfolios_evaluation <- function(portfolio_ret){
    
    # portfolio_return: a dataframe consisting of: Date, ret, size
    
    df <- portfolio_ret
    
    n_na <- sum(is.na(df$size))
    
    size <- data.frame(Average = mean(df$size),
                       Min = min(df$size),
                       Max = max(df$size),
                       empty_rate = n_na/nrow(df) # empty rate
    )
    
    df_temp <- factors %>% left_join(df,by = c('Date')) %>% mutate(ER = ret - RF) 
    
    df_temp_unsmoothed <- df_temp %>% 
      mutate(ret_unsmoothed = ma_unsmooth(ret,H = 3)) %>%
      mutate(ER = ret_unsmoothed - RF)
    
    port_wealth <- df_temp_unsmoothed %>%
      mutate(cumret = (cumprod(1+ret_unsmoothed/100)-1)*100) %>%
      mutate(wealth = cumprod(1+ret_unsmoothed/100)) %>%
      mutate(Date = as.Date(paste0(Date,'01'),format = '%Y%m%d'))
    
    df_unsmoothed <- df_temp_unsmoothed %>% select(ER,any_of(factors_list))
    df_unsmoothed$ER[is.na(df_unsmoothed$ER)] <- 0
    
    res_unsmoothed <- reg_metrics(dfi = df_unsmoothed)
    
    res_unsmoothed <- data.frame(res_unsmoothed,size)
    
    list(port_wealth = port_wealth, 
         perf_report = res_unsmoothed, 
         ports_unsmoothed_ret = df_unsmoothed)
    
  }
  
  # equally weighted sample of eligible funds:
  
  sample_ew_port <- sample_portfolios_evaluation(portfolio_ret = sample_ew_ret)
  
  # equally weighted sample of eligible funds having positive in-sample alpha:
  
  sample_ewp_port <- sample_portfolios_evaluation(portfolio_ret = sample_ewp_ret)
  
  # past alpha  weighted sample of eligible funds having positive in-sample alpha:
  
  sample_AlphaWP_port <- sample_portfolios_evaluation(portfolio_ret = sample_past_alpha_wp_ret)
  
  # combine results of all sample portfolios:
  
  res_ews_umsmoothed <- bind_rows(sample_ew_port$perf_report,
                                  sample_ewp_port$perf_report,
                                  sample_AlphaWP_port$perf_report) %>% select(-pvalue)
  
  row.names(res_ews_umsmoothed)  <-  c('EW','EWP','AlphaWP')
  
  file_name_tex <- paste0('Performance of portfolios EW EWP and aWP',
                          ' rebalancing every ',Holding_horizon,' year.html')
  latex_table <- xtable(res_ews_umsmoothed,digits = c(0, 1, 1, 1, 1, 0, 0, 0))
  
  print(latex_table,type = 'html',file = file_name_tex,include.rownames = TRUE)
  
  
  ####################################################################
  # wealth evolution plot
  
  windowsFonts(A = windowsFont("Times New Roman"))
  
  png(filename = paste0('Wealth Evolution of portfolios under ',n_factors,' factor-model.png'),
      width = 2720 ,height = 1300,res = 200)
  par(family = "A",mar=c(4,4,1,1))  # Apply Times New Roman to all elements
  
  ports <- perf_report_baseline_ports$port_wealth %>% bind_rows()
  ports_weighted <- perf_report_past_alpha_weighted_ports$port_wealth %>% bind_rows()
  
  ymax <- ceiling(max(ports_weighted$wealth,ports$wealth))
  
  wealth_baseline <- perf_report_baseline_ports$port_wealth
  wealth_past_alpha_weighted <- perf_report_past_alpha_weighted_ports$port_wealth
  
  sew_port <- sample_ew_port$port_wealth
  sewp_port <- sample_ewp_port$port_wealth
  sawp_port <- sample_AlphaWP_port$port_wealth
  
  grid_col <- adjustcolor("#c6dbef", alpha.f = 0.5)
  
  plot(sew_port$Date, sew_port$wealth, type = "l",
       xlab = "Year", ylab = "Wealth ($)",ylim = c(0,ymax),
       #main = "Cumulative Return over Time",
       col = "black", lwd = 1,lty=2)
  
  hs <- c(0:ymax)
  
  for (i in 1:length(hs)) {
    abline(h = hs[i],col = grid_col)  
  }
  
  
  lines(y = sew_port$wealth,x = sew_port$Date,col = 'black',lty = 2)
  lines(y = sewp_port$wealth,x = sewp_port$Date,col = 'black',lty = 1)
  lines(y = sawp_port$wealth,x = sawp_port$Date,col = 'black',lty = 3)
  cls <- c('red','green','blue','purple')# rainbow(4) # blues <- c("#1f77b4", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef")
  
  for (i in 1:4) {
    lines(x = wealth_baseline[[i]]$Date, y = wealth_baseline[[i]]$wealth,col = cls[[i]])  
  }
  
  for (i in 1:4) {
    lines(x = wealth_past_alpha_weighted[[i]]$Date, y = wealth_past_alpha_weighted[[i]]$wealth,col = cls[[i]],lty = 3)  
  }
  
  legend_labels <- c(
    "FWER0.1%",
    "FWER0.5%",
    "FWER1%",
    "FWER5%",
    "EW Plus",
    "EW",
    expression(Past~alpha~weighted~FWER0.1*"%"),
    expression(Past~alpha~weighted~FWER0.5*"%"),
    expression(Past~alpha~weighted~FWER1*"%"),
    expression(Past~alpha~weighted~FWER5*"%"),
    expression(Past~alpha~weighted~Plus)
  )
  
  legend("topleft",
         legend = legend_labels,
         
         col = c(cls[1:4],          # FWER
                 "black", "black", # EW+ and EW
                 cls[1:4],
                 "black"),
         lty = c(rep(1, 4),         # FWER
                 1, 2,              # EW+, EW
                 rep(3, 4),      # weighted FWER
                 3),                # Past alpha W+
         x.intersp = 0.5,
         ncol = 2)
  
  dev.off()
  
  ####################################################################
  # alpha evolution plot
  
  alpha_evolution <- function(df){
    
    df <- df %>% 
      mutate(ER = ret_unsmoothed-RF) %>% 
      filter(Date >= as.Date('2000-01-01')) %>% 
      arrange(Date)
    
    date_set <- df$Date
    
    a <- c()
    
    for (k in 60:length(date_set)) {
      
      dfk <- df %>% 
        filter(Date <= as.Date(date_set[k]) ) %>% 
        select(ER,any_of(factors_list))
      mk <- reg_metrics(dfi = dfk)
      a <- c(a,mk$Alpha)
    }
    data.frame(Date = date_set[-1:-59], alpha = a) %>% as_tibble()
  }
  
  dfa <- lapply(wealth_baseline, function(x) alpha_evolution(x))
  
  dfaw <- lapply(wealth_past_alpha_weighted, function(x) alpha_evolution(x))
  
  ewa <- alpha_evolution(df = sew_port)
  
  ewpa <- alpha_evolution(df = sewp_port)
  
  awpa <- alpha_evolution(df = sawp_port)
  
  windowsFonts(A = windowsFont("Times New Roman"))
  
  png(filename = paste0('Alpha Evolution of portfolios using covariates ', n_factors,' factor-model.png'),
      width = 2720 ,height = 1300,res = 200)
  par(family = "A",mar=c(4,4,1,1))  # Apply Times New Roman to all elements
  
  dfas <- dfa %>% bind_rows()
  dfasw <- dfaw %>% bind_rows()
  
  ymax <- ceiling(max(dfasw$alpha,dfas$alpha)*1.01)
  ymin <- floor(min(dfasw$alpha,dfas$alpha,0)*1.01)
  plot(x = ewa$Date, y = ewa$alpha, type = "l",
       xlab = "Year", ylab = "Alpha (%)",
       ylim = c(ymin,ymax),
       col = "black", lwd = 1,lty = 2)
  
  hs <- c(ymin:ymax)
  
  for (i in 1:length(hs)) {
    abline(h = hs[i],col = grid_col)  
  }
  
  lines(y = ewa$alpha,x = ewa$Date,col = 'black',lty = 2)
  lines(y = ewpa$alpha,x = ewpa$Date,col = 'black',lty = 1)
  lines(y = awpa$alpha,x = awpa$Date,col = 'black',lty = 3)
  
  for (i in 1:4) {
    lines(y = dfa[[i]]$alpha,x = dfa[[i]]$Date,col = cls[[i]],lty = 1)  
  }
  
  for (i in 1:4) {
    lines(y = dfaw[[i]]$alpha,x = dfaw[[i]]$Date,col = cls[[i]],lty = 3)  
  }
  
  legend("topright",
         legend = legend_labels,
         
         col = c(cls[1:4],          # FWER
                 "black", "black", # EW+ and EW
                 cls[1:4],
                 "black"),
         lty = c(rep(1, 4),         # FWER
                 1, 2,              # EW+, EW
                 rep(3, 4),      # weighted FWER
                 3),                # Past alpha W+
         x.intersp = 0.5,
         ncol = 2)
  
  dev.off()
  
  ################################################################################
  # Fund style analysis
  # we do mot unsmooth the fund return for this experiment
  # since some funds are delisted in the middle of the year
  
  alpha_fund_selected <- c()
  category_selected <- c()
  
  for (j in 1:length(target_set)) {
    
    category <- rets <- list()
    
    for (i in seq(1,length(portfolios),Holding_horizon)) {
      
      portfolio_j <- portfolios[[i]]$funds[[j]]
      
      if (length(portfolio_j) > 0) { 
        
        # in our experiment with real data, we always have length(portfolio_j) > 0  
        
        portfolio_j_ret <- df_ret %>% 
          filter(Name %in% portfolio_j) %>% 
          filter(Date > ES[i] & Date <= ES[min(i+Holding_horizon,length(ES))]) %>% 
          select(Date,Name,Ret,PrimaryCategory) %>%
          group_split(Name)
        
        perf_holding_period <- list() # to store performance of each fund in portfolio_j
        
        category[[i]] <- lapply(portfolio_j_ret, function(x) x[1,] %>% select(Date,PrimaryCategory)) %>% 
          bind_rows()
        
        n_delist <- 0 # count number of delisted funds over the holding period
        
        # for each fund in the portfolio_j:
        
        for (m in 1:length(portfolio_j_ret)) {
          
          f <- portfolio_j_ret[[m]]
          
          f_temp <- f %>% 
            left_join(factors,by = 'Date') %>% 
            arrange(Date)  %>% 
            mutate(ERet = Ret - RF) %>% 
            select(ERet,any_of(factors_list))
          
          md <- lm(ERet ~.,data = f_temp)
          
          full_period <- sum((factors$Date > ES[i]) & 
                               (factors$Date <= ES[min(i+Holding_horizon,length(ES))]))
          
          if (nrow(f) < full_period) {
            n_delist <- n_delist + 1
            alpha_ex_delist <- NA # so the fund is excluded from metrics that calculate among non-delisted funds
          }else{
            alpha_ex_delist <- md$coefficients[1]
          }
          
          perf_holding_period[[m]] <- data.frame(Year = substring(f$Date[1],first = 1,last = 4), 
                                                 Name = f$Name[1],
                                                 alpha = md$coefficients[1],
                                                 AveRet = mean(f$Ret),AveERet = mean(f_temp$ERet),
                                                 alpha_ex_delist = alpha_ex_delist
          )
        }
        
      }else{
        # if the fwer+ does not select any funds, return a null for the list perf_holding_period
        perf_holding_period <- NULL 
      }
      
      # if perf_holding_period is not null, calculate metrics:
      
      if (!is.null(perf_holding_period)) {
        
        rets[[i]] <- perf_holding_period %>% 
          bind_rows() %>% 
          mutate(S = (alpha > 0)) %>% 
          group_by(Year) %>% 
          summarise(positive_alpha_rate = mean(alpha > 0),
                    positive_alpha = sum(alpha > 0),
                    positive_alpha_rate_ex_delist = mean(alpha_ex_delist > 0,na.rm = TRUE),
                    positive_alpha_ex_delist = sum(alpha_ex_delist > 0,na.rm = TRUE),
                    n_delist = n_delist,
                    positive_eret = sum(AveERet>0),
                    n_fund = n()
          )   
      }
    }
    
    alpha_fund_selected[[j]] <- rets %>% bind_rows() %>% mutate(target = target_set[j])
    
    category_selected[[j]] <- category
  }
  
  df_whole <- list()
  
  for (i in 1:(length(alpha_fund_selected))) {
    
    df_temp <- alpha_fund_selected[[i]]
    
    df_whole[[i]] <- df_temp %>% 
      
      mutate(retention_prop = 1 - n_delist/n_fund,
             positive_eret_prop = positive_eret/n_fund) %>%
      summarise(target = target[1]*100,
                n_fund = mean(n_fund),
                retention_rate = mean(retention_prop)*100,
                # out-performance rate among funds with full return 
                pos_alpha_ex_delist_rate = mean(positive_alpha_rate_ex_delist,na.rm = TRUE)*100,
                # out-performance rate among all funds selected one 
                pos_eret_rate_all = mean(positive_eret_prop)*100
      )
    
  }
  
  setwd(tables_figures_folder)
  
  df_whole <- df_whole %>% bind_rows()
  
  names(df_whole) <- c('FWER taget',
                       'Average number of selected funds',
                       'OOS retention rate',
                       'OOS outperforming rate in retained funds',
                       'OOS outperforming rate in selected funds')
  
  library(xtable)
  
  latex_table <- xtable(df_whole,digits = c(0,1, 0, 0, 0, 0))
  
  file_name_tex <- paste0('Outperformance and retion rates of funds selected by fwer plus using all covariates over holding period of ', 
                          Holding_horizon,' years.html')
  
  print(latex_table,type = 'html',file = file_name_tex,include.rownames = FALSE)
  
  ################################################################################
  # selected funds by style plot
  
  # we report for only one FWER target, results for other targets look similar:
  
  tau <- which(target_set == 0.05) # change here to report for other fwer target in target_set
  
  category_tau <- category_selected[[tau]] %>% bind_rows()
  
  
  style_decomposition_plot <- function(category,h_gap,b_gap){
    
    category$PrimaryCategory[category$PrimaryCategory == 'Event-Driven'] <- 'Event Driven'
    category$PrimaryCategory[category$PrimaryCategory == 'Macro'] <- 'Global Macro'
    
    category$PrimaryCategory <- factor(category$PrimaryCategory)
    
    
    plot_data <- category %>%
      mutate(
        Year = as.numeric(substring(Date,first = 1,last = 4)), # get the year
        subperiod_id = floor((Year - min(Year)) / 8), # sub-period ID
        year_start = min(Year) + 8 * subperiod_id, 
        year_end = pmin(year_start + 7, max(Year)),
        Subperiod = paste0(year_start, "-", year_end)
      ) %>% select(Subperiod,PrimaryCategory, Year)
    
    plot_data <- plot_data %>% group_by(Subperiod, PrimaryCategory) %>%
      summarise(Freq = n(), # how many funds of a category appear in a period
                YearsInBlock = n_distinct(Year), # how many year the category appears over the sub-period 
                .groups = "drop") %>%
      mutate(AvgFreq = Freq / YearsInBlock) # average number of funds per year
    
    subperiods <- unique(plot_data$Subperiod)
    categories <- levels(category$PrimaryCategory)
    
    x_vals <- 1:(length(subperiods))
    
    line_types <- (1:length(categories))+1
    point_types <- rep(c(16, 18), length.out = length(categories))
    cols <- rainbow(length(categories))
    
    
    y_max <- round(1.3*max(plot_data$AvgFreq))
    
    plot(x_vals, rep(NA, length(x_vals)),
         type = "n", xaxt = "n",xaxs = "r",
         xlab = "Subperiod", ylab = "Average number of funds per year",
         ylim = c(0, y_max))
    
    axis(2)
    
    text(x_vals,y = par("usr")[3] - b_gap,
         labels = subperiods,
         srt = 45,
         adj = 1,
         xpd = TRUE)
    
    hs <- seq(from=0,to = y_max, by = h_gap)
    
    for (h in hs) {
      abline(h = h,col = grid_col)  
    }
    
    for(v in x_vals) {
      abline(v = v, col = grid_col)
    }
    
    for(i in 1:length(categories)) {
      
      temp <- plot_data %>% filter(PrimaryCategory == categories[i])
      
      y_vals <- rep(NA, length(subperiods))
      
      y_vals[match(temp$Subperiod, subperiods)] <- temp$AvgFreq
      
      lines(x_vals, y_vals, lty = line_types[i], col = cols[i], lwd = 1)
      
      points(x_vals, y_vals, pch = point_types[i], col = cols[i], cex = 1)
    }
    
    legend("topleft",
           legend = categories,
           col = cols,
           lty = line_types,
           pch = point_types,
           ncol = 2,
           bty = "o",
           bg = "white")
  }
  
  setwd(tables_figures_folder)
  
  windowsFonts(A = windowsFont("Times New Roman"))
  file_name <- paste0('Funds Selected by fwer plus by style for fwer target ',target_set[tau],' - ', n_factors,
                      ' factor-model.png')
  
  png(filename = file_name, width = 2720 ,height = 1300,res = 200)
  
  par(family = "A", mar = c(5,5,1,1), mgp = c(3.75,1,0))
  style_decomposition_plot(category = category_tau,h_gap = 5,b_gap = 0.4)
  dev.off()
  
  ############################################################################
  # Fund style decomposition: sample of all eligible funds
  
  #category <- pvalue %>% select(Name,ISEnd) %>% filter(ISEnd > 199912 & ISEnd <= 202312) %>%
  #  left_join(df_ret[,c('Date','Name','PrimaryCategory')],by = c('ISEnd'='Date','Name')) %>% 
  #  mutate(Date = ISEnd) %>%
  #  select(-Name ,- ISEnd)
  
  category_ew <- list() 
  for (i in 1:length(portfolios_ew)) {
    portfolio_ew_ret <- df_ret %>% 
      filter(Name %in% portfolios_ew[[i]]$funds) %>% 
      filter(Date > ES[i] & Date <= ES[min(i+Holding_horizon,length(ES))]) %>% 
      select(Date,Name,Ret,PrimaryCategory) %>%
      group_split(Name)
    category_ew[[i]] <- lapply(portfolio_ew_ret, function(x) x[1,] %>% select(Date,PrimaryCategory)) %>% 
      bind_rows()
  }
  
  category_ew <- category_ew %>% bind_rows()
  
  file_name <- paste0('Eligible funds by Style.png')
  png(filename = file_name, width = 2720 ,height = 1300,res = 200)
  par(family = "A", mar = c(5,5,1,1), mgp = c(3.75,1,0))
  style_decomposition_plot(category = category_ew,h_gap = 100,b_gap = 10)
  dev.off()
  
  ################################################################################
  # portfolio returns decomposition by style:
  
  portfolios_j <- list() # portfolio corresponding to FWER = tau
  
  for (i in 1:length(portfolios)) {
    portfolios_j[[i]] <- portfolios[[i]]$funds[[tau]]
  }
  
  return_decomposition_by_style <- function(ports){
    # ports: a list of funds in the portfolio for each year
    
    port_ret <- list()
    
    for (i in 1:length(ports)) {
      
      port_ret[[i]] <- df_ret %>% filter(Name %in% ports[[i]]) %>% 
        filter(Date > ES[i] & Date <= ES[i+Holding_horizon])
    }
    
    port_ret <- port_ret %>% bind_rows() %>% select(Date,Name,Ret,PrimaryCategory)
    
    port_ret$PrimaryCategory[port_ret$PrimaryCategory == 'Event-Driven'] <- 'Event Driven'
    port_ret$PrimaryCategory[port_ret$PrimaryCategory == 'Macro'] <- 'Global Macro'
    
    portfolios_return_style <- port_ret %>% 
      select(Date,PrimaryCategory,Ret) %>% 
      group_by(Date,PrimaryCategory) %>%
      summarise(Ret = mean(Ret,na.rm = TRUE),size = n()) %>% 
      ungroup()
    
    df_style <- portfolios_return_style %>% group_by(PrimaryCategory) %>% group_split()
    
    df_report <- list()
    
    # for each style
    
    for (j in 1:length(df_style)) {
      
      df_temp_unsmoothed <- df_style[[j]] %>% arrange(Date) %>%
        mutate(Ret = ma_unsmooth(ret = Ret,H = 3)) %>%
        left_join(factors, by = 'Date') %>% mutate(ER = Ret-RF)
      
      report <- function(df_temp_unsmoothed,d1,d2){
        
        S <- df_temp_unsmoothed %>% 
          filter(Date >= d1 & Date < d2)
        
        if (nrow(S) < 8) {
          res <- NULL
        }else{
          # S$ER[is.na(S$ER)] <- 0
          s <- S$size
          
          present_time = nrow(S)#/sum(ES >= d1 & ES < d2)/12
          
          S <- S  %>% 
            select(ER,any_of(factors_list))
          
          res <-  data.frame(Period = paste(d1,'-',d2,sep = '',collapse = ''), 
                             style = df_temp_unsmoothed$PrimaryCategory[1],
                             reg_metrics(dfi = S), 
                             Presence = present_time,
                             AverageSize = mean(s),MinSize = min(s), MaxSize = max(s)
          )
        }
        
        res
      }
      
      df_report[[j]] <- report(df_temp_unsmoothed = df_temp_unsmoothed,d1 = 200001,d2 = 202401)
    }
    
    df_report %>% bind_rows() %>% select(style, Alpha, tstat, ER, Presence,AverageSize)
    
  }
  
  setwd(tables_figures_folder)
  
  # report for portfolio corresponding to FWER = tau
  
  report_table <- return_decomposition_by_style(ports = portfolios_j) 
  
  latex_table <- xtable(report_table,digits = c(0,0,1, 1, 1, 0, 1))
  
  file_name_tex <- paste0('Performance by style of funds selected by fwer plus - rebalancing every ',
                          Holding_horizon,' year with fwer target ',target_set[tau],'.html')
  
  print(latex_table,type = 'html',file = file_name_tex,include.rownames = FALSE)
  
  
  # report for the portfolio of eligible funds:  
  
  sample_portfolios_ew <- list()
  
  for (i in 1: length(portfolios_ew)) {
    sample_portfolios_ew[[i]] <- portfolios_ew[[i]]$funds
  }
  
  report_table <- return_decomposition_by_style(ports = sample_portfolios_ew)
  
  latex_table <- xtable(report_table,digits = c(0,0,1, 1, 1, 0, 1))
  
  file_name_tex <- paste0('Performance by style of porfolio of all eligible funds - rebalancing every ',
                          Holding_horizon,'-year.html')
  
  print(latex_table,type = 'html',file = file_name_tex,include.rownames = FALSE)
  
} # condition Holding_horizon == 1 ends

setwd(Empirical_folder)

