BS <- seq(199701,202101,100) # start dates of in-sample
ES <- seq(199912,202312,100) # end dates of in-sample

IS_length <- 36 # in-sample window length (in months)

setwd(Empirical_folder)

source('FWER Procedure.R')

source('Unsmooth return function.R')


if (n_factors == 'Chen_9Factor') {
  
  factors_list <- c("MKT", "AGR", "BAB", "LRSK", "ROA", "TSMOM", "X10Y", "CRDT", "TRM")
  
  factors <- readRDS('Fake Chen Nine factors.RDS')%>% as_tibble() %>% 
    na.omit() %>% select(Date,any_of(factors_list),RF)
  
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
  
  factors <- readRDS('Fake Factors.RDS') %>% select(Date,any_of(factors_list),RF)
}


df <- readRDS('Fake Data for Running Replication Code.RDS')  

df <- df %>% filter(!is.na(Ret)) %>% filter(Date >= BS[1])

df$Date <- df$Date %>% as.double()

factors$Date <- factors$Date %>% as.double()

pvalue_calculator <- function(bg,ed){
  
  # bg: In-sample start date
  # ed: In-sample end date
  
  dfs0 <- df %>% 
    filter(Date >= bg & Date <= ed) %>% 
    group_by(Name) %>% 
    group_split()
  
  # consider only funds with full observations over the in-sample period
  
  full_obs_list <- which(sapply(dfs0, function(x) nrow(x)) >= IS_length) 
  
  df1 <- lapply(full_obs_list, function(k) dfs0[[k]]) %>% 
    bind_rows() 
  
  dfs2 <- df1 %>%
    arrange(Name,Date) %>%
    group_split(Name)
  
  dfs_temp <- list()
  
  for (j in 1:length(dfs2)) {
    
    df_temp <- dfs2[[j]] %>% arrange(Date)
    
    ret_unsmoothed <- ma_unsmooth(ret = df_temp$Ret,H = 3)      
    
    df_temp$Ret <- ret_unsmoothed
    
    dfs_temp[[j]] <- df_temp
  }
  
  dfs2 <- dfs_temp
  
  
  # collect only Date, Name, ER and factors
  
  dfs <- dfs2 %>% bind_rows() %>%
    left_join(factors,by='Date') %>%
    mutate(ER = Ret - RF) %>%
    select(Date,Name,ER,any_of(factors_list)) %>%
    group_split(Name)
  
  res_i <- function(dfi){
    
    dfi_ret <- dfi %>% select(-Date,-Name) # retain only ER and factors for regression
    
    model_reg  <- lm(ER ~., data = dfi_ret)
    
    model_sum  <- summary(model_reg)
    
    model_coef <- model_sum$coefficients
    
    # lag for HAC, chosen based on Newey-West 1987 paper
    
    L <- floor(4*(IS_length/100)^(2/9))
    
    NW <- coeftest(model_reg,vcov=NeweyWest(model_reg,
                                            prewhite = FALSE,
                                            adjust = TRUE,
                                            lag = L))
    
    res <- data.frame(Name = dfi$Name[1],
                      ISStart = dfi$Date[1],
                      ISEnd = tail(dfi$Date,1),
                      sample_size = dim(dfi)[1],
                      pvalue = model_coef[1,4],
                      pvalueNW = NW[1,4],
                      tstat = NW[1,3],
                      alpha  = model_coef[1,1],
                      model_rsq  = model_sum$r.squared)
    
    res
  }
  
  ouput <- lapply(dfs, function(dfi) res_i(dfi))
  
  ouput <- ouput %>% bind_rows() %>% as_tibble()
  
  ouput
}


Results <- list()

for (i in 1:length(ES)) {
  print(i)
  Results[[i]] <- pvalue_calculator(bg = BS[i],ed = ES[i])
}

setwd(pvalue_covariates_folder)

pvalue_file_name <- paste0('Pvalues and R2 for ', n_factors, ' factor-model.RDS',sep = '',collapse = '')

saveRDS(Results %>% bind_rows(),file = pvalue_file_name)

setwd(Empirical_folder)