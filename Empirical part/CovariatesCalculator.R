BS <- seq(199701,202101,100) # begin YYYYMM of in-sample periods
ES <- seq(199912,202312,100) # end YYYYMM of in-sample periods

setwd(Empirical_folder)

# import factors

factors <- readRDS('Fake Factors.RDS') %>% filter(Date >= BS[1])

factors$Date <- as.double(factors$Date)

six_factor_Bali <- readRDS('Fake Bali six factor.RDS') %>% filter(Date >= BS[1])

Chen_nine_factor <- readRDS('Fake Chen Nine factors.RDS') %>% filter(Date >= BS[1])

if (n_factors == 'Chen_9Factor') {
  factors_list <- c("MKT", "AGR", "BAB", "LRSK", "ROA", "TSMOM", "X10Y", "CRDT", "TRM")
}else if (n_factors == '7') {
  factors_list <- c('EQ','Size','CS','DGS10C','PTFSBD','PTFSFX','PTFSCOM')
}else if (n_factors == '4') {
  factors_list <- c('MktRF','SMB','MOM','HML')
}else if (n_factors == '6') {
  factors_list <- c('MktRF','SMB','MOM','HML','CS','DGS10C')
}else if (n_factors == '9') {
  factors_list <- c('MktRF','SMB','MOM','HML','CS','DGS10C','PTFSBD','PTFSFX','PTFSCOM')
}

# import data

df0 <- readRDS('Fake Data for Running Replication Code.RDS')

df <- df0 %>% bind_rows() %>% 
  select(Date,Name,Ret) %>%
  left_join(factors[,c('Date','RF')],by='Date') %>%
  mutate(ERet = Ret - RF)

IVolBaliWeigert <- function(dfj,Window = 36){
  
  df_temp <- tail(dfj,n = Window) %>% # consider only the most recent Window months
    select(Date,eret) %>% 
    left_join(six_factor_Bali,by = 'Date') %>% 
    select(-Date, -RF)
  md <- lm(eret ~ .,data = df_temp)
  s <- summary(md)
  
  IVol <- s$sigma^2
  
  IVol
}

IVol <- function(dfj,n_factors,Window = 36){

  if (n_factors != 'Chen_9Factor') {
    ft <- factors %>% select(Date,any_of(factors_list))
  }else{
    ft <- Chen_nine_factor
  }
  
  ft$Date <- as.numeric(ft$Date)
  df_temp <- tail(dfj,n = Window) %>% # consider only the most recent Window months
    select(Date,eret) %>% 
    left_join(ft,by = 'Date') %>% 
    select(-Date)
  
  md <- lm(eret ~ .,data = df_temp)
  
  s <- summary(md)
  
  IVol <- s$sigma^2
  
  IVol
}

variances <- function(x,window = 12){
  
  x <- tail(x,window) # consider only the most recent Window months
  
  x <- na.omit(x)
  
  n <- length(x)
  if (n >= window) {
    res <- var(x) 
  }else{
    res <- NA
  }
  
  res
}

kurtoses <- function(x,window = 12){
  
  x <- tail(x,window) # consider only the most recent Window months
  
  x <- na.omit(x)
  
  n <- length(x)
  if (n >= window) {
    if (sd(x) != 0) {
      res <- e1071::kurtosis(x)
    }else{
      res <- NA
    }
  }else{
    res <- NA
  }
  res
}


skewness <- function(x,window = 12){
  
  x <- tail(x,window) # consider only the most recent Window months
  
  x <- na.omit(x)
  
  n <- length(x)
  
  if (n >= window) {
    if (sd(x) != 0) {
      res <- e1071::skewness(x)
    }else{
      res <- NA
    }
  }else{
    res <- NA
  }
  
  res
}


SR <- function(x,window = 12){
  
  x <- tail(x,window) # consider only the most recent Window months
  
  x <- na.omit(x)
  
  n <- length(x)
  
  if (n >= window) {
    if (sd(x) != 0) {
      res <- mean(x)/sd(x)
    }else{
      res <- NA
    }
  }else{
    res <- NA
  }
  
  res
}

cumret <- function(x,window = 12){
  
  x <- tail(x,window)
  
  x <- na.omit(x)
  
  n <- length(x)
  
  if (n >= window) {
    res <- sum(x)
  }else{
    res <- NA
  }
  
  res
}

return_based_covariates <- list()

for (i in 1:length(ES)) {
  
  dfs <- df %>% filter(Date <= ES[i]) %>% group_split(Name)
  
  eligible_list <- which(sapply(dfs, function(x) nrow(x)) >= 36) 
  
  df_temp <- lapply(eligible_list, function(k) dfs[[k]]) %>% bind_rows()
  
  dfs0 <- df_temp %>% group_split(Name)
  
  covariate <- list()
  
  for (j in 1:length(dfs0)) {
    
    dfj <- dfs0[[j]] %>% arrange(Date) %>% rename('eret' = 'ERet')
    
    n <- length(na.omit(dfj$Ret))
    
    if (n >= 12) {
      
      var12 <- variances(x = dfj$eret,window = 12)
      
      cumret12 <- cumret(x = dfj$eret,window = 12)
      
      skewness12 <- skewness(x = dfj$eret,window = 12)
      
      ShR12 <- SR(x = dfj$eret,window = 12)
      
      kurtosis12 <- kurtoses(x = dfj$eret,window = 12)
      
      IVolBaliBaliWeigert <- IVolBaliWeigert(dfj = dfj,Window = 36)
      
      IVol_model <- IVol(dfj = dfj,n_factors = n_factors,Window = 36)
      
      covariate[[j]] <- data.frame(ISEnd = ES[i], Name = dfj$Name[1],
                                   std12 = sqrt(var12),
                                   cumret12,
                                   skewness12,
                                   kurtosis12,
                                   ShR12 = ShR12, 
                                   IVolBaliWeigert = IVolBaliBaliWeigert,
                                   IVol = IVol_model
      )
    }
  }
  return_based_covariates[[i]] <- covariate %>% bind_rows()
}

df_return_based_covariates <- return_based_covariates %>% bind_rows() %>% as_tibble()

################################################################################

macro_covs <- c('EPU',
                'Inflation',
                'MacroRisk',
                'InvSenChgBeta',
                'DefaultSpread',
                'EMU',
                'TermSpread',
                'Liquidity'
)

macro_variables <- readRDS('Fake macro variables.RDS') %>% select(Date,any_of(macro_covs))

Sentiment_variables <- readRDS('Fake Varibales for Calculating Investor Sentiment Change.RDS')

ISC <- function(end_sample,df_sentiment){
  
  df_pca <- df_sentiment %>% 
    select(yearmo,pdnd,ripo,nipo,cefd,s) %>% 
    arrange(yearmo) %>% 
    mutate(pdnd = pdnd - lag(pdnd),
           ripo = ripo - lag(ripo),
           nipo = nipo - lag(nipo),
           cefd = cefd - lag(cefd),
           s = s - lag(s)) %>% 
    na.omit() %>%
    filter(yearmo <= end_sample)
  
  pca_result <- prcomp(df_pca[,-1], scale. = TRUE)
  
  InvestorSentiment <- pca_result$x[,1]
  
  df_pca$InvestorSentiment <- InvestorSentiment
  
  df_macro <- df_sentiment %>% 
    select(yearmo,consdur,consserv,recess,employ,cpi) 
  
  df1 <- df_pca %>% select(yearmo,InvestorSentiment) %>% 
    left_join(df_macro,by = 'yearmo')
  
  md <- lm(InvestorSentiment ~.,data = df1[,-1])
  
  InvestorSentimentChanges <- md$residuals
  
  df1$InvestorSentimentChanges <- InvestorSentimentChanges
  
  df1 %>% 
    mutate(Date = yearmo) %>% 
    select(Date,InvestorSentimentChanges)
}

res_sentiment_change <- list()

if (n_factors == 'Chen_9factor') {
  factors_df <- Chen_nine_factor
}else{
  factors_df <- factors
}

for (i in 1:length(BS)) {
  
  df_temp <- df %>% 
    filter(Date <= ES[i] & Date >= BS[i]) %>% 
    select(Name,Date,Ret) %>%
    group_split(Name)
  
  df_isc <- ISC(end_sample = ES[i],df_sentiment =  Sentiment_variables)
  
  isc <- list()
  
  for (j in 1:length(df_temp)) {
    
    dfj <- df_temp[[j]] %>% 
      filter(Date >= BS[i] & Date <= ES[i]) %>% 
      left_join(df_isc,by = 'Date') %>%
      left_join(factors_df %>% select(Date,any_of(factors_list),RF),
                by = 'Date') %>%
      mutate(ERet = Ret - RF) %>%
      select(-Ret,-RF)
    
    if (nrow(dfj) >= 36) {
      
      md <- lm(ERet ~. , data = dfj %>% select(-Name,-Date,))
      
      ivsc <- md$coefficients['InvestorSentimentChanges'] %>% as.numeric()
      
      isc[[j]] <- data.frame(Name = dfj$Name[1],
                             ISEnd = ES[i],
                             InvSenChgBeta = ivsc)
    }
  }
  
  res_sentiment_change[[i]] <- isc %>% bind_rows() %>% as_tibble()
}

df_InSenChg <- res_sentiment_change %>% bind_rows()

loading_covariates <- list()

for (i in 1:length(BS)) {
  print(i)
  
  df_temp <- df %>% 
    filter(Date >= BS[i] & Date <= ES[i]) %>% 
    mutate(ERet = Ret - RF) %>%
    select(Name,Date,ERet) %>%
    left_join(macro_variables,by = 'Date') %>%
    group_split(Name)
  
  res_temps <- list()
  
  for (f in 1:length(df_temp)) {
    
    dff <-  df_temp[[f]] 
    
    #An example of dff
    # A tibble: 36 × 10
    #Name    Date  ERet   EPU Inflation MacroRisk DefaultSpread   EMU TermSpread  Liquidity
    #<chr>  <dbl> <dbl> <dbl>     <dbl>     <dbl>         <dbl> <dbl>      <dbl>      <dbl>
    
    if (nrow(dff) >= 36) {
      res_temp <- data.frame(ISEnd = ES[i],Name = dff$Name[1])
      
      for (j in 4:ncol(dff)) {
        
        if (nrow(na.omit(dff[,c(3,j)])) >= 36) {
          
          if (names(dff)[j] %in% c('Liquidity')) {
            
            dff_temp <- dff[,c(3,j)] %>% as.data.frame() 
            x <- dff_temp[,2] %>% as.vector()
            x2 <- x^2
            
            md <- lm(dff_temp$ERet ~  x+ x2)
            b <- md$coefficients[3] %>% as_tibble()
            
          }else{
            md <- lm(ERet ~ . , data = dff[,c(3,j)])
            b <- md$coefficients[2] %>% as_tibble()
          }
          
        }else{
          b <- data.frame(X = NA)
        }
        names(b) <- names(dff)[j]
        res_temp <- data.frame(res_temp,b)
      }
      
      res_temps[[f]] <- res_temp
    }
  }
  
  loading_covariates[[i]] <- res_temps %>% bind_rows() %>% as_tibble()
}

df_loading_covariates <- loading_covariates %>% bind_rows()


df_covariates <- df_return_based_covariates %>% 
  left_join(df_InSenChg,by = c('Name','ISEnd')) %>% 
  left_join(df_loading_covariates,c('Name','ISEnd'))

df_covariates$ISEnd <- as.numeric(df_covariates$ISEnd)

file_name <- paste0('Covariates ',n_factors,' factor-model.RDS')

setwd(pvalue_covariates_folder)

saveRDS(df_covariates,file = file_name)
