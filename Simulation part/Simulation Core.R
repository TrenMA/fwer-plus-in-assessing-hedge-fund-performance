################################################################################
# number of iterations
iter <- 5
ptm <- 2.5 # n_positive/n_nagative alpha which is fixed in our paper

# list of seven factors
factors_list <- c('EQ','Size','CS','DGS10C','PTFSBD','PTFSFX','PTFSCOM')

# File Factors.RDS contains
# Date: YYYYMM
# EQ: Equity Market Factor
# Size: The Size Spread Factor
# DGS10C: The Bond Market Factor
# CS: The Credit Spread Factor
# PTFSBD: Return of PTFS Bond lookback straddle
# PTFSFX: Return of PTFS Currency Lookback Straddle
# PTFSCOM:Return of PTFS Commodity Lookback Straddle
# see https://people.duke.edu/~dah7/HFRFData.htm
# RF: monthly risk-free rate obtained from Kenneth R. French - Data Library
# https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html

Factors <- readRDS('Factors.RDS') %>% 
  as_tibble() %>% 
  filter(Date <= '202312' & Date >= '199701') %>%
  select(Date, all_of(factors_list),RF)

Factors$Date <- as.double(Factors$Date)

# merge with factors and calculte excess return
df_eret <- df_ret %>% bind_rows() %>%
  left_join(Factors,by='Date') %>% 
  mutate(ERet = Ret - RF) %>%
  group_split(Name)

# calculate empirical coefficients which will be used to generate simulation data
model_coefficients <- function(dfi){
  model_reg  <- lm(dfi$ERet ~ dfi$EQ + dfi$Size + dfi$CS + dfi$DGS10C + dfi$PTFSBD + dfi$PTFSFX + dfi$PTFSCOM)  
  model_sum  <- summary(model_reg)
  model_coef <- model_sum$coefficients
  
  coef <- model_coef[,1]# %>% as.vector()
  coef <- matrix(coef,nrow = 1)
  colnames(coef) <- c('alpha',"eq","s","cs","dgs","bond","fx","com")
  
  
  coef <- cbind(coef,alpha_std   = model_coef[1,2], 
                model_arsq  = model_sum$adj.r.squared,
                model_sigma = model_sum$sigma,
                sample_size = nrow(dfi))
  coef %>% as_tibble()
  
}

EmpiricalMetrics <- lapply(df_eret, function(x) model_coefficients(x)) %>% bind_rows()

FH <- Factors %>% select(-Date,-RF)
V_Factors <- cov(FH)
m_Factors <- colMeans(FH)
beta_hat <- EmpiricalMetrics %>% select("eq","s","cs","dgs","bond","fx","com")
V_beta <- cov(beta_hat)
m_beta <- colMeans(beta_hat)
sigma_hat <- EmpiricalMetrics$model_sigma
m_sigma_hat <- median(sigma_hat)
Obs_number <- EmpiricalMetrics$sample_size
mean_obs_number <- round(mean(Obs_number))
median_obs_number <- round(median(Obs_number))
delta <- mean(sapply(beta_hat, function(x) mean(x))[2:7])
sigma_G <- mean(sapply(FH, function(x) sd(x))[2:7])

for (Z_type in Z_type_set) {
  
  if (Z_type == 'independent') {
    rho_set <- 0; k <- NULL
  }
  
  if (Z_type == 'dependent') {
    rho_set <- 0.5; k <- NULL
  }
  
  if (Z_type == 'noniid') {
    rho_set <- c(0,0.25,0.5,0.75,0.9)
  }
  
  for (b_case in B) {
    for (n_funds in N_funds) {
      for (TimeLength in n_obs) {
        for (rho in rho_set) {
          
          st <- Sys.time()
          b <- b_list[[b_case]] # beta for pi0Z
          
          n_Z <- length(b) - 1
          
          for (alpha in fund_alpha) {
            
            alpha_positives <- alpha
            alpha_negatives <- -alpha_positives
            
            target_set <- c(0.001, 0.005,seq(0.01,0.2,0.01))
            
            obs_per_fund <- TimeLength # in a balanced data setting, all time series has the same length as the factors
            use_empirical_alpha <- FALSE
            
            output_file_name <- paste0('Simualtion_with_b_',b_case,'_Zrho_',rho,
                                      '_obs_per_fund_',TimeLength,'_n_funds_',n_funds,
                                      '_alpha_',alpha,'_StepSPA_',include_spa,
                                      '_stepM_',include_stepM,include_noise,'.RDS')
            
            if (!is.null(k)) {
              output_file_name <- paste0('Simualtion_with_b_',b_case,'_Zrho_',rho,
                                        '_obs_per_fund_',TimeLength,'_n_funds_',n_funds,
                                        '_alpha_',alpha,'_StepSPA_',include_spa,
                                        '_stepM_',include_stepM,include_noise,
                                        '_noniid_k_',k,'.RDS')
            }
            
            Z_generator <- function(type = 'independent', k = NULL){
              
              Z <- matrix(NA,nrow = n_funds,ncol = n_Z)
              
              if (type == 'independent') {
                
                for (i in 1:n_Z) {
                  Z[,i] <- rnorm(n = n_funds,mean = 0,sd = 1)
                }
                
              }else if (type == 'dependent'){
                
                Z <- MASS::mvrnorm(n = n_funds,mu = c(0,0),
                                   Sigma = matrix(c(1,rho,rho,1),nrow = 2,byrow = TRUE))
              }else if (type == 'noniid'){
                
                
                Sigma_block <- matrix(rho,nrow =  k,ncol =  k)
                diag(Sigma_block) <- 1
                
                Z1 <- unlist(
                  lapply(1:(n_funds/k), function(x) MASS::mvrnorm(n = 1,mu = rep(0,k),Sigma = Sigma_block))
                )
                
                Z2 <- unlist(
                  lapply(1:(n_funds/k), function(x) MASS::mvrnorm(n = 1,mu = rep(0,k),Sigma = Sigma_block))
                )
                
                # Combine into data frame
                Z <- data.frame(Z1 = Z1, Z2 = Z2) %>% as.matrix()
              }
              
              Z
            }
            
            simulation_results <- list()
            
            set.seed(1)
            
            for (i in 1:iter) {
              
              
              pi0Z <- function(Z,b){
                
                if (length(b) <= 2) {
                  
                  Z <- as.vector(Z)
                  
                  pi0z <- 1/(1+exp(-b[1] - Z * b[-1] ))
                  
                }else{
                  b1 <- matrix(b[-1],ncol = 1)
                  
                  pi0z <- 1/(1+exp(-b[1] - (Z %*% b1 ) ) )
                }
                
                pi0z
              }
              
              Z <- Z_generator(type = Z_type,k = k)
              
              pi0z <- pi0Z(Z = Z,b = b)
              
              h <- sapply(pi0z, function(p) rbinom(n = 1,size = 1,prob = 1 - p)) # zero and non-zero alphas
              
              alts <- which(h == 1) # non-zero alpha index
              
              n_pos <- round(ptm*sum(h==1)/(ptm+1) )
              
              pos <- sample(x = alts,size = n_pos,replace = FALSE) # index where funds having positive alpha
              
              # true alphas for funds
              a <- h # mainly to keep track the zero alpha, zero alphas will be updated below
              
              a[pos] <- alpha_positives 
              
              nes <- setdiff(alts,pos) # index where funds having negative alpha
              
              a[nes] <- alpha_negatives
              
              pvalue_alpha_ds <- Simulated_Data_Generator(number_of_funds = n_funds,
                                                          alpha_real = a,
                                                          TimeLength = TimeLength,
                                                          m_sigma_hat = m_sigma_hat,
                                                          PA_only = FALSE)
              
              pvalue_alpha <- pvalue_alpha_ds$PnA
              
              PA <- data.frame(pvalue_alpha,alpha_true = a)
              
              df <- list(PA = PA,Z = as.data.frame(Z))
              
              ds_ret <- pvalue_alpha_ds$ds_ret
              
              
              each_iteration_fwer_metrics <- function(df,ds,target_set = seq(0.005,0.1,0.005)){
                Z <- df$Z
                PA <- df$PA
                
                # fwer+ procedure starts
                fwer <- fwer_twosided(pvals = PA$pvalues,pi0.var = Z,alpha = target_set)
                
                fwer_selected <- lapply(fwer$rejection, function(x) PA[x,] %>% filter(alphas_est > 0))
                # fwer+ procedure ended
                
                fwe <- sapply(fwer_selected, function(x) any(x$alpha_true <= 0))
                fwe_detected_rate <- sapply(fwer_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                
                ##################################################################
                # One covariate case
                fwer1 <- fwer_twosided(pvals = PA$pvalues,pi0.var = Z[,1],alpha = target_set)
                fwer1_selected <- lapply(fwer1$rejection, function(x) PA[x,] %>% filter(alphas_est > 0))
                fwe1 <- sapply(fwer1_selected, function(x) any(x$alpha_true <= 0))
                fwe1_detected_rate <- sapply(fwer1_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                
                seed_current <- .Random.seed
                Z_uninfo <- rnorm(n = n_funds,mean = 0,sd = 1)
                
                fwer_uninfo <- fwer_twosided(pvals = PA$pvalues,pi0.var = Z_uninfo,alpha = target_set)      
                fwer_uninfo_selected <- lapply(fwer_uninfo$rejection, function(x) PA[x,] %>% filter(alphas_est > 0))
                fwe_uninfo <- sapply(fwer_uninfo_selected, function(x) any(x$alpha_true <= 0))
                fwe_uninfo_detected_rate <- sapply(fwer_uninfo_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                
                if (!is.null(include_noise)) {
                  res_noises <- list()
                  for (si in 1:3) {
                    s <- c(1.5,3.0,4.5)[si]
                    
                    eta05 <- rnorm(n = n_funds,mean = 0,sd = s)
                    xi05 <- rnorm(n = n_funds,mean = 0,sd = s)
                    
                    Z_noise <- data.frame(V1 = Z[,1]+eta05, V2 = Z[,1]+xi05)
                    
                    fwer_noise <- fwer_twosided(pvals = PA$pvalues,pi0.var = Z_noise,alpha = target_set)      
                    fwer_noise_selected <- lapply(fwer_noise$rejection, function(x) PA[x,] %>% filter(alphas_est > 0))
                    fwe_noise <- sapply(fwer_noise_selected, function(x) any(x$alpha_true <= 0))
                    fwe_noise_detected_rate <- sapply(fwer_noise_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                    
                    res_2noise <- data.frame(fwe_2noise = fwe_noise,fwe_2noise_detected_rate = fwe_noise_detected_rate)
                    
                    names(res_2noise) <- c(paste('fwe_2noise_',s,sep = '',collapse = ''),
                                           paste('fwe_2noise_detected_rate_',s,sep = '',collapse = '')
                    )
                    
                    fwer_noise <- fwer_twosided(pvals = PA$pvalues,pi0.var = Z_noise[,1],alpha = target_set)      
                    fwer_noise_selected <- lapply(fwer_noise$rejection, function(x) PA[x,] %>% filter(alphas_est > 0))
                    fwe_noise <- sapply(fwer_noise_selected, function(x) any(x$alpha_true <= 0))
                    fwe_noise_detected_rate <- sapply(fwer_noise_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                    
                    res_1noise <- data.frame(fwe_1noise = fwe_noise,fwe_1noise_detected_rate = fwe_noise_detected_rate)
                    
                    names(res_1noise) <- c(paste('fwe_1noise_',s,sep = '',collapse = ''),
                                           paste('fwe_1noise_detected_rate_',s,sep = '',collapse = '')
                    )
                    
                    pca <- prcomp(Z_noise,scale. = TRUE)
                    Z_pca <- pca$x[,1]
                    
                    fwer_noise <- fwer_twosided(pvals = PA$pvalues,pi0.var = Z_pca,alpha = target_set)      
                    fwer_noise_selected <- lapply(fwer_noise$rejection, function(x) PA[x,] %>% filter(alphas_est > 0))
                    fwe_noise <- sapply(fwer_noise_selected, function(x) any(x$alpha_true <= 0))
                    fwe_noise_detected_rate <- sapply(fwer_noise_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                    
                    res_pca <- data.frame(fwe_pca = fwe_noise,fwe_pca_detected_rate = fwe_noise_detected_rate)
                    
                    names(res_pca) <- c(paste('fwe_pca_',s,sep = '',collapse = ''),
                                        paste('fwe_pca_detected_rate_',s,sep = '',collapse = '')
                    )
                    
                    res_noises[[si]] = data.frame(res_1noise,res_2noise,res_pca)
                    
                  }
                  
                  if (length(res_noises) <= 1) {
                    res_noise <- res_noises[[1]]
                  }else{
                    res_noise <- res_noises[[1]]
                    for (si in 2:length(res_noises)) {
                      res_noise <- data.frame(res_noise,res_noises[[si]])
                    }
                  }
                  
                  for (s in c(0.5,1.0,1.5)) {
                    eta05 <- rnorm(n = n_funds,mean = 0,sd = s)
                    xi05 <- rnorm(n = n_funds,mean = 0,sd = s)
                    
                    
                    Z_noise <- data.frame(V1 = Z[,1]+eta05, V2 = Z[,2]+xi05)
                    
                    fwer_noise <- fwer_twosided(pvals = PA$pvalues,pi0.var = Z_noise,alpha = target_set)      
                    fwer_noise_selected <- lapply(fwer_noise$rejection, function(x) PA[x,] %>% filter(alphas_est > 0))
                    fwe_noise <- sapply(fwer_noise_selected, function(x) any(x$alpha_true <= 0))
                    fwe_noise_detected_rate <- sapply(fwer_noise_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                    
                    res_noise_new <- data.frame(fwe_noise,fwe_noise_detected_rate)
                    
                    names(res_noise_new) <- c(paste('fwer_noise_',s,sep = '',collapse = ''),
                                              paste('fwer_noise_detected_rate_',s,sep = '',collapse = '')
                    )
                    res_noise <- data.frame(res_noise,res_noise_new)
                  }
                }            
                
                if (include_spa) {
                  
                  spa <- step_spa(ds = ds_ret,targets = target_set,obs_per_fund = obs_per_fund,q = q_spa)
                  
                  spa_selected <- lapply(spa$rejection, function(x) PA[x,])
                  spa_fwe <- sapply(spa_selected, function(x) any(x$alpha_true <= 0))
                  spa_detected_rate <- sapply(spa_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                  
                }
                if (include_stepM) {
                  
                  t <- Sys.time()
                  SM <- stepM(ds = ds_ret,targets = target_set,obs_per_fund = obs_per_fund,q = q_stepM)
                  difftime(Sys.time(),t)
                  
                  SM_selected <- lapply(SM$rejection, function(x) PA[x,])
                  stepM_fwe <- sapply(SM_selected, function(x) any(x$alpha_true <= 0))
                  stepM_detected_rate <- sapply(SM_selected, function(x)  sum(x$alpha_true > 0)/ sum(PA$alpha_true > 0))
                }
                .Random.seed <<- seed_current 
                
                # some population stats that might be useful: 
                # true proportion of zero alpha fund (pi0)
                pi0 <- sum(PA$alpha_true == 0)/length(PA$alpha_true)
                # proportion of positive alpha funds (pip)
                pip <- sum(PA$alpha_true > 0)/length(PA$alpha_true)
                # proportion of negative alpha funds (pim)
                pim <- sum(PA$alpha_true < 0)/length(PA$alpha_true)
                
                
                #return baseline dataframe
                res <- data.frame(target = target_set, 
                                  fwe,fwe_detected_rate,
                                  fwe1,fwe1_detected_rate,
                                  fwe_uninfo,fwe_uninfo_detected_rate,
                                  pip=pip,pi0=pi0,pim=pim)
                
                if (include_spa) {
                  # add result for step_spa
                  res <- data.frame(res,spa_fwe,spa_detected_rate)
                }
                
                if (include_stepM) {
                  # add result for stepM
                  res <- data.frame(res,stepM_fwe,stepM_detected_rate)
                } 
                
                if (!is.null(include_noise)) {
                  # add result for the case covariates with noise
                  res <- data.frame(res,res_noise)
                } 
                
                res %>% as_tibble()
              }
              
              simulation_results[[i]] <- each_iteration_fwer_metrics(df = df,ds = ds_ret,target_set = target_set)
              
            }
            
            print(difftime(Sys.time(),st))
            
            # save results in the result folder
            
            setwd(folder_simulation_results)
            saveRDS(simulation_results,file = output_file_name)
            
            # returning to code folder
            setwd(code_folder)
          }
        }
      }
    }
  }
}