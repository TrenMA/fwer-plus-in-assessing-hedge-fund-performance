GenerateBootsIndex <- function(n,q){
  ID <- c()
  i <- 1
  ID[i] <- sample(x = 1:n,size = 1,replace = FALSE)
  i <- i+1
  while (i <= n) {
    U <- runif(n = 1,min = 0,max = 1)
    if (U < q) {
      ID[i] <- sample(x = 1:n,size = 1,replace = FALSE)
    }else{
      ID[i] <- ID[i-1] + 1
      if (ID[i] > n) {
        ID[i] <- 1
      }
    }
    i <- i + 1
  }
  return(ID)
}

pi0Z <- function(Z,b,form = 'logistic',pi0 = NULL){
  if (length(b) <= 2) {
    Z <- as.vector(Z)
    if (form == 'logistic') {
      pi0z <- 1/(1+exp(-b[1] - Z * b[-1] ))
    }else if (form == 'sine') {
      if (!is.null(pi0)) {
        pi0z <- (pi0 + 0.4*sin(Z))
      }else{
        pi0z <- (0.6 + 0.4*sin(Z))
      }
    }
  }else{
    b1 <- matrix(b[-1],ncol = 1)
    if (form == 'logistic') {
      pi0z <- 1/(1+exp(-b[1] - (Z %*% b1 ) ) )
    }else if (form == 'sine') {
      if (!is.null(pi0)) {
        pi0z <- (pi0 + 0.4*(sin(Z[,1]+Z[,2])) )
      }else{
        pi0z <- (0.6 + 0.4*(sin(Z[,1]+Z[,2])) )
      }
    }
  }
  pi0z
}

Simulated_Data_Generator <- function(number_of_funds,alpha_real,
                                     m_sigma_hat,
                                     TimeLength,
                                     PA_only = TRUE){
  # Return: simulated estimate alpha and its p-value
  #         for all funds
  # Inputs: obs_per_fund could be given or will be assigned as unbalanced type (will consider in later steps)
  
  PnA <- list()
  
  fb <- Factors_beta_simulate(TimeLength = TimeLength,number_of_funds = number_of_funds)
  Factors_simulate <- fb$Factors_simulate
  beta_simulate <- fb$beta_simulate
  
  ds <- list()
  ds_ret <- list()
  for (i in 1:number_of_funds) {
    epsilon <- rnorm(obs_per_fund,mean = 0,sd = m_sigma_hat)
    ERet <- alpha_real[i] +  as.matrix(Factors_simulate) %*% t(as.matrix(beta_simulate[i,])) + epsilon
    dfi <- data.frame(ERet,Factors_simulate)
    ds_ret[[i]] <- dfi
    # since in this simulation, epsilon is iid, we use simple regression p-value 
    reg_model <- lm(ERet ~ ., data =  dfi)  
    sum_model <- summary(reg_model)
    
    ds[[i]] <- reg_model$residuals + reg_model$coefficients[1]
    
    PnA[[i]] <- data.frame(pvalues = sum_model$coefficients[1,4],alphas_est = sum_model$coefficients[1,1])
  }
  if (PA_only) {
    res <- bind_rows(PnA) %>% as_tibble()
  }else{
    res <- list(PnA = bind_rows(PnA) %>% as_tibble(),
                ds = ds, ds_ret = ds_ret)
  }
  
  return(res)  
}

Factors_beta_simulate <- function(TimeLength,number_of_funds){
  Factors_simulate <- replicate(n = TimeLength,expr =  m_Factors +t(t(chol(V_Factors)) %*% rnorm(length(m_Factors))),simplify = TRUE)
  row.names(Factors_simulate) <- names(m_Factors)
  Factors_simulate <- as.data.frame(t(Factors_simulate)) %>% as_tibble()
  
  beta_simulate <- replicate(n = number_of_funds,expr = m_beta+t(t(chol(V_beta)) %*% rnorm(length(m_beta))),simplify = TRUE)
  row.names(beta_simulate) <- names(m_beta)
  beta_simulate <- as.data.frame(t(beta_simulate)) %>% as_tibble()
  return(list(Factors_simulate=Factors_simulate,beta_simulate=beta_simulate))
}


omega_mu_hat <- function(d,q,studentized = TRUE){
  # returns the std of a vector d and inverse of block length
  #  inputs: d
  n <- length(d)
  g <- sapply(0:(n-1), function(i) sum((d -  mean(d))[1:(n-i)]*  (d -  mean(d))[(i+1):n]) / n) # length of n
  i <- 1:(n-1) 
  K <- ((n-i)/n) * (1-q)^i + (i/n) * (1-q)^(n-i) # length of n-1
  om <- sqrt(g[1] + 2*sum(K*g[2:n]))
  
  m <- mean(d)
  mu_hat = m*(sqrt(n)*m <= - om*sqrt(2*log(log(n))))
  
  if (studentized) {
    statistic <- sqrt(n)*m / om
  }else{
    statistic <- sqrt(n)*m
  }
  
  data.frame(omega = om, mu_hat = mu_hat,d_bar = m,statistic = statistic)
}


omega_hat <- function(d,q){
  # returns the std of a vector d and inverse of block length
  #  inputs: d
  n <- length(d)
  g <- sapply(0:(n-1), function(i) sum((d -  mean(d))[1:(n-i)]*  (d -  mean(d))[(i+1):n]) / n) # length of n
  i <- 1:(n-1) 
  K <- ((n-i)/n) * (1-q)^i + (i/n) * (1-q)^(n-i) # length of n-1
  om <- sqrt(g[1] + 2*sum(K*g[2:n]))
  
  om
}

step_spa <- function(ds,q,targets = seq(0.005,0.1,0.005),obs_per_fund,studentized = TRUE,parallel = TRUE){
  
  dt_calculator <- function(ds,Bindex){
    
    X1 <- ds[[1]][Bindex,] %>% select(-ERet) %>% as.matrix()
    X <- cbind(rep(1,nrow(X1)),X1)
    X_sandwich <- chol2inv(chol( t(X)%*%X ) )%*% t(X)
    d_i <- function(df_ret,Bindex){
      df <- df_ret[Bindex,]
      y <- df$ERet
      coe <- X_sandwich %*% as.matrix(y)
      y_hat <- X %*% coe
      
      y - y_hat + coe[1,1]
    } 
    
    lapply(ds, function(x) d_i(Bindex=Bindex,df = x))
    
  }
  
  dt_list <- dt_calculator(ds = ds,Bindex = 1:obs_per_fund) 
  
  if (parallel) {
    omega_mu_hats <- parallel::mclapply(dt_list, function(x) omega_mu_hat(d = x,q = q,studentized = studentized),mc.cores = mc.cores ) %>% 
      bind_rows() %>% as_tibble()    
  }else{
    omega_mu_hats <- lapply(dt_list, function(x) omega_mu_hat(d = x,q = q,studentized = studentized) ) %>% 
      bind_rows() %>% as_tibble()
  }
  
  Bindex_list <- replicate(n = 1000, GenerateBootsIndex(n = TimeLength,q = q),simplify = FALSE)
  
  d_bar_b_function <- function(Bindex){
    dtb <- dt_calculator(ds = ds,Bindex = Bindex)
    d_bar_b <- sapply(dtb, function(x) mean(x) )
    d_bar_b
  }
  if (parallel) {
    d_bar_b <- parallel::mclapply(Bindex_list, function(x) d_bar_b_function(Bindex=x),mc.cores = mc.cores)     
  }else{
    d_bar_b <- lapply(Bindex_list, function(x) d_bar_b_function(Bindex=x)) 
  }
  
  
  d_bar_b <- matrix(d_bar_b %>% unlist() ,ncol = length(Bindex_list),byrow = FALSE)
  
  dim(d_bar_b)
  
  stats_b_set <- function(omega_mu_hats,d_bar_b,sub_index,obs_per_fund,studentized){
    # returns set of bootstrapped statistics distribution, given
    # inputs of the set of index of included funds' sub_index 
    omega_mu_hats_sub <- omega_mu_hats[sub_index,]
    d_bar_b_sub <- d_bar_b[sub_index,]
    if (studentized) {
      omega <- omega_mu_hats_sub$omega
      res <- sapply(1:ncol(d_bar_b_sub), function(x) sqrt(obs_per_fund) * max( 
        (d_bar_b_sub[,x] - omega_mu_hats_sub$d_bar + omega_mu_hats_sub$mu_hat)/omega 
      )
      )
      
    }else{
      res <- sapply(1:ncol(d_bar_b_sub), function(x) sqrt(obs_per_fund) * max(d_bar_b_sub[,x] - 
                                                                                omega_mu_hats_sub$d_bar + 
                                                                                omega_mu_hats_sub$mu_hat) )
    }
    res
  }
  
  rejections_list <- rejections_stepM_list <- list()
  for (tau in 1:length(targets)) {
    
    target <- targets[tau]
    
    # whole sample
    stats_b_dis <- stats_b_set(omega_mu_hats = omega_mu_hats,d_bar_b = d_bar_b,
                               sub_index = 1:nrow(omega_mu_hats),
                               obs_per_fund = obs_per_fund,
                               studentized = studentized)
    
    c_star <- max(0,quantile(x = stats_b_dis,probs = 1-target))
    
    
    Stats <- data.frame(statistic = omega_mu_hats$statistic,id = 1:length(omega_mu_hats$statistic))
    
    new_rejection <- Stats$id[Stats$statistic > c_star]
    rejections <- new_rejection
    
    
    # sub-sample
    sub_index <- setdiff(Stats$id, new_rejection)
    
    
    sub_stats <- Stats[Stats$id %in% sub_index,]
    
    sub_stats_b_dis <- stats_b_set(sub_index = sub_index,omega_mu_hats = omega_mu_hats,
                                   d_bar_b = d_bar_b,obs_per_fund = obs_per_fund,
                                   studentized = studentized)
    
    while (length(new_rejection) > 0) {
      c_s <- max(0,quantile(x = sub_stats_b_dis,probs = 1-target))
      new_rejection <- sub_stats$id[sub_stats$statistic > c_s]
      
      if (length(new_rejection) > 0) {
        #update for the sub-sample
        sub_index <- setdiff(sub_stats$id, new_rejection)
        rejections <- c(rejections,new_rejection)
        
        sub_stats <- sub_stats[sub_stats$id %in% sub_index,]
        
        sub_stats_b_dis <- stats_b_set(sub_index = sub_index,omega_mu_hats = omega_mu_hats,
                                       d_bar_b = d_bar_b,obs_per_fund = obs_per_fund,
                                       studentized = studentized)
        
      }
    }
    
    rejections_list[[tau]] <- rejections
  }
  list(rejection = rejections_list)
}



library("sandwich")
library(lmtest)

stepM <- function(ds,q,targets = seq(0.005,0.1,0.005),obs_per_fund,studentized = TRUE,parallel = TRUE){
  
  dt_calculator <- function(ds,Bindex){
    
    X1 <- ds[[1]][Bindex,] %>% select(-ERet) %>% as.matrix()
    X <- cbind(rep(1,nrow(X1)),X1)
    X_sandwich <- chol2inv(chol( t(X)%*%X ) )%*% t(X)
    d_i <- function(df_ret,Bindex){
      df <- df_ret[Bindex,]
      y <- df$ERet
      coe <- X_sandwich %*% as.matrix(y)
      y_hat <- X %*% coe
      
      y - y_hat + coe[1,1]
    } 
    
    lapply(ds, function(x) d_i(Bindex=Bindex,df = x))
    
  }
  
  dt_list <- dt_calculator(ds = ds,Bindex = 1:obs_per_fund) 
  if (parallel) {
    omega_mu_hats <- parallel::mclapply(dt_list, function(x) omega_mu_hat(d = x,q = q,studentized = studentized),mc.cores = mc.cores ) %>% 
      bind_rows() %>% as_tibble()
  }else{
    omega_mu_hats <- lapply(dt_list, function(x) omega_mu_hat(d = x,q = q,studentized = studentized)) %>% 
      bind_rows() %>% as_tibble()
  }
  
  Bindex_list <- replicate(n = 1000, GenerateBootsIndex(n = TimeLength,q = q),simplify = FALSE)
  
  d_bar_omega_b_function <- function(Bindex){
    dtb <- dt_calculator(ds = ds,Bindex = Bindex)
    
    d_bar_b <- sapply(dtb, function(x) mean(x) )
    omega_b <- sapply(dtb, function(x) omega_hat(d = x,q = q) )
    data.frame(d_bar_b,omega_b)
  }
  if (parallel) {
    d_bar_omega_b <- parallel::mclapply(Bindex_list, function(x) d_bar_omega_b_function(Bindex=x),mc.cores = mc.cores) 
  }else{
    d_bar_omega_b <- lapply(Bindex_list, function(x) d_bar_omega_b_function(Bindex=x)) 
  }
  
  d_bar_b <- lapply(d_bar_omega_b, function(x) x$d_bar_b)
  d_bar_b <- matrix(d_bar_b %>% unlist() ,ncol = length(Bindex_list),byrow = FALSE)
  
  dim(d_bar_b)
  # for step M
  omega_b <- sapply(d_bar_omega_b, function(x) x$omega_b)
  
  stats_b_set <- function(omega_mu_hats,d_bar_b,omega_b,sub_index,obs_per_fund){
    # returns set of bootstrapped statistics distribution, given
    # inputs of the set of index of included funds' sub_index 
    omega_mu_hats_sub <- omega_mu_hats[sub_index,]
    d_bar_b_sub <- d_bar_b[sub_index,]
    omega_b_sub <- omega_b[sub_index,]
    
    sapply(1:ncol(d_bar_b_sub), function(x) max( (d_bar_b_sub[,x] - omega_mu_hats_sub$d_bar)/omega_b_sub[,x] ) )
  }
  
  rejections_list <- rejections_stepM_list <- list()
  for (tau in 1:length(targets)) {
    
    target <- targets[tau]
    
    # whole sample
    stats_b_dis <- stats_b_set(omega_mu_hats = omega_mu_hats,omega_b = omega_b,d_bar_b = d_bar_b,
                               sub_index = 1:nrow(omega_mu_hats),
                               obs_per_fund = obs_per_fund)
    
    c_star <- quantile(x = stats_b_dis,probs = 1-target)
    
    Stats <- data.frame(statistic = omega_mu_hats$d_bar/omega_mu_hats$omega,id = 1:length(omega_mu_hats$omega))
    
    new_rejection <- Stats$id[Stats$statistic > c_star]
    rejections <- new_rejection
    
    
    # sub-sample
    sub_index <- setdiff(Stats$id, new_rejection)
    length(sub_index)
    
    sub_stats <- Stats[Stats$id %in% sub_index,]
    dim(sub_stats)
    
    sub_stats_b_dis <- stats_b_set(sub_index = sub_index,omega_b = omega_b,omega_mu_hats = omega_mu_hats,
                                   d_bar_b = d_bar_b,obs_per_fund = obs_per_fund)
    
    while (length(new_rejection) > 0) {
      c_s <- quantile(x = sub_stats_b_dis,probs = 1-target)
      new_rejection <- sub_stats$id[sub_stats$statistic > c_s]
      
      if (length(new_rejection) > 0) {
        #update for the sub-sample
        sub_index <- setdiff(sub_stats$id, new_rejection)
        rejections <- c(rejections,new_rejection)
        
        sub_stats <- sub_stats[sub_stats$id %in% sub_index,]
        
        sub_stats_b_dis <- stats_b_set(sub_index = sub_index,omega_mu_hats = omega_mu_hats,
                                       d_bar_b = d_bar_b,omega_b = omega_b,obs_per_fund = obs_per_fund)
        
      }
    }
    
    sub_stats_b_dis_stepM <- stats_b_set(sub_index = sub_index,omega_b = omega_b,omega_mu_hats = omega_mu_hats,
                                         d_bar_b = d_bar_b,obs_per_fund = obs_per_fund)
    rejections_list[[tau]] <- rejections
  }
  list(rejection = rejections_list)
}
