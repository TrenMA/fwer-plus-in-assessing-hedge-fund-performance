ma_unsmooth <- function(ret,H,optimal_lags = FALSE){
  # ret is a fund's return time series 
  # H is number of lags
  
  # we assume that ret_unsmooth = ret_mean + eta
  # And thus, ret = MA(H) of ret_unsmooth = sum(theta*ret_unsmooth) 
  #               = ret_mean + sum(theta * eta)
  
  ret_demean <- ret - mean(ret)
  
  ma <- tryCatch(arima(x = ret_demean,order = c(0,0,H),
                       include.mean = FALSE,transform.pars = TRUE,
                       init = rep(-0.2,H), method = "ML",
                       optim.control = list(maxit = 5000)), 
                 error = function(e) NULL)
  
  if (is.null(ma)) {
    
    theta0 <- 1
    theta_lags <- rep(0,H)
    ma_null <- 1
    
  } else {
    zeros <- H - length(ma$coef)
    
    # add ma lags that are zeros
    theta_lags <- c(ma$coef,rep(0,zeros)) # it is (ma1, ma2, ..., maH)
    
    # normalized
    thetas_sum <- 1 + sum(theta_lags)
    theta_lags <- theta_lags/thetas_sum
    theta0 <- 1/thetas_sum
    
    ma_null <- 0
  }
  
  # check if the thetas are too big or too small
  
  thetas <- c(theta0,theta_lags)
  
  check <- sum(thetas > 1.25 | thetas < -0.45)
  
  if (check > 0) {
    theta0 <- 1
    theta_lags <- rep(0,H)
  }
  
  # eta
  eta <- c()
  
  # initiate first error terms in MA model: first H current errors are those
  # from demean smoothed return
  
  eta[1:H] <- ret_demean[1:H]
  
  # from H+1,  
  for (t in (H+1):length(ret_demean)) {
    ret_demean_t <- ret_demean[t] 
    
    # current error term
    eta[t] <- (ret_demean_t - sum(theta_lags*eta[(t-1):(t-H)]) )/theta0
    
  }
  # correct for eta so that its mean is zero
  # and add it into mean return to recover economic/unsmoothed return
  
  ret_unsmooth <- (eta - mean(eta)) + mean(ret)  
  
  ret_unsmooth
}
