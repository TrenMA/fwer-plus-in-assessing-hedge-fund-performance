# FWER procedure for two-sided tests: fwer_twosided

# This code follows Jun Chen's GitHub repository for the FWER component
# of the CAMT package:
# https://github.com/jchen1981/CAMT
# Modifications:
# - The code has been adapted to return results for a set of FWER targets
#   (parameter alpha) in a single run, rather than a single target (see lines 217–223).
# - When solving for covariate coefficients, if the unconstrained solver
#   `nlm` fails, the constrained solver `nlminb` is used as a fallback (see lines 179–188).

solvek.fwer <- function(pvals, pi0, init = 0.25, nlm.iter = 10, pvals.cutoff = 1e-15) {
  
  pvals [pvals  <= pvals.cutoff] <- pvals.cutoff
  pvals [pvals >= 1 - pvals.cutoff] <- 1 - pvals.cutoff
  
  fun <- function(k) {
    - sum(log(pi0 + (1 - pi0) * k * pvals ^ (k - 1))) 
  }
  
  err <- try(k <- nlminb(init, fun, control = list(iter.max = nlm.iter), 
                         lower = 0, upper = 1)$par, silent = TRUE)
  if (class(err) == 'try-error') {
    k <- 0.25
  }
  
  return(k)
}

fwer_twosided <- function(pvals, pi0.var, data = data, alpha, 
                          pvals.cutoff = 1e-15, k.fix = NULL,
                          EM.paras = list(iterlim = 500, tol = 1e-5, 
                                          k.init = NULL, pi0.init = NULL, 
                                          nlm.iter = 50), 
                          trace = FALSE, return.model.matrix = FALSE) {
  
  m0 <- length(pvals)
  
  if (sum(is.na(pvals)) > 0) stop("The pvalues contain NA! Please remove!\n")
  
  m <- m0 <- length(pvals)
  
  qvals <- qvalue(pvals, pi0.method = "bootstrap")
  gamma <- qvals$lambda[which.min(abs(qvals$pi0.lambda - qvals$pi0))]
  Y <- (pvals > gamma)
  b0 <- (1 - gamma) * Y + gamma * (1 - Y)
  
  # Initialize pi0
  if (is.null(EM.paras$pi0.init)) {
    if (trace) cat('Estimating initial pi0 values ...\n')
    pi0.init <- qvals$pi0
    if (pi0.init >= 0.99)  pi0.init <- 0.99
    EM.paras$pi0.init <- pi0.init
  } else {
    pi0.init <- EM.paras$pi0.init
  }
  
  # Initialize k
  if (is.null(EM.paras$k.init)) {
    if (trace) cat('Estimating initial k values ...\n')
    k.init <- solvek.fwer(pvals, pi0.init, pvals.cutoff = pvals.cutoff)
    EM.paras$k.init <- k.init
  } else {
    k.init <- EM.paras$k.init
  }
  
  if (is.null(EM.paras$iterlim)) {
    iterlim <- 50
    EM.paras$iterlim <- iterlim
  } else {
    iterlim <- EM.paras$iterlim
  }
  
  if (is.null(EM.paras$nlm.iter)) {
    nlm.iter <- 5
    EM.paras$nlm.iter <- nlm.iter
  } else {
    nlm.iter <- EM.paras$nlm.iter
  }
  
  if (is.null(EM.paras$tol)) {
    tol <- 1e-5
    EM.paras$tol <- tol
  } else {
    tol <- EM.paras$tol
  }
  
  if (is.null(pi0.var)) {
    pi0.var <- matrix(rep(1, m))
  } else {
    if (length(intersect(class(pi0.var), c('matrix', 'numeric', 'character', 'factor', 'data.frame', 'formula'))) == 0) {
      stop('Currently do not support the class of pi0.var!\n')
    }
    if (length(intersect('matrix', class(pi0.var))) != 0) {
      if (mode(pi0.var) == 'numeric') {
        if (sum(pi0.var[, 1] == 1) == m0) {
          pi0.var <- pi0.var
        } else {
          pi0.var <- cbind(rep(1, m0), pi0.var)
        }
        
      } else {
        warning('The "pi0.var" matrix contains non-numeric values! Will treat it as a data frame!\n')
        pi0.var <- model.matrix(~., data.frame(pi0.var))
      }
    } else {
      if (length(intersect(c('character', 'factor'), class(pi0.var))) != 0 ) {
        pi0.var <- model.matrix(~., data.frame(pi0.var))
      }  else {
        if (length(intersect(c('data.frame'), class(pi0.var))) != 0) {
          pi0.var <- model.matrix(~., pi0.var)
        }  else {
          if (length(intersect('numeric', class(pi0.var))) != 0) {
            pi0.var <- cbind(rep(1, m0), pi0.var)
            
          } else {
            if (length(intersect(c('formula'), class(pi0.var))) != 0) {
              if (is.null(data)) {
                stop('"pi0.var" is a formula and "data" should be a data frame!\n')
              } else {
                pi0.var <- model.matrix(pi0.var, data)
              }
              
            }
          }
        }
      }
    }
  }
  
  # Threshold the p-value for more robustness
  pvals[pvals <= pvals.cutoff] <- pvals.cutoff
  pvals[pvals >= 1 - pvals.cutoff] <- 1 - pvals.cutoff
  
  len.theta <- ncol(pi0.var)
  theta <- c(binomial()$linkfun(pi0.init), rep(0, len.theta - 1))
  
  if (is.null(k.fix)) {
    k <- k.init
  } else {
    k <- k.fix
  }
  
  fun1 <- function(theta, q0) {
    q1 <- 1 - q0
    exp.eta.theta <- exp(as.vector(pi0.var %*% theta))
    pi0 <- exp.eta.theta / (1 + exp.eta.theta)
    - sum(q0 * log(pi0) + q1 * log(1 - pi0))
  }
  fun2 <- function(k, q0) {
    q1 <- 1 - q0  
    gammak <- gamma ^ k
    - sum(q1 * (Y * log(1 - gammak) + k * (1 - Y) * log(gamma)))
  }
  
  iter <- 1
  loglik1 <- 1
  
  if (trace) cat('Run EM algorithm...\n')
  while (iter <= iterlim) {
    
    if (trace) cat('.')
    # E step
    exp.eta.theta <- exp(as.vector(pi0.var %*% theta))
    pi0 <- exp.eta.theta / (1 + exp.eta.theta)
    
    gammak <- gamma ^ k
    b1 <- (1 - gammak) * Y + gammak * (1 - Y)
    pi0b0 <- pi0 * b0
    pib <- pi0b0 + (1 - pi0) * b1
    q0 <- pi0b0 / pib
    
    loglik2 <- sum(log(pib))
    
    if ((abs((loglik2 - loglik1) / loglik1) < tol) | (iter >= iterlim)) {
      break
    } else {
      
      # If nlm (un-constraint) counters errors then use nlminb (constraint)
      theta_temp <- tryCatch(suppressWarnings(nlm(fun1, theta, q0 = q0, iterlim = nlm.iter)),error=function(e) NA)
      
      if (length(theta_temp) <= 1) {
        if (is.na(theta_temp)) {
          theta <- suppressWarnings( nlminb(start = theta,objective = fun1,q0 = q0))$par
        }
      }else{
        theta <- theta_temp$estimate  
      }
      
      if (is.null(k.fix)) {
        k <- suppressWarnings(nlminb(k, fun2, q0 = q0, control = list(iter.max = nlm.iter), 
                                     lower = 0, upper = 1))$par
      }
    }
    
    loglik1 <- loglik2
    iter <- iter + 1
  }
  
  if (iter == iterlim) warning('Maximum number of iterations reached!')
  
  
  fwer.obj <- list(call = match.call(), pvals = pvals, pi0.var = pi0.var, alpha = alpha,
                        pi0 = pi0, k = k, pi0.coef = theta, gamma = gamma,
                        loglik = loglik2, EM.paras = EM.paras, EM.iter = iter) 
  

  fwer.rejection <- function (fwer.obj, alpha){
    
    pvals <- fwer.obj$pvals
    pi0 <- fwer.obj$pi0
    gamma <- fwer.obj$gamma
    k <- fwer.obj$k
    Y <- (pvals > gamma)
    pi1 <- 1 - pi0
    
    # Enable results for a set of fwer targets each run:     
    rej <- list()
    for (i in 1:length(alpha)) {
      tau <- k * (alpha[i] * (1 - gamma))^(k - 1) * (sum(Y * (pi1/pi0)^(1/(1 - k))))^(1 - k)
      t <- (pi1 * k/(pi0 * tau))^(1/(1 - k))
      rej[[i]] <- which(pvals < pmin(t, gamma))
    }
    
    names(rej) <- alpha
    return(rej)
  }
  
  fwer.obj$rejection <- fwer.rejection(fwer.obj, alpha)
  
  if (return.model.matrix) {
    fwer.obj$pi0.var <- fwer.obj$pi0.var
  }
  
  return(fwer.obj)
}
