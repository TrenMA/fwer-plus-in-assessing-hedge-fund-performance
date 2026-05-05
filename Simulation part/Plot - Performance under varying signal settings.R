################################################################################
# Performance of fwer+ under varying signal settings: 
# independent and correlated covariates 
# all alpha cases, all beta cases
################################################################################

include_noise <- NULL
n_funds <- 1000
TimeLength <- 36
fund_alpha <- c(0.5,1,1.5)
include_stepM <- include_spa <- FALSE;

# collect simulation results for independent covariates

l <- 1
dfs0 <- list()
rho <- 0

for (alpha in fund_alpha) {
  k <- 1
  dfs <- list()
  for (b_case in 1:3) {
    
    setwd(folder_simulation_results)
    file_name <- paste('Simualtion_with_b_',b_case,'_Zrho_',rho,
                       '_obs_per_fund_',TimeLength,'_n_funds_',n_funds,
                       '_alpha_',alpha,'_StepSPA_',include_spa,
                       '_stepM_',include_stepM,include_noise,
                       '.RDS',sep = '',collapse = '')
    
    df <- readRDS(file = file_name) %>% bind_rows()  
    dfs[[k]] <- df %>% group_by(target) %>% 
      summarise(fwer = mean(fwe),
                power = mean(fwe_detected_rate),
                fwer1 = mean(fwe1),
                power1 = mean(fwe1_detected_rate),
                fwer_un = mean(fwe_uninfo),
                power_un = mean(fwe_uninfo_detected_rate),
                pi0 = mean(pi0),
                pip = mean(pip),
                pim = mean(pim))%>% filter(target <= 0.1)
    
    k <- k + 1
  }
  dfs0[[l]] <- dfs; l <- l + 1
}

# collect simulation results for correlated covariates

rho <- 0.5
l <- 1
dfs05 <- list()

for (alpha in fund_alpha) {
  setwd(folder_simulation_results)
  
  dfs <- list()
  k <- 1
  for (b_case in 1:3) {
    file_name <- paste('Simualtion_with_b_',b_case,'_Zrho_',rho,
                       '_obs_per_fund_',TimeLength,'_n_funds_',n_funds,
                       '_alpha_',alpha,'_StepSPA_',include_spa,
                       '_stepM_',include_stepM,include_noise,
                       '.RDS',sep = '',collapse = '')
    df <- readRDS(file = file_name) %>% bind_rows()  
    dfs[[k]] <- df %>% group_by(target) %>% 
      summarise(fwer = mean(fwe),
                power = mean(fwe_detected_rate),
                fwer1 = mean(fwe1),
                power1 = mean(fwe1_detected_rate),
                fwer_un = mean(fwe_uninfo),
                power_un = mean(fwe_uninfo_detected_rate),
                pi0 = mean(pi0),
                pip = mean(pip),
                pim = mean(pim)) %>% filter(target <= 0.1)
    
    k <- k + 1
  }
  dfs05[[l]] <- dfs; l <- l + 1
}


plot_fwer <- function(df0,df05, main) {
  
  df <- df0[[1]] * 100
  
  plot(df$fwer ~ df$target, type = "n", col = "blue", xlim = c(0, 10 ),
       ylim = c(0, 10 ), ylab = "Actual FWER (%)", xlab = "FWER target (%)", 
       main = main, family = "A")
  
  # Transparent grid lines
  grid_col <- rgb(0.5, 0.5, 0.5, alpha = 0.1)  # light gray with transparency
  abline(h = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  abline(v = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  
  lines(x = c(0, 100), y = c(0, 100), col = "black", lty = 2)
  colors <- c( 'red3', 'green4','blue')
  
  for (i in 1:3) {
    df <- df0[[i]] *100
    lines(x = df$target, y = df$fwer, col = colors[i])    
  }
  
  for (i in 1:3) {
    df <- df05[[i]] *100
    lines(x = df$target, y = df$fwer, col = colors[i],lty = 3)    
  }
  
  leg <- list(
    'Weakly informative Z',
    "Moderately informative Z",
    "Strongly informative Z"
  )
  
  legend('topleft', legend = leg, col = colors, lty = 1, ncol = 1, cex = 0.8)
}

plot_power <- function(df0,df05, main) {
  
  df1 <- df0 %>% bind_rows()
  df2 <- df05 %>% bind_rows()
  
  df <- bind_rows(df1,df2) * 100
  
  plot(df$fwer ~ df$target, type = "n", col = "blue", xlim = c(0, 10 ),
       ylim = c(0, max(df$power)), ylab = "Power (%)", xlab = "FWER target (%)", 
       main = main, family = "A")
  
  # Transparent grid lines
  grid_col <- rgb(0.5, 0.5, 0.5, alpha = 0.1)  # light gray with transparency
  abline(h = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  abline(v = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  
  colors <- c( 'red3', 'green4','blue')
  
  for (i in 1:3) {
    df <- df0[[i]] *100
    lines(x = df$target, y = df$power, col = colors[i])    
  }
  
  for (i in 1:3) {
    df <- df05[[i]] *100
    lines(x = df$target, y = df$power, col = colors[i],lty = 3)    
  }
  
  leg <- list(
    'Weakly informative Z',
    "Moderately informative Z",
    "Strongly informative Z"
  )
  
  legend('bottomright', legend = leg, col = colors, lty = 1, ncol = 1, cex = 0.8)
}

# plot figure

windowsFonts(A = windowsFont("Times New Roman"))
file_name <- paste('Performance of fwer plus under varying signal settings.png',sep = '',collapse = '')
png(filename = file_name,width = 3520 ,height = 2400,res = 450)
par(mfrow = c(2,3),mar=c(4,4,3,1),family = "A")

m1 <- bquote(alpha == 0.5)
m2 <- bquote(alpha == 1.0)
m3 <- bquote(alpha == 1.5)

mains <- c(m1,m2,m3)
for (k in 1:3) {
  
  main <- mains[k]
  plot_fwer(df0=dfs0[[k]],df05 = dfs05[[k]],main = main)
  
}

for (k in 1:3) {
  main <- mains[k]
  plot_power(df0=dfs0[[k]],df05 = dfs05[[k]],main = main)
  
}

dev.off()