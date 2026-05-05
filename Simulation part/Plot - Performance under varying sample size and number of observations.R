
################################################################################
# Plot for performance of fwer+ under varying sample size and number of observations
################################################################################

b_case <- 2
rho <- 0
l <- 1
dfs0 <- list()
l <- 1
alpha <- 1

for (TimeLength in c(24,36,48)) {
  k <- 1
  dfs <- list()
  for (n_funds in c(500,1000,2000,3000)) {
    
    setwd(folder_simulation_results)
    file_name <- paste0('Simualtion_with_b_',b_case,'_Zrho_',rho,
                        '_obs_per_fund_',TimeLength,'_n_funds_',n_funds,
                        '_alpha_',alpha,'_StepSPA_',include_spa,
                        '_stepM_',include_stepM,include_noise,'.RDS')
    
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
  dfs0[[l]] <- dfs; l <- l + 1
}


plot_fwer_multiple_n_T <- function(dfs, main) {
  
  df <- dfs[[1]] * 100
  
  plot(df$fwer ~ df$target, type = "n", col = "blue", xlim = c(0, 10 ),
       ylim = c(0, 10 ), ylab = "Actual FWER (%)", xlab = "FWER target (%)", 
       main = main, family = "A")
  
  # Transparent grid lines
  grid_col <- rgb(0.5, 0.5, 0.5, alpha = 0.1)  # light gray with transparency
  abline(h = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  abline(v = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  
  lines(x = c(0, 100), y = c(0, 100), col = "black", lty = 2)
  colors <- c( 'red3', 'green4','blue','black','violet')
  
  for (i in 1:length(dfs)) {
    df <- dfs[[i]] *100
    lines(x = df$target, y = df$fwer, col = colors[i])    
  }
  
  
  leg <- list(
    bquote(n == 500),
    bquote(n == 1000),
    bquote(n == 2000),
    bquote(n == 3000)
  )
  
  legend('topleft', legend = leg, col = colors, lty = 1, ncol = 1, cex = 0.8)
}


plot_power_multiple_n_T <- function(dfs, main) {
  
  df <- dfs %>% bind_rows() * 100
  
  plot(df$power ~ df$target, type = "n", col = "blue", xlim = c(0, 10 ),
       ylim = c(0, max(df$power)), ylab = "Power (%)", xlab = "FWER target (%)", 
       main = main, family = "A")
  
  # Transparent grid lines
  grid_col <- rgb(0.5, 0.5, 0.5, alpha = 0.1)  # light gray with transparency
  abline(h = seq(0,max(df$power),length = 5), col = grid_col, lty = 1)
  abline(v = seq(0,10,length = 5), col = grid_col, lty = 1)
  
  colors <- c( 'red3', 'green4','blue','black','violet')
  
  for (i in 1:length(dfs)) {
    df <- dfs[[i]] *100
    lines(x = df$target, y = df$power, col = colors[i])    
  }
  
  
  leg <- list(
    bquote(n == 500),
    bquote(n == 1000),
    bquote(n == 2000),
    bquote(n == 3000)
  )
  
  legend('topleft', legend = leg, col = colors, lty = 1, ncol = 1, cex = 0.8)
}


windowsFonts(A = windowsFont("Times New Roman"))

file_name <- paste('Performance of fwer plus under varying sample size and number of observations.png',sep = '',collapse = '')
png(filename = file_name,width = 3520 ,height = 2400,res = 450)
par(mfrow = c(2,3),mar=c(4,4,3,1),family = "A")

m1 <- bquote(T == 24)
m2 <- bquote(T == 36)
m3 <- bquote(T == 48)

mains <- c(m1,m2,m3)
for (k in 1:3) {
  
  main <- mains[k]
  plot_fwer_multiple_n_T(dfs = dfs0[[k]],main = main)
  
}

for (k in 1:3) {
  main <- mains[k]
  plot_power_multiple_n_T(dfs = dfs0[[k]],main = main)
  
}

dev.off()
