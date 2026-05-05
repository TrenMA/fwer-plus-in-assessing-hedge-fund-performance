
include_spa <- TRUE; q_spa <- 0.9
include_stepM <- TRUE; q_stepM <- 0.9

rho <- 0
n_funds <- 1000
TimeLength <- 36
b_case <- 2

dfs <- res <- list()
k <- 1

for (alpha in fund_alpha[1:3]) {
  
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
              fwer_stepM = mean(stepM_fwe),
              power_stepM = mean(stepM_detected_rate),
              fwer_spa = mean(spa_fwe),
              power_spa = mean(spa_detected_rate),
              pi0 = mean(pi0),
              pip = mean(pip),
              pim = mean(pim)) %>% filter(target <= 0.1)
  
  k <- k + 1
}

df <- dfs[[1]]

plot_fwer_steps <- function(df, main) {
  
  df <- df * 100
  
  plot(df$fwer ~ df$target, type = "n", col = "blue", xlim = c(0, 10),
       ylim = c(0, 10), ylab = "Actual FWER (%)", xlab = "FWER target (%)", 
       main = main, family = "A")
  
  # Transparent grid lines
  grid_col <- rgb(0.5, 0.5, 0.5, alpha = 0.1)  # light gray with transparency
  abline(h = seq(0,10,5), col = grid_col, lty = 1)
  abline(v = seq(0,10,5), col = grid_col, lty = 1)
  
  
  colors <- c('blue', 'red3', 'green4')
  
  lines(x = df$target, y = df$fwer, col = colors[1])
  lines(x = df$target, y = df$fwer_spa, col = colors[2])
  lines(x = df$target, y = df$fwer_stepM, col = colors[3])
  
  lines(x = c(0, 100), y = c(0, 100), col = "black", lty = 2)
  
  leg <- list(
    bquote(FWER^"+"),
    "StepSPA",
    "StepM"
  )
  
  legend('topleft', legend = leg, col = colors, lty = 1, ncol = 1, cex = 0.8)
}


plot_power_steps <- function(df, main,power_max) {
  df <- df * 100
  
  plot(df$power ~ df$target, type = "n", col = "blue", xlim = c(0, 10),
       ylim = c(0, max(df$power)), ylab = "Power (%)", xlab = "FWER target (%)", 
       main = main, family = "A")
  
  # Transparent grid lines
  grid_col <- rgb(0.5, 0.5, 0.5, alpha = 0.1)  # light gray with transparency
  abline(h = seq(0,max(df$power),length = 5), col = grid_col, lty = 1)
  abline(v = seq(0,10,5), col = grid_col, lty = 1)
  
  colors <- c('blue', 'red3', 'green4')
  
  lines(x = df$target, y = df$power, col = colors[1])
  lines(x = df$target, y = df$power_spa, col = colors[2])
  lines(x = df$target, y = df$power_stepM, col = colors[3])
  
  #lines(x = c(0, 100), y = c(0, 100), col = "black", lty = 2)
  
  leg <- list(
    bquote(FWER^"+"),
    "StepSPA",
    "StepM"
  )
  
  legend('topleft', legend = leg, col = colors, lty = 1, ncol = 1, cex = 0.8)
}

windowsFonts(A = windowsFont("Times New Roman"))
file_name <- paste('Performance comparison among three methods under varying alpha.png',sep = '',collapse = '')
png(filename = file_name,width = 3520 ,height = 2400,res = 450)
par(mfrow = c(2,3),mar=c(4,4,3,1),family = "A")
m1 <- bquote(alpha == 0.5)
m2 <- bquote(alpha == 1.0)
m3 <- bquote(alpha == 1.5)

mains <- c(m1,m2,m3)

for (k in 1:3) {
  main <- mains[k]
  plot_fwer_steps(df = dfs[[k]],main = main)
}
dff <- dfs %>% bind_rows()
power_max <- max(dff$power)*100
for (k in 1:3) {
  main <- mains[k]
  plot_power_steps(df = dfs[[k]],main = main,power_max = power_max)
}
dev.off()

