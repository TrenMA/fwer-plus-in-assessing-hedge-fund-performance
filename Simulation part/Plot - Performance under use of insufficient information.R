
#######################################################################
# noise

include_noise <- '_Noise_Covariates'

n_funds <- 1000
TimeLength <- 36

rho <- 0
l <- 1
dfs0 <- list()

alpha <- fund_alpha[2]
k <- 1
dfs <- list()
for (b_case in 1:3) {
  
  setwd(folder_simulation_results)
  file_name <-  paste0('Simualtion_with_b_',b_case,'_Zrho_',rho,
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
              fwer_noise_0.5 = mean(fwer_noise_0.5),
              power_noise_0.5 = mean(fwer_noise_detected_rate_0.5),
              fwer_noise_1 = mean(fwer_noise_1),
              power_noise_1 = mean(fwer_noise_detected_rate_1),
              fwer_noise_1.5 = mean(fwer_noise_1.5),
              power_noise_1.5 = mean(fwer_noise_detected_rate_1.5),
              fwer_1noise = mean(fwe_1noise_4.5),
              fwer_2noise = mean(fwe_2noise_4.5),
              fwer_pca = mean(fwe_pca_4.5),
              power_1noise = mean(fwe_1noise_detected_rate_4.5),
              power_2noise = mean(fwe_2noise_detected_rate_4.5),
              power_pca = mean(fwe_pca_detected_rate_4.5),
              pi0 = mean(pi0),
              pip = mean(pip),
              pim = mean(pim)) %>% filter(target <= 0.1)
  
  k <- k + 1
}

df <- dfs[[1]] * 100
plot_fwer_noise <- function(df, main) {
  df <- df*100
  plot(df$fwer ~ df$target, type = "n", col = "blue", xlim = c(0, 10 ),
       ylim = c(0, 10 ), ylab = "Actual FWER (%)", xlab = "FWER target (%)", 
       main = main, family = "A")
  
  # Transparent grid lines
  grid_col <- rgb(0.5, 0.5, 0.5, alpha = 0.1)  # light gray with transparency
  abline(h = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  abline(v = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  
  lines(x = c(0, 100), y = c(0, 100), col = "black", lty = 2)
  colors <- c( 'red3', 'green4','blue4','black','violet','green','blue')
  
  
  
  lines(x = df$target, y = df$fwer_noise_0.5, col = colors[1])
  lines(x = df$target, y = df$fwer_noise_1, col = colors[2])
  lines(x = df$target, y = df$fwer_noise_1.5, col = colors[3])
  lines(x = df$target, y = df$fwer, col = colors[4],lty = 2)
  lines(x = df$target, y = df$fwer1, col = colors[5],lty = 2)
  lines(x = df$target, y = df$fwer_un, col = colors[6],lty = 2)
  
  
  leg <- list(
    bquote(Z^"*" ~ "," ~ sigma == 0.5),
    bquote(Z^"*" ~ "," ~ sigma == 1.0),
    bquote(Z^"*" ~ "," ~ sigma == 1.5),
    bquote(Z == paste("(", U, ",", V, ")")),
    bquote(Z == U),
    "Noise"
  )
  
  legend('topleft', legend = leg, col = colors, lty = c(1,1,1,2,2,2), ncol = 1, cex = 0.8)
}

plot_power_noise <- function(df, main) {
  
  df <- df*100
  
  plot(df$power ~ df$target, type = "n", col = "blue", xlim = c(0, 10 ),
       ylim = c(0, max(df$power)+0.5), ylab = "Power (%)", xlab = "FWER target (%)", 
       main = main, family = "A")
  
  # Transparent grid lines
  grid_col <- rgb(0.5, 0.5, 0.5, alpha = 0.1)  # light gray with transparency
  abline(h = seq(0,max(df$power)+0.5,length = 5), col = grid_col, lty = 1)
  abline(v = seq(0,10 ,length = 5), col = grid_col, lty = 1)
  
  
  colors <- c( 'red3', 'green4','blue4','black','violet','green','blue')
  
  lines(x = df$target, y = df$power_noise_0.5, col = colors[1])
  lines(x = df$target, y = df$power_noise_1, col = colors[2])
  lines(x = df$target, y = df$power_noise_1.5, col = colors[3])
  lines(x = df$target, y = df$power, col = colors[4],lty = 2)
  lines(x = df$target, y = df$power1, col = colors[5],lty = 2)
  lines(x = df$target, y = df$power_un, col = colors[6],lty = 2)
  
  leg <- list(
    bquote(Z^"*" ~ "," ~ sigma == 0.5),
    bquote(Z^"*" ~ "," ~ sigma == 1.0),
    bquote(Z^"*" ~ "," ~ sigma == 1.5),
    bquote(Z == paste("(", U, ",", V, ")")),
    bquote(Z == U),
    "Noise"
  )
  
  legend('topleft', legend = leg, col = colors, lty = c(1,1,1,2,2,2), ncol = 1, cex = 0.8)
}


windowsFonts(A = windowsFont("Times New Roman"))

file_name <- paste('Performance of fwer plus under use of insufficient information.png',sep = '',collapse = '')
png(filename = file_name,width = 3520 ,height = 2400,res = 450)
par(mfrow = c(2,3),mar=c(4,4,3,1),family = "A")

m1 <- bquote(paste("Weakly informative Z"))
m2 <- bquote(paste("Moderately informative Z"))
m3 <- bquote(paste("Strongly informative Z"))

mains <- c(m1,m2,m3)

for (k in 1:3) {
  
  main <- mains[k]
  plot_fwer_noise(df = dfs[[k]],main = main)
  
}

for (k in 1:3) {
  main <- mains[k]
  plot_power_noise(df = dfs[[k]],main = main)
  
}

dev.off()
