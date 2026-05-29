library(dplyr)
source('Unsmooth return function.R')

df <- readRDS('Fake Data for Running Replication Code.RDS')  

dfs <- df %>% 
  select(Date,Ret,Name) %>%
  group_split(Name)

length(dfs)

res <- list()

for (i in 1:length(dfs)) {
  df_temp <- dfs[[i]] %>% arrange(Date)
  Ret_unsmoothed <- ma_unsmooth(ret = df_temp$Ret,H = 3)
  df_temp$Ret <- Ret_unsmoothed
  res[[i]] <- df_temp
}

saveRDS(res %>% bind_rows(), file = 'Unsmoothed Return for Simulation.RDS')

