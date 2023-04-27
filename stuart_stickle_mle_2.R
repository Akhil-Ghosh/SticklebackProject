library(tidyverse)
library(tmvtnorm)

dat <- read.csv("stuart_stickle_imp_100_raw.csv")
M <- max(dat$imputation_number)

for(m in 1:M){
  print(paste0("Running Imputation ",m))
  tmp <- dat %>% filter(imputation_number == m)

  T <- max(tmp$time)
  gentime_summary <- tmp %>%
    group_by(gender,time) %>%
    summarise(mu=mean(stl)) %>%
    group_by(gender) %>%
    mutate(theta = mean(mu)) %>%
    mutate(u = mu - theta) %>%
    ungroup
}
