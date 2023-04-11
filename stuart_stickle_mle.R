library(tidyverse)
library(tmvtnorm)

dat <- read.csv("stuart_stickle_imp_raw.csv")
dat2 <- read.csv("stickleback_fossil_vs_extant.csv")

dat_small <- dat %>% filter(imputation_number == 1)

mvec <- function(x_m,x_f,theta_m,theta_f,kappa){
  t <- dt$time
  sex <- model.matrix(~ 0 + gender,dt)
  mvec <- (sex %*% c(x_f,x_m))*exp(-kappa*t) +
    (sex %*% c(theta_f,theta_m))*(1 - exp(-kappa*t))
  return(mvec)
}

vcovmat <- function(kappa,sigma){
  t <- dt$time
  sex <- model.matrix(~ 0 + gender,dt)
  vars <- diag(sigma^2/(2*kappa) * (1 - exp(-2*kappa*t)))
  time_diff <- abs(matrix(t,nrow=nrow(dt),ncol=nrow(dt),byrow=TRUE) -
                   matrix(t,nrow=nrow(dt),ncol=nrow(dt)))
  time_add <- matrix(t,nrow=nrow(dt),ncol=nrow(dt),byrow=TRUE) +
              matrix(t,nrow=nrow(dt),ncol=nrow(dt))
  covs <- sigma^2/(2*kappa)*(exp(-2*kappa*time_diff) * exp(-2*kappa*time_add))
  vcovmat <- vars+covs
  return(vcovmat)
}

minusll <- function(x_m,x_f,theta_m,theta_f,kappa,sigma){
  mm <- mvec(x_m,x_f,theta_m,theta_f,kappa)
  vv <- vcovmat(kappa,sigma)
  y <- dt$stl
  return(-dmvnorm(y,mean=mm,sigma=vv,log=TRUE) %>% sum)
}

final <- NULL

for(i in 1:5){
  print(paste0("Running Imputation ",i))
  dt <- dat %>% filter(imputation_number == i)
  test <- mle2(minusll,
               start=list(x_m=66,x_f=61,theta_m=66,theta_f=61,kappa=0.44,sigma=11))
  savs <- c(coef(test)[3:4],vcov(test)[3:4,3:4] %>% as.vector)
  final <- rbind(final,savs)
}
