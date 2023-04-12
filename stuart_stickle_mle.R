library(tidyverse)
library(tmvtnorm)
library(bbmle)

dat <- read.csv("stuart_stickle_imp_raw.csv")
M <- max(dat$imputation_number)

minusll <- function(theta_m,theta_f,kappa,sigma,tau){
  y <- dt$stl
  t <- dt$time
  sex <- model.matrix(~ 0 + gender,dt)
  mm <- (sex %*% c(theta_f,theta_m))
  vars <- diag(sigma^2,nrow=nrow(dt))
  time_diff <- abs(matrix(t,nrow=nrow(dt),ncol=nrow(dt),byrow=TRUE) -
                   matrix(t,nrow=nrow(dt),ncol=nrow(dt)))
  covs <- tau^2/(1 - kappa^2)*kappa^time_diff
  vv <- vars + covs
  return(-dmvnorm(y,mean=mm,sigma=vv,log=TRUE) %>% sum)
}

minus_ll <- function(theta_m,theta_f,kappa_m,kappa_f,sigma_m,sigma_f,tau_m,tau_f){

}

final <- NULL

for(m in 1:M){
  print(paste0("Running Imputation ",m))
  dt <- dat %>% filter(imputation_number == m)
  test <- mle2(minusll,
               start=list(theta_m=56,theta_f=53,
                          kappa=0.5,
                          sigma=9,
                          tau=5))
  savs <- c(coef(test)[1:2],vcov(test)[1:2,1:2] %>% as.vector)
  final <- rbind(final,savs)
}

thetahat <- apply(final[,1:2],2,mean)
Bm <- apply(final[,3:6],2,mean) %>% matrix(nrow=2)
Vm <- cov(final[,1:2])
Tm <- (1 + 1/M)*Bm + Vm
tstat <- (c(1,-1) %*% thetahat) / sqrt(c(1,-1) %*% Tm %*% c(1,-1))
p <- 2*pnorm(tstat,lower.tail=FALSE)

