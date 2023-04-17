library(tidyverse)
library(tmvtnorm)
library(bbmle)

dat <- read.csv("stuart_stickle_imp_raw.csv")
M <- max(dat$imputation_number)

minusll <- function(theta_m,theta_f,kappa,sigma,tau,rho){ #m0_m,m0_f,theta_m,theta_f,kappa,sigma,tau,rho){
  y <- dt$stl
  t <- dt$time
  sex <- model.matrix(~ 0 + gender,dt)
  # mm <- (sex %*% c(m0_f,m0_m))*exp(-kappa*t) +
  #       (sex %*% c(theta_f,theta_m))*(1 - exp(-kappa*t))
  mm <- sex %*% c(theta_f,theta_m)
  vars <- diag(sigma^2,nrow=nrow(dt))
  rho_mat <- matrix(c(1,rho,rho,1),nrow=2)
  time_diff <- abs(matrix(t,nrow=nrow(dt),ncol=nrow(dt),byrow=TRUE) -
                   matrix(t,nrow=nrow(dt),ncol=nrow(dt)))
  # time_abs  <- abs(matrix(t,nrow=nrow(dt),ncol=nrow(dt),byrow=TRUE) +
  #                  matrix(t,nrow=nrow(dt),ncol=nrow(dt)))
  # covs <- tau^2/(2*kappa)*(sex %*% rho_mat %*% t(sex))*
  #  (exp(-kappa*time_diff) - exp(-kappa*time_abs))
  covs <- tau^2/(1 - kappa^2)*(sex %*% rho_mat %*% t(sex))*kappa^time_diff
  vv <- vars + covs
  return(-dmvnorm(y,mean=mm,sigma=vv,log=TRUE) %>% sum)
}

final <- NULL

for(m in 1:M){
  print(paste0("Running Imputation ",m))
  dt <- dat %>% filter(imputation_number == m)
  test <- mle2(minusll,
               start=list(# m0_m=70,m0_f=69,
                          theta_m=53,theta_f=50,
                          kappa=0.4,
                          sigma=9,
                          tau=5,rho=0.8))
  savs <- c(coef(test)[1:2],vcov(test)[1:2,1:2] %>% as.vector)
  # savs <- c(coef(test)[3:4],vcov(test)[3:4,3:4] %>% as.vector)
  final <- rbind(final,savs)
}

thetahat <- apply(final[,1:2],2,mean)
Bm <- apply(final[,3:6],2,mean) %>% matrix(nrow=2)
Vm <- cov(final[,1:2])
Tm <- (1 + 1/M)*Bm + Vm
tstat <- (c(1,-1) %*% thetahat) / sqrt(c(1,-1) %*% Tm %*% c(1,-1))
p <- 2*pnorm(tstat,lower.tail=FALSE)

