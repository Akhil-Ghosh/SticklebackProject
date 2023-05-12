rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)

dat <- read.csv("Data/stuart_stickle_100imputations_all.csv")

mod <- "
data{
  int<lower = 0> N;
  int<lower = 0> T;
  real<lower = 0> y[N];
  int<lower = 0,upper = 1> gender[N];
  int<lower = 1,upper = T> time[N];
}
parameters{
  real<lower=0> theta_f;
  real<lower=0> theta_m;
  real u_f[T];
  real u_m[T];
  real<lower=-1,upper=1> kappa;
  real<lower=0> sigma;
  real<lower=0> tau;
}
transformed parameters{
  real mu_f[T];
  real mu_m[T];
  for (t in 1:T){
    mu_f[t] = theta_f + u_f[t];
    mu_m[t] = theta_m + u_m[t];
  }
}
model{
  for (i in 1:N){
    y[i] ~ normal(gender[i]*(mu_f[time[i]]) +
                  (1 - gender[i])*(mu_m[time[i]]),sigma);
  }
  for (t in 2:T){
    u_f[t] ~ normal(kappa*u_f[t-1],tau);
    u_m[t] ~ normal(kappa*u_m[t-1],tau);
  }
  u_f[1] ~ normal(0,tau/sqrt(1-kappa*kappa));
  u_m[1] ~ normal(0,tau/sqrt(1-kappa*kappa));

  theta_f ~ normal(54,20);
  theta_m ~ normal(54,20);
  kappa ~ normal(0.5,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
}
"
sm <- stan_model(model_code = mod, verbose = TRUE)
M <- max(dat$imp)

for (m in 1:M){
  tmp <- dat %>% filter(imp == m)
  stl <- tmp$stl
  cle <- tmp$cle.sc
  pmx <- tmp$pmx.sc
  gender <- ifelse(tmp$gender == "Female",1,0)
  time <- tmp$time

  N <- length(stl)
  T <- max(time)

  print(paste("Running Imputation",m,"for stl"))
  fit.stl <- sampling(sm,
                      list(N=N,T=T,y=stl,gender=gender,time=time),
                      iter=10000,chains=1,warmup=5000,thin=5,
                      pars = c("theta_f","theta_m","mu_f","mu_m","kappa","sigma","tau"))
  print(paste("Running Imputation",m,"for cle"))
  fit.cle <- sampling(sm,
                      list(N=N,T=T,y=cle,gender=gender,time=time),
                      iter=10000,chains=1,warmup=5000,thin=5,
                      pars = c("theta_f","theta_m","mu_f","mu_m","kappa","sigma","tau"))
  print(paste("Running Imputation",m,"for pmx"))
  fit.pmx <- sampling(sm,
                      list(N=N,T=T,y=pmx,gender=gender,time=time),
                      iter=10000,chains=1,warmup=5000,thin=5,
                      pars = c("theta_f","theta_m","mu_f","mu_m","kappa","sigma","tau"))

  if (m == 1){
    samps.stl <- as.data.frame(fit.stl)
    samps.cle <- as.data.frame(fit.cle)
    samps.pmx <- as.data.frame(fit.pmx)
  } else {
    samps.stl <- rbind(samps.stl,as.data.frame(fit.stl))
    samps.cle <- rbind(samps.cle,as.data.frame(fit.cle))
    samps.pmx <- rbind(samps.stl,as.data.frame(fit.pmx))
  }
}

saveRDS(samps.stl,"posterior_samples_stl.RDS")
saveRDS(samps.cle,"posterior_samples_cle.RDS")
saveRDS(samps.pmx,"posterior_samples_pmx.RDS")

