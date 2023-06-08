rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)

dat <- read.csv("Data/stuart_stickle_100imputations_all.csv")
dat$group <- factor(paste(dat$gender,letters[dat$time]))

cont_ind <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  matrix[N,2*T] X;
  vector[N] y;
}
parameters{
  real<lower=0> theta_f;
  real<lower=0> theta_m;
  real<lower=-1,upper=1> kappa;
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[2*T] u;
}
transformed parameters{
  vector[2*T] mu;
  mu[1:T] = theta_f + u[1:T];
  mu[(T+1):(2*T)] = theta_m + u[(T+1):(2*T)];
}

model{
  // Data Model
  y ~ normal(X * mu,sigma);

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(kappa*u[1:(T-1)],tau);
  u[(T+2):(2*T)] ~ normal(kappa*u[(T+1):(2*T-1)],tau);

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau/sqrt(1-kappa*kappa));

  // Priors
  theta_f ~ normal(0,100);
  theta_m ~ normal(0,100);
  kappa ~ normal(0,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
}
"

cont_dep <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  matrix[N,2*T] X;
  vector[N] y;
  vector[N] stl;
}
parameters{
  real<lower=0> theta_f;
  real<lower=0> theta_m;
  real<lower=-1,upper=1> kappa;
  real<lower=0> sigma;
  real<lower=0> tau;
  real beta;
  vector[2*T] u;
}
transformed parameters{
  vector[2*T] mu;
  mu[1:T] = theta_f + u[1:T];
  mu[(T+1):(2*T)] = theta_m + u[(T+1):(2*T)];
}

model{
  // Data Model
  y ~ normal(X * mu + beta * stl,sigma);

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(kappa*u[1:(T-1)],tau);
  u[(T+2):(2*T)] ~ normal(kappa*u[(T+1):(2*T-1)],tau);

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau/sqrt(1-kappa*kappa));

  // Priors
  theta_f ~ normal(0,100);
  theta_m ~ normal(0,100);
  kappa ~ normal(0,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
  beta ~ normal(0,5);
}
"

disc_ind <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  matrix[N,2*T] X;
  int<lower=0>[N] y;
}
parameters{
  real<lower=0> theta_f;
  real<lower=0> theta_m;
  real<lower=-1,upper=1> kappa;
  real<lower=0> tau;
  vector[2*T] u;
}
transformed parameters{
  vector[2*T] mu;
  mu[1:T] = theta_f + u[1:T];
  mu[(T+1):(2*T)] = theta_m + u[(T+1):(2*T)];
}

model{
  // Data Model
  y ~ poisson_log(X * mu);

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(kappa*u[1:(T-1)],tau);
  u[(T+2):(2*T)] ~ normal(kappa*u[(T+1):(2*T-1)],tau);

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau/sqrt(1-kappa*kappa));

  // Priors
  theta_f ~ normal(0,100);
  theta_m ~ normal(0,100);
  kappa ~ normal(0,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
}
"

disc_dep <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  matrix[N,2*T] X;
  int<lower=0>[N] y;
  vector[N] stl;
}
parameters{
  real<lower=0> theta_f;
  real<lower=0> theta_m;
  real<lower=-1,upper=1> kappa;
  real<lower=0> tau;
  real beta;
  vector[2*T] u;
}
transformed parameters{
  vector[2*T] mu;
  mu[1:T] = theta_f + u[1:T];
  mu[(T+1):(2*T)] = theta_m + u[(T+1):(2*T)];
}

model{
  // Data Model
  y ~ poisson_log(X * mu + beta * stl);

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(kappa*u[1:(T-1)],tau);
  u[(T+2):(2*T)] ~ normal(kappa*u[(T+1):(2*T-1)],tau);

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau/sqrt(1-kappa*kappa));

  // Priors
  theta_f ~ normal(0,100);
  theta_m ~ normal(0,100);
  kappa ~ normal(0,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
  beta ~ normal(0,5);
}
"
stan_cont_ind <- stan_model(model_code = cont_ind)
stan_cont_dep <- stan_model(model_code = cont_dep)
stan_disc_ind <- stan_model(model_code = disc_ind)
stan_disc_dep <- stan_model(model_code = disc_dep)
M <- max(dat$imp)

vars <- c("stl","lps","ect","tpg",
          "cle","pmx","ds1","ds2",
          "ds3","lpt","mdf","mav",
          "maf","mcv","mds","mpt")

for (m in 1:M){
  tmp <- dat %>% filter(imp == m)
  stl <- tmp$stl
  mdf <- tmp$mdf
  maf <- tmp$maf
  lps <- tmp$lps
  ect <- tmp$ect
  tpg <- tmp$tpg
  mav <- tmp$mav
  mcv <- tmp$mcv
  mds <- tmp$mds
  mpt <- tmp$mpt
  cle <- tmp$cle
  pmx <- tmp$pmx
  ds1 <- tmp$ds1
  ds2 <- tmp$ds2
  ds3 <- tmp$ds3
  lpt <- tmp$lpt
  X <- model.matrix(~0+group,tmp)

  N <- length(stl)
  T <- max(tmp$time)

  print(paste("Running Imputation",m))
  fit <- sampling(sm,
                  list(N=N,T=T,X=X,mds=mds,
                       stl=stl,lpt=lpt,mdf=mdf,
                       maf=maf,lps=lps,ect=ect,
                       tpg=tpg,mav=mav,mcv=mcv,
                       mpt=mpt,cle=cle,pmx=pmx,
                       ds1=ds1,ds2=ds2,ds3=ds3),
                  iter=1500,chains=1,warmup=1000,thin=5,
                  pars = c("theta_f","theta_m","kappa","sigma","tau","beta",
                           "mu_stl","mu_lps","mu_ect","mu_tpg",
                           "mu_cle","mu_pmx","mu_ds1","mu_ds2",
                           "mu_ds3","mu_lpt","mu_mdf","mu_mav",
                           "mu_maf","mu_mcv","mu_mds","mu_mpt"))
  if (m == 1){
    samps <- as.data.frame(fit)
  } else {
    samps <- rbind(samps,as.data.frame(fit))
  }
}

# saveRDS(samps,"posterior_samples.RDS")
saveRDS(samps,"posterior_samples_mav_adj.RDS")
