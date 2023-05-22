rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)

dat <- read.csv("Data/stuart_stickle_100imputations_all.csv")
dat$group <- factor(paste(dat$gender,letters[dat$time]))

mod <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  matrix[N,2*T] X;
  vector[N] stl;
  vector[N] lps;
  vector[N] ect;
  vector[N] tpg;
  vector[N] cle;
  vector[N] pmx;
  vector[N] ds1;
  vector[N] ds2;
  vector[N] ds3;
  vector[N] lpt;
  int<lower=0> mds[N];
  int<lower=0> mdf[N];
  int<lower=0> mav[N];
  int<lower=0> maf[N];
  int<lower=0> mcf[N];
  int<lower=0> mpt[N];
}
parameters{
  real<lower=0> theta_f[16];
  real<lower=0> theta_m[16];
  real<lower=-1,upper=1> kappa[16];
  real<lower=0> sigma[10];
  real<lower=0> tau[16];
  real beta[9];
  vector[2*T] u_stl;
  vector[2*T] u_lps;
  vector[2*T] u_ect;
  vector[2*T] u_tpg;
  vector[2*T] u_cle;
  vector[2*T] u_pmx;
  vector[2*T] u_ds1;
  vector[2*T] u_ds2;
  vector[2*T] u_ds3;
  vector[2*T] u_lpt;
  vector[2*T] u_mds;
  vector[2*T] u_mdf;
  vector[2*T] u_mav;
  vector[2*T] u_maf;
  vector[2*T] u_mcf;
  vector[2*T] u_mpt;
}
transformed parameters{
  vector[2*T] mu_stl;
  vector[2*T] mu_lps;
  vector[2*T] mu_ect;
  vector[2*T] mu_tpg;
  vector[2*T] mu_cle;
  vector[2*T] mu_pmx;
  vector[2*T] mu_ds1;
  vector[2*T] mu_ds2;
  vector[2*T] mu_ds3;
  vector[2*T] mu_lpt;
  vector[2*T] mu_mds;
  vector[2*T] mu_mdf;
  vector[2*T] mu_mav;
  vector[2*T] mu_maf;
  vector[2*T] mu_mcf;
  vector[2*T] mu_mpt;

  mu_stl[1:T] = theta_f[1] + u_stl[1:T];
  mu_stl[(T+1):(2*T)] = theta_m[1] + u_stl[(T+1):(2*T)];
  mu_lps[1:T] = theta_f[2] + u_lps[1:T];
  mu_lps[(T+1):(2*T)] = theta_m[2] + u_lps[(T+1):(2*T)];
  mu_ect[1:T] = theta_f[3] + u_ect[1:T];
  mu_ect[(T+1):(2*T)] = theta_m[3] + u_ect[(T+1):(2*T)];
  mu_tpg[1:T] = theta_f[4] + u_tpg[1:T];
  mu_tpg[(T+1):(2*T)] = theta_m[4] + u_tpg[(T+1):(2*T)];
  mu_cle[1:T] = theta_f[5] + u_cle[1:T];
  mu_cle[(T+1):(2*T)] = theta_m[5] + u_cle[(T+1):(2*T)];
  mu_pmx[1:T] = theta_f[6] + u_pmx[1:T];
  mu_pmx[(T+1):(2*T)] = theta_m[6] + u_pmx[(T+1):(2*T)];
  mu_ds1[1:T] = theta_f[7] + u_ds1[1:T];
  mu_ds1[(T+1):(2*T)] = theta_m[7] + u_ds1[(T+1):(2*T)];
  mu_ds2[1:T] = theta_f[8] + u_ds2[1:T];
  mu_ds2[(T+1):(2*T)] = theta_m[8] + u_ds2[(T+1):(2*T)];
  mu_ds3[1:T] = theta_f[9] + u_ds3[1:T];
  mu_ds3[(T+1):(2*T)] = theta_m[9] + u_ds3[(T+1):(2*T)];
  mu_lpt[1:T] = theta_f[10] + u_lpt[1:T];
  mu_lpt[(T+1):(2*T)] = theta_m[10] + u_lpt[(T+1):(2*T)];
  mu_mds[1:T] = theta_f[11] + u_mds[1:T];
  mu_mds[(T+1):(2*T)] = theta_m[11] + u_mds[(T+1):(2*T)];
  mu_mdf[1:T] = theta_f[12] + u_mdf[1:T];
  mu_mdf[(T+1):(2*T)] = theta_m[12] + u_mdf[(T+1):(2*T)];
  mu_mav[1:T] = theta_f[13] + u_mav[1:T];
  mu_mav[(T+1):(2*T)] = theta_m[13] + u_mav[(T+1):(2*T)];
  mu_maf[1:T] = theta_f[14] + u_maf[1:T];
  mu_maf[(T+1):(2*T)] = theta_m[14] + u_maf[(T+1):(2*T)];
  mu_mcf[1:T] = theta_f[15] + u_mcf[1:T];
  mu_mcf[(T+1):(2*T)] = theta_m[15] + u_mcf[(T+1):(2*T)];
  mu_mpt[1:T] = theta_f[16] + u_mpt[1:T];
  mu_mpt[(T+1):(2*T)] = theta_m[16] + u_mpt[(T+1):(2*T)];
}
model{
  stl ~ normal(X * mu_stl,sigma[1]);
  lps ~ normal(X * mu_lps + beta[1] * (stl - X * mu_stl),sigma[2]);
  ect ~ normal(X * mu_ect + beta[2] * (stl - X * mu_stl),sigma[3]);
  tpg ~ normal(X * mu_tpg + beta[3] * (stl - X * mu_stl),sigma[4]);
  cle ~ normal(X * mu_cle + beta[4] * (stl - X * mu_stl),sigma[5]);
  pmx ~ normal(X * mu_pmx + beta[5] * (stl - X * mu_stl),sigma[6]);
  ds1 ~ normal(X * mu_ds1 + beta[6] * (stl - X * mu_stl),sigma[7]);
  ds2 ~ normal(X * mu_ds2 + beta[7] * (stl - X * mu_stl),sigma[8]);
  ds3 ~ normal(X * mu_ds3 + beta[8] * (stl - X * mu_stl),sigma[9]);
  lpt ~ normal(X * mu_lpt + beta[9] * (stl - X * mu_stl),sigma[10]);
  mds ~ poisson(exp(X * mu_mds));
  mdf ~ poisson(exp(X * mu_mdf));
  mav ~ poisson(exp(X * mu_mav));
  maf ~ poisson(exp(X * mu_maf));
  mcf ~ poisson(exp(X * mu_mcf));
  mpt ~ poisson(exp(X * mu_mpt));

  u_stl[2:T] ~ normal(kappa[1]*u_stl[1:(T-1)],tau[1]);
  u_stl[(T+2):(2*T)] ~ normal(kappa[1]*u_stl[(T+1):(2*T-1)],tau[1]);
  u_lps[2:T] ~ normal(kappa[2]*u_lps[1:(T-1)],tau[2]);
  u_lps[(T+2):(2*T)] ~ normal(kappa[2]*u_lps[(T+1):(2*T-1)],tau[2]);
  u_ect[2:T] ~ normal(kappa[3]*u_ect[1:(T-1)],tau[3]);
  u_ect[(T+2):(2*T)] ~ normal(kappa[3]*u_ect[(T+1):(2*T-1)],tau[3]);
  u_tpg[2:T] ~ normal(kappa[4]*u_tpg[1:(T-1)],tau[4]);
  u_tpg[(T+2):(2*T)] ~ normal(kappa[4]*u_tpg[(T+1):(2*T-1)],tau[4]);
  u_cle[2:T] ~ normal(kappa[5]*u_cle[1:(T-1)],tau[5]);
  u_cle[(T+2):(2*T)] ~ normal(kappa[5]*u_cle[(T+1):(2*T-1)],tau[5]);
  u_pmx[2:T] ~ normal(kappa[6]*u_pmx[1:(T-1)],tau[6]);
  u_pmx[(T+2):(2*T)] ~ normal(kappa[6]*u_pmx[(T+1):(2*T-1)],tau[6]);
  u_ds1[2:T] ~ normal(kappa[7]*u_ds1[1:(T-1)],tau[7]);
  u_ds1[(T+2):(2*T)] ~ normal(kappa[7]*u_ds1[(T+1):(2*T-1)],tau[7]);
  u_ds2[2:T] ~ normal(kappa[8]*u_ds2[1:(T-1)],tau[8]);
  u_ds2[(T+2):(2*T)] ~ normal(kappa[8]*u_ds2[(T+1):(2*T-1)],tau[8]);
  u_ds3[2:T] ~ normal(kappa[9]*u_ds3[1:(T-1)],tau[9]);
  u_ds3[(T+2):(2*T)] ~ normal(kappa[9]*u_ds3[(T+1):(2*T-1)],tau[9]);
  u_lpt[2:T] ~ normal(kappa[10]*u_lpt[1:(T-1)],tau[10]);
  u_lpt[(T+2):(2*T)] ~ normal(kappa[10]*u_lpt[(T+1):(2*T-1)],tau[10]);
  u_mds[2:T] ~ normal(kappa[11]*u_mds[1:(T-1)],tau[11]);
  u_mds[(T+2):(2*T)] ~ normal(kappa[11]*u_mds[(T+1):(2*T-1)],tau[11]);
  u_mdf[2:T] ~ normal(kappa[12]*u_mdf[1:(T-1)],tau[12]);
  u_mdf[(T+2):(2*T)] ~ normal(kappa[12]*u_mdf[(T+1):(2*T-1)],tau[12]);
  u_mav[2:T] ~ normal(kappa[13]*u_mav[1:(T-1)],tau[13]);
  u_mav[(T+2):(2*T)] ~ normal(kappa[13]*u_mav[(T+1):(2*T-1)],tau[13]);
  u_maf[2:T] ~ normal(kappa[14]*u_maf[1:(T-1)],tau[14]);
  u_maf[(T+2):(2*T)] ~ normal(kappa[14]*u_maf[(T+1):(2*T-1)],tau[14]);
  u_mcf[2:T] ~ normal(kappa[15]*u_mcf[1:(T-1)],tau[15]);
  u_mcf[(T+2):(2*T)] ~ normal(kappa[15]*u_mcf[(T+1):(2*T-1)],tau[15]);
  u_mpt[2:T] ~ normal(kappa[16]*u_mpt[1:(T-1)],tau[16]);
  u_mpt[(T+2):(2*T)] ~ normal(kappa[16]*u_mpt[(T+1):(2*T-1)],tau[16]);

  u_stl[{1,T+1}] ~ normal(0,tau[1]/sqrt(1-kappa[1]*kappa[1]));
  u_lps[{1,T+1}] ~ normal(0,tau[2]/sqrt(1-kappa[2]*kappa[2]));
  u_ect[{1,T+1}] ~ normal(0,tau[3]/sqrt(1-kappa[3]*kappa[3]));
  u_tpg[{1,T+1}] ~ normal(0,tau[4]/sqrt(1-kappa[4]*kappa[4]));
  u_cle[{1,T+1}] ~ normal(0,tau[5]/sqrt(1-kappa[5]*kappa[5]));
  u_pmx[{1,T+1}] ~ normal(0,tau[6]/sqrt(1-kappa[6]*kappa[6]));
  u_ds1[{1,T+1}] ~ normal(0,tau[7]/sqrt(1-kappa[7]*kappa[7]));
  u_ds2[{1,T+1}] ~ normal(0,tau[8]/sqrt(1-kappa[8]*kappa[8]));
  u_ds3[{1,T+1}] ~ normal(0,tau[9]/sqrt(1-kappa[9]*kappa[9]));
  u_lpt[{1,T+1}] ~ normal(0,tau[10]/sqrt(1-kappa[10]*kappa[10]));
  u_mds[{1,T+1}] ~ normal(0,tau[11]/sqrt(1-kappa[11]*kappa[11]));
  u_mdf[{1,T+1}] ~ normal(0,tau[12]/sqrt(1-kappa[12]*kappa[12]));
  u_mav[{1,T+1}] ~ normal(0,tau[13]/sqrt(1-kappa[13]*kappa[13]));
  u_maf[{1,T+1}] ~ normal(0,tau[14]/sqrt(1-kappa[14]*kappa[14]));
  u_mcf[{1,T+1}] ~ normal(0,tau[15]/sqrt(1-kappa[15]*kappa[15]));
  u_mpt[{1,T+1}] ~ normal(0,tau[16]/sqrt(1-kappa[16]*kappa[16]));

  theta_f ~ normal(0,100);
  theta_m ~ normal(0,100);
  kappa ~ normal(0,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
  beta ~ normal(0,5);
}
"
sm <- stan_model(model_code = mod, verbose = TRUE)
M <- max(dat$imp)

for (m in 1:M){
  tmp <- dat %>% filter(imp == m)
  stl <- tmp$stl
  mds <- tmp$mds
  mdf <- tmp$mdf
  maf <- tmp$maf
  lps <- tmp$lps
  ect <- tmp$ect
  tpg <- tmp$tpg
  mav <- tmp$mav
  mcv <- tmp$mcv
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
                  list(N=N,T=T,X=X,stl=stl,
                       lpt=lpt,mds=mds,mdf=mdf,
                       maf=maf,lps=lps,ect=ect,
                       tpg=tpg,mav=mav,mcv=mcv,
                       mpt=mpt,cle=cle,pmx=pmx,
                       ds1=ds1,ds2=ds2,ds3=ds3),
                  iter=6000,chains=1,warmup=5000,thin=10,
                  pars = c("theta_f","theta_m","kappa","sigma","tau","beta",
                           "mu_stl","mu_lps","mu_ect","mu_tpg",
                           "mu_cle","mu_pmx","mu_ds1","mu_ds2",
                           "mu_ds3","mu_lpt","mu_mds","mu_mdf",
                           "mu_mav","mu_maf","mu_mcf","mu_mpt"))
  if (m == 1){
    samps <- as.data.frame(fit)
  } else {
    samps <- rbind(samps,as.data.frame(fit))
  }
}

saveRDS(samps,"posterior_samples.RDS")
