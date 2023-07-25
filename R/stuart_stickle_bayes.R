rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)

dat <- read.csv("Data/updated_stickle_imputations_100.csv")
times <- read.csv("RawData_ReadOnly/KSampleMeanTimes.csv")
dat$group <- factor(paste(dat$gender,letters[dat$time]))

cont <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  vector[N] y;
  vector[N] stl;
  matrix[N,2*T] X;
  vector[T] time;
  vector[T-1] Delta;
  int<lower=0,upper=3> mod;
  int<lower=0,upper=3> mod_stl;
}
parameters{
  real beta0_f;
  real beta0_m;
  real beta0_f_stl;
  real beta0_m_stl;
  real beta1_f;
  real beta1_m;
  real beta1_f_stl;
  real beta1_m_stl;
  real<lower=0> kappa;
  real<lower=0> kappa_stl;
  real<lower=0> sigma;
  real<lower=0> sigma_stl;
  real<lower=0> tau;
  real<lower=0> tau_stl;
  real<lower=0> tau_0;
  real<lower=0> tau_0_stl;
  real gamma;
  vector[2*T] u;
  vector[2*T] u_stl;
}
transformed parameters{
  vector[2*T] mu;
  vector[2*T] mu_stl;

  if (mod == 0){
    mu[1:T] = beta0_f + beta1_f*time + u[1:T];
    mu[(T+1):(2*T)] = beta0_m + beta1_m*time + u[(T+1):(2*T)];
  } else if (mod == 1) {
    mu[1:T] = beta0_f + beta1_f*time;
    mu[(T+1):(2*T)] = beta0_m + beta1_m*time;
  } else if (mod == 2){
    mu[1:T] = beta0_f + u[1:T];
    mu[(T+1):(2*T)] = beta0_m + u[(T+1):(2*T)];
  } else {
    mu[1:T] = beta0_f + 0*time;
    mu[(T+1):(2*T)] = beta0_m + 0*time;
  }

  if (mod_stl==0){
    mu_stl[1:T] = beta0_f + beta1_f*time + u_stl[1:T];
    mu_stl[(T+1):(2*T)] = beta0_m + beta1_m*time + u_stl[(T+1):(2*T)];
  } else if (mod_stl == 1) {
    mu_stl[1:T] = beta0_f + beta1_f*time;
    mu_stl[(T+1):(2*T)] = beta0_m + beta1_m*time;
  } else if (mod_stl == 2){
    mu_stl[1:T] = beta0_f + u_stl[1:T];
    mu_stl[(T+1):(2*T)] = beta0_m + u_stl[(T+1):(2*T)];
  } else {
    mu_stl[1:T] = beta0_f + 0*time;
    mu_stl[(T+1):(2*T)] = beta0_m + 0*time;
  }

  vector[T-1] OU_mean;
  vector[T-1] OU_var;
  vector[T-1] OU_mean_stl;
  vector[T-1] OU_var_stl;

  OU_mean = exp(-kappa*Delta);
  OU_var = pow(tau,2)*(1 - exp(-2*kappa*Delta))/(2*kappa);
  OU_mean_stl = exp(-kappa_stl*Delta);
  OU_var_stl = pow(tau_stl,2)*(1 - exp(-2*kappa_stl*Delta))/(2*kappa_stl);
}

model{
  // Data Model
  y ~ normal(X * mu + gamma * (stl - X * mu_stl),sigma);
  stl ~ normal(X * mu_stl,sigma_stl);

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(OU_mean .* u[1:(T-1)],sqrt(OU_var));
  u[(T+2):(2*T)] ~ normal(OU_mean .* u[(T+1):(2*T-1)],sqrt(OU_var));

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau_0);

  // OU Process for Gender and Time Specific Means STL
  u_stl[2:T] ~ normal(OU_mean_stl .* u_stl[1:(T-1)],sqrt(OU_var_stl));
  u_stl[(T+2):(2*T)] ~ normal(OU_mean_stl .* u_stl[(T+1):(2*T-1)],sqrt(OU_var_stl));

  // Initial Gender and Time Mean Model STL
  u_stl[{1,T+1}] ~ normal(0,tau_0_stl);

  // Priors
  beta0_f ~ normal(0,100);
  beta0_m ~ normal(0,100);
  beta1_f ~ normal(0,3);
  beta1_m ~ normal(0,3);
  kappa ~ normal(0,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
  tau_0 ~ normal(0,20);
  beta0_f_stl ~ normal(0,100);
  beta0_m_stl ~ normal(0,100);
  beta1_f_stl ~ normal(0,3);
  beta1_m_stl ~ normal(0,3);
  kappa_stl ~ normal(0,1);
  sigma_stl ~ normal(0,10);
  tau_stl ~ normal(0,10);
  tau_0_stl ~ normal(0,20);
  gamma ~ normal(0,5);
}
"

disc <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  int y[N];
  vector[N] stl;
  matrix[N,2*T] X;
  vector[T] time;
  vector[T-1] Delta;
  int<lower=0,upper=3> mod;
  int<lower=0,upper=3> mod_stl;
  int<lower=0,upper=1> ind;
}
parameters{
  real beta0_f;
  real beta0_m;
  real beta0_f_stl;
  real beta0_m_stl;
  real beta1_f;
  real beta1_m;
  real beta1_f_stl;
  real beta1_m_stl;
  real<lower=0> kappa;
  real<lower=0> kappa_stl;
  real<lower=0> sigma_stl;
  real<lower=0> tau;
  real<lower=0> tau_stl;
  real<lower=0> tau_0;
  real<lower=0> tau_0_stl;
  real gamma;
  vector[2*T] u;
  vector[2*T] u_stl;
}
transformed parameters{
  vector[2*T] mu;
  vector[2*T] mu_stl;

  if (mod == 0){
    mu[1:T] = beta0_f + beta1_f*time + u[1:T];
    mu[(T+1):(2*T)] = beta0_m + beta1_m*time + u[(T+1):(2*T)];
  } else if (mod == 1) {
    mu[1:T] = beta0_f + beta1_f*time;
    mu[(T+1):(2*T)] = beta0_m + beta1_m*time;
  } else if (mod == 2){
    mu[1:T] = beta0_f + u[1:T];
    mu[(T+1):(2*T)] = beta0_m + u[(T+1):(2*T)];
  } else {
    mu[1:T] = beta0_f + 0*time;
    mu[(T+1):(2*T)] = beta0_m + 0*time;
  }

  if (mod_stl==0){
    mu_stl[1:T] = beta0_f + beta1_f*time + u_stl[1:T];
    mu_stl[(T+1):(2*T)] = beta0_m + beta1_m*time + u_stl[(T+1):(2*T)];
  } else if (mod_stl == 1) {
    mu_stl[1:T] = beta0_f + beta1_f*time;
    mu_stl[(T+1):(2*T)] = beta0_m + beta1_m*time;
  } else if (mod_stl == 2){
    mu_stl[1:T] = beta0_f + u_stl[1:T];
    mu_stl[(T+1):(2*T)] = beta0_m + u_stl[(T+1):(2*T)];
  } else {
    mu_stl[1:T] = beta0_f + 0*time;
    mu_stl[(T+1):(2*T)] = beta0_m + 0*time;
  }

  vector[T-1] OU_mean;
  vector[T-1] OU_var;
  vector[T-1] OU_mean_stl;
  vector[T-1] OU_var_stl;

  OU_mean = exp(-kappa*Delta);
  OU_var = pow(tau,2)*(1 - exp(-2*kappa*Delta))/(2*kappa);
  OU_mean_stl = exp(-kappa_stl*Delta);
  OU_var_stl = pow(tau_stl,2)*(1 - exp(-2*kappa_stl*Delta))/(2*kappa_stl);
}

model{
  // Data Model
  if (ind == 1){
    y ~ poisson_log(X * mu);
  } else {
    y ~ poisson_log(X * mu + gamma * (stl - X * mu_stl) - pow(gamma*sigma_stl,2)/2);
  }
  stl ~ normal(X * mu_stl,sigma_stl);

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(OU_mean .* u[1:(T-1)],sqrt(OU_var));
  u[(T+2):(2*T)] ~ normal(OU_mean .* u[(T+1):(2*T-1)],sqrt(OU_var));

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau_0);

  // OU Process for Gender and Time Specific Means STL
  u_stl[2:T] ~ normal(OU_mean_stl .* u_stl[1:(T-1)],sqrt(OU_var_stl));
  u_stl[(T+2):(2*T)] ~ normal(OU_mean_stl .* u_stl[(T+1):(2*T-1)],sqrt(OU_var_stl));

  // Initial Gender and Time Mean Model STL
  u_stl[{1,T+1}] ~ normal(0,tau_0_stl);

  // Priors
  beta0_f ~ normal(0,100);
  beta0_m ~ normal(0,100);
  beta1_f ~ normal(0,3);
  beta1_m ~ normal(0,3);
  kappa ~ normal(0,1);
  tau ~ normal(0,10);
  tau_0 ~ normal(0,20);
  beta0_f_stl ~ normal(0,100);
  beta0_m_stl ~ normal(0,100);
  beta1_f_stl ~ normal(0,3);
  beta1_m_stl ~ normal(0,3);
  kappa_stl ~ normal(0,1);
  sigma_stl ~ normal(0,10);
  tau_stl ~ normal(0,10);
  tau_0_stl ~ normal(0,20);
  gamma ~ normal(0,5);
}
"
stan_cont <- stan_model(model_code = cont)
stan_disc <- stan_model(model_code = disc)
M <- max(dat$imp)

vars <- c("lps","ect","tpg",
          "cle","pmx","ds1","ds2",
          "ds3","lpt","mdf","mav",
          "maf","mcv","mds","mpt")

gender <- matrix((dat$gender=="Female")*1,ncol=M)
stl <- matrix(dat$stl,ncol=M)
index <- dat$time[1:nrow(gender)]

N <- nrow(gender)
T <- max(dat$time)

time <- times$Inverted.Year/1000
Delta <- diff(times$Inverted.Year)/1000

for (mod in 0:3){
  for (mod_stl in 0:3){
    for (var in vars){
      samps <- NULL
      y <- matrix(dat[,var],ncol=M)

      for (m in 1:M){
        print(paste0("Running Model ",mod,"-",mod_stl," for variable ",var," for Iteration ",m))
        X <- model.matrix(~0+group,dat %>% filter(imp==m))
        if(var == "mav" | var == "mcv"){
          fit <- sampling(stan_disc,
                          list(N=N,T=T,X=X,
                               y=y[,m],stl=stl[,m],
                               index=index,time=time,Delta=Delta,
                               ind=0,mod=mod,mod_stl=mod_stl),
                          iter=1100,chains=1,warmup=1000,
                          pars = c("beta0_f","beta0_m",
                                   "beta1_f","beta1_m",
                                   "beta0_f_stl","beta0_m_stl",
                                   "beta1_f_stl","beta1_m_stl",
                                   "kappa","kappa_stl","sigma_stl",
                                   "tau","tau_stl","tau_0","tau_0_stl",
                                   "gamma","mu","mu_stl"))
          samps <- rbind(samps,as.data.frame(fit))
        } else if (substr(var,1,1) == "m") {
          if (mod_stl == 0){
            fit <- sampling(stan_disc,
                            list(N=N,T=T,X=X,
                                 y=y[,m],stl=stl[,m],
                                 index=index,time=time,Delta=Delta,
                                 ind=1,mod=mod,mod_stl=mod_stl),
                            iter=1100,chains=1,warmup=1000,
                            pars = c("beta0_f","beta0_m",
                                     "beta1_f","beta1_m",
                                     "beta0_f_stl","beta0_m_stl",
                                     "beta1_f_stl","beta1_m_stl",
                                     "kappa","kappa_stl","sigma_stl",
                                     "tau","tau_stl","tau_0","tau_0_stl",
                                     "gamma","mu","mu_stl"))
            samps <- rbind(samps,as.data.frame(fit))
          }
        } else {
          fit <- sampling(stan_cont,
                          list(N=N,T=T,X=X,
                               y=y[,m],stl=stl[,m],
                               index=index,time=time,Delta=Delta,
                               mod=mod,mod_stl=mod_stl),
                          iter=1100,chains=1,warmup=1000,
                          pars = c("beta0_f","beta0_m",
                                   "beta1_f","beta1_m",
                                   "beta0_f_stl","beta0_m_stl",
                                   "beta1_f_stl","beta1_m_stl",
                                   "kappa","kappa_stl","sigma","sigma_stl",
                                   "tau","tau_stl","tau_0","tau_0_stl",
                                   "gamma","mu","mu_stl"))
          samps <- rbind(samps,as.data.frame(fit))
        }
      }
      samps_name <- paste0("Figures/",var,"/posterior_samples_",var,"_",mod,mod_stl,".RDS")
      saveRDS(samps,file=samps_name)
    }
  }
}
