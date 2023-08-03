rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)

dat <- read.csv("Data/updated_stickle_imputations_100.csv")
times <- read.csv("RawData_ReadOnly/KSampleMeanTimes.csv")
times$time <- c(1:18)
dat <- dat %>% left_join(times %>% select(time,Inverted.Year))
dat$Inverted.Year <- dat$Inverted.Year/1000
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
  int<lower=0,upper=1> ind;
}
parameters{
  real beta0_f;
  real beta0_m;
  real beta1_f;
  real beta1_m;
  real<lower=0> kappa;
  real<lower=0> sigma;
  real<lower=0> tau;
  real<lower=0> tau_0;
  real gamma;
  vector[2*T] u;
}
transformed parameters{
  vector[2*T] mu;
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

  vector[T-1] OU_mean;
  vector[T-1] OU_var;

  OU_mean = exp(-kappa*Delta);
  OU_var = pow(tau,2)*(1 - exp(-2*kappa*Delta))/(2*kappa);
}

model{
  // Data Model
  if (ind == 1){
    y ~ normal(X * mu,sigma);
  } else {
    y ~ normal(X * mu + gamma * stl,sigma);
  }

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(OU_mean .* u[1:(T-1)],sqrt(OU_var));
  u[(T+2):(2*T)] ~ normal(OU_mean .* u[(T+1):(2*T-1)],sqrt(OU_var));

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau_0);

  // Priors
  beta0_f ~ normal(0,100);
  beta0_m ~ normal(0,100);
  beta1_f ~ normal(0,3);
  beta1_m ~ normal(0,3);
  kappa ~ normal(0,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
  tau_0 ~ normal(0,20);
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
  int<lower=0,upper=1> ind;
}
parameters{
  real beta0_f;
  real beta0_m;
  real beta1_f;
  real beta1_m;
  real beta1_f_stl;
  real beta1_m_stl;
  real<lower=0> kappa;
  real<lower=0> tau;
  real<lower=0> tau_0;
  real gamma;
  vector[2*T] u;
}
transformed parameters{
  vector[2*T] mu;

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

  vector[T-1] OU_mean;
  vector[T-1] OU_var;

  OU_mean = exp(-kappa*Delta);
  OU_var = pow(tau,2)*(1 - exp(-2*kappa*Delta))/(2*kappa);
}

model{
  // Data Model
  if (ind == 1){
    y ~ poisson_log(X * mu);
  } else {
    y ~ poisson_log(X * mu + gamma * stl);
  }

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(OU_mean .* u[1:(T-1)],sqrt(OU_var));
  u[(T+2):(2*T)] ~ normal(OU_mean .* u[(T+1):(2*T-1)],sqrt(OU_var));

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau_0);

  // Priors
  beta0_f ~ normal(0,100);
  beta0_m ~ normal(0,100);
  beta1_f ~ normal(0,3);
  beta1_m ~ normal(0,3);
  kappa ~ normal(0,1);
  tau ~ normal(0,10);
  tau_0 ~ normal(0,20);
  gamma ~ normal(0,5);
}
"
stan_cont <- stan_model(model_code = cont)
stan_disc <- stan_model(model_code = disc)
M <- max(dat$imp)

vars <- c("stl","lps","ect","tpg",
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
  for (var in vars){
    samps <- NULL
    y <- matrix(dat[,var],ncol=M)

    for (m in 1:M){
      print(paste0("Running Model ",mod," for variable ",var," for Iteration ",m))
      X <- model.matrix(~0+group,dat %>% filter(imp==m))
      if(var == "mav" | var == "mcv"){
        fit <- sampling(stan_disc,
                        list(N=N,T=T,X=X,
                             y=y[,m],stl=stl[,m],
                             time=time,Delta=Delta,
                             ind=0,mod=mod),
                        iter=1500,chains=1,warmup=1000,thin=5,
                        pars = c("beta0_f","beta0_m",
                                 "beta1_f","beta1_m",
                                 "kappa",
                                 "tau","tau_0",
                                 "gamma","mu"))
        samps <- rbind(samps,as.data.frame(fit))
      } else if (substr(var,1,1) == "m") {
        fit <- sampling(stan_disc,
                        list(N=N,T=T,X=X,
                             y=y[,m],stl=stl[,m],
                             time=time,Delta=Delta,
                             ind=1,mod=mod),
                        iter=1500,chains=1,warmup=1000,thin=5,
                        pars = c("beta0_f","beta0_m",
                                 "beta1_f","beta1_m",
                                 "kappa",
                                 "tau","tau_0",
                                 "gamma","mu"))
        samps <- rbind(samps,as.data.frame(fit))
      } else if (var=="stl") {
        fit <- sampling(stan_cont,
                        list(N=N,T=T,X=X,
                             y=y[,m],stl=stl[,m],
                             time=time,Delta=Delta,
                             ind=1,mod=mod),
                        iter=1500,chains=1,warmup=1000,thin=5,
                        pars = c("beta0_f","beta0_m",
                                 "beta1_f","beta1_m",
                                 "kappa","sigma",
                                 "tau","tau_0",
                                 "gamma","mu"))
        samps <- rbind(samps,as.data.frame(fit))
      } else {
        fit <- sampling(stan_cont,
                        list(N=N,T=T,X=X,
                             y=y[,m],stl=stl[,m],
                             time=time,Delta=Delta,
                             ind=0,mod=mod),
                        iter=1500,chains=1,warmup=1000,thin=5,
                        pars = c("beta0_f","beta0_m",
                                 "beta1_f","beta1_m",
                                 "kappa","sigma",
                                 "tau","tau_0",
                                 "gamma","mu"))
        samps <- rbind(samps,as.data.frame(fit))
      }
    }
    samps_name <- paste0("Figures/",var,"/posterior_samples_",var,"_",mod,".RDS")
    saveRDS(samps,file=samps_name)
  }
}
