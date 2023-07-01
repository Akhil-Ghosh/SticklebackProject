rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)

dat <- read.csv("Data/updated_stickle_imputations_100.csv")
dat$group <- factor(paste(dat$gender,letters[dat$time]))

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
  real<lower=0> theta_f_stl;
  real<lower=0> theta_m_stl;
  real<lower=-1,upper=1> kappa;
  real<lower=-1,upper=1> kappa_stl;
  real<lower=0> sigma;
  real<lower=0> sigma_stl;
  real<lower=0> tau;
  real<lower=0> tau_stl;
  real beta;
  vector[2*T] u;
  vector[2*T] u_stl;
}
transformed parameters{
  vector[2*T] mu;
  mu[1:T] = theta_f + u[1:T];
  mu[(T+1):(2*T)] = theta_m + u[(T+1):(2*T)];
  vector[2*T] mu_stl;
  mu_stl[1:T] = theta_f_stl + u_stl[1:T];
  mu_stl[(T+1):(2*T)] = theta_m_stl + u_stl[(T+1):(2*T)];
}

model{
  // Data Model
  y ~ normal(X * mu + beta * (stl - X * mu_stl),sigma);
  stl ~ normal(X * mu_stl,sigma_stl);

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(kappa*u[1:(T-1)],tau);
  u[(T+2):(2*T)] ~ normal(kappa*u[(T+1):(2*T-1)],tau);
  u_stl[2:T] ~ normal(kappa_stl*u_stl[1:(T-1)],tau_stl);
  u_stl[(T+2):(2*T)] ~ normal(kappa_stl*u_stl[(T+1):(2*T-1)],tau_stl);

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau/sqrt(1-kappa*kappa));
  u_stl[{1,T+1}] ~ normal(0,tau_stl/sqrt(1-kappa_stl*kappa_stl));

  // Priors
  theta_f ~ normal(0,100);
  theta_m ~ normal(0,100);
  kappa ~ normal(0,1);
  sigma ~ normal(0,10);
  tau ~ normal(0,10);
  theta_f_stl ~ normal(0,100);
  theta_m_stl ~ normal(0,100);
  kappa_stl ~ normal(0,1);
  sigma_stl ~ normal(0,10);
  tau_stl ~ normal(0,10);
  beta ~ normal(0,5);
}
"

disc_ind <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  matrix[N,2*T] X;
  int<lower=0> y[N];
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
  vector[2*T] lambda = exp(mu);
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
  tau ~ normal(0,10);
}
"

disc_dep <- "
data{
  int<lower=0> N;
  int<lower=0> T;
  matrix[N,2*T] X;
  int<lower=0> y[N];
  vector[N] stl;
}
parameters{
  real<lower=0> theta_f;
  real<lower=0> theta_m;
  real<lower=0> theta_f_stl;
  real<lower=0> theta_m_stl;
  real<lower=-1,upper=1> kappa;
  real<lower=-1,upper=1> kappa_stl;
  real<lower=0> sigma_stl;
  real<lower=0> tau;
  real<lower=0> tau_stl;
  real beta;
  vector[2*T] u;
  vector[2*T] u_stl;
}
transformed parameters{
  vector[2*T] mu;
  mu[1:T] = theta_f + u[1:T];
  mu[(T+1):(2*T)] = theta_m + u[(T+1):(2*T)];
  vector[2*T] mu_stl;
  mu_stl[1:T] = theta_f_stl + u_stl[1:T];
  mu_stl[(T+1):(2*T)] = theta_m_stl + u_stl[(T+1):(2*T)];
  vector[2*T] lambda = exp(mu);
}
model{
  // Data Model
  y ~ poisson_log(X * mu + beta * (stl - X * mu_stl) - pow(beta*sigma_stl,2)/2);
  stl ~ normal(X * mu_stl,sigma_stl);

  // OU Process for Gender and Time Specific Means
  u[2:T] ~ normal(kappa*u[1:(T-1)],tau);
  u[(T+2):(2*T)] ~ normal(kappa*u[(T+1):(2*T-1)],tau);
  u_stl[2:T] ~ normal(kappa_stl*u_stl[1:(T-1)],tau_stl);
  u_stl[(T+2):(2*T)] ~ normal(kappa_stl*u_stl[(T+1):(2*T-1)],tau_stl);

  // Initial Gender and Time Mean Model
  u[{1,T+1}] ~ normal(0,tau/sqrt(1-kappa*kappa));
  u_stl[{1,T+1}] ~ normal(0,tau_stl/sqrt(1-kappa_stl*kappa_stl));

  // Priors
  theta_f ~ normal(0,100);
  theta_m ~ normal(0,100);
  kappa ~ normal(0,1);
  tau ~ normal(0,10);
  theta_f_stl ~ normal(0,100);
  theta_m_stl ~ normal(0,100);
  kappa_stl ~ normal(0,1);
  sigma_stl ~ normal(0,10);
  tau_stl ~ normal(0,10);
  beta ~ normal(0,5);
}
"
stan_cont_dep <- stan_model(model_code = cont_dep)
stan_disc_ind <- stan_model(model_code = disc_ind)
stan_disc_dep <- stan_model(model_code = disc_dep)
M <- max(dat$imp)

vars <- c("lps","ect","tpg",
          "cle","pmx","ds1","ds2",
          "ds3","lpt","mdf","mav",
          "maf","mcv","mds","mpt")

for (m in 1:M){
  tmp <- dat %>% filter(imp == m)
  X <- model.matrix(~0+group,tmp)

  N <- nrow(tmp)
  T <- max(tmp$time)

  for (var in vars){
    print(paste("Running Imputation",m,"for variable",var))

    y <- tmp[,var]
    stl <- tmp$stl

    if(var == "mav" | var == "mcv"){
      fit <- sampling(stan_disc_dep,
                      list(N=N,T=T,X=X,y=y,stl=stl),
                      iter=1500,chains=1,warmup=1000,thin=5,
                      pars = c("theta_f","theta_m",
                               "kappa","tau",
                               "beta","lambda"))
    } else if (substr(var,1,1) == "m") {
      fit <- sampling(stan_disc_ind,
                      list(N=N,T=T,X=X,y=y),
                      iter=1500,chains=1,warmup=1000,thin=5,
                      pars = c("theta_f","theta_m",
                               "kappa","tau",
                               "lambda"))
    } else if (var == "lps") {
      fit <- sampling(stan_cont_dep,
                      list(N=N,T=T,X=X,y=y,stl=stl),
                      iter=1500,chains=1,warmup=1000,thin=5,
                      pars = c("theta_f","theta_m",
                               "kappa","sigma",
                               "tau","beta","mu",
                               "theta_f_stl","theta_m_stl",
                               "kappa_stl","sigma_stl",
                               "tau_stl","mu_stl"))
    } else {
      fit <- sampling(stan_cont_dep,
                      list(N=N,T=T,X=X,y=y,stl=stl),
                      iter=1500,chains=1,warmup=1000,thin=5,
                      pars = c("theta_f","theta_m",
                               "kappa","sigma",
                               "tau","beta","mu"))
    }

    if (var == "lps"){
      imp_fit <- as.data.frame(fit)
      imp_fit <- imp_fit %>% select(!'lp__')
      names(imp_fit)[names(imp_fit) == "theta_f"] <- paste0("theta_f_",var)
      names(imp_fit)[names(imp_fit) == "theta_m"] <- paste0("theta_m_",var)
      names(imp_fit)[names(imp_fit) == "kappa"] <- paste0("kappa_",var)
      names(imp_fit)[names(imp_fit) == "sigma"] <- paste0("sigma_",var)
      names(imp_fit)[names(imp_fit) == "tau"] <- paste0("tau_",var)
      names(imp_fit)[names(imp_fit) == "beta"] <- paste0("beta_",var)
      names(imp_fit)[grep("mu\\[",names(imp_fit))] <- paste0("mu_",var,"[",1:36,"]")
      names(imp_fit)[grep("lambda\\[",names(imp_fit))] <- paste0("lambda_",var,"[",1:36,"]")
    } else {
      imp_fit <- cbind(imp_fit,as.data.frame(fit))
      names(imp_fit)[names(imp_fit) == "theta_f"] <- paste0("theta_f_",var)
      names(imp_fit)[names(imp_fit) == "theta_m"] <- paste0("theta_m_",var)
      names(imp_fit)[names(imp_fit) == "kappa"] <- paste0("kappa_",var)
      names(imp_fit)[names(imp_fit) == "sigma"] <- paste0("sigma_",var)
      names(imp_fit)[names(imp_fit) == "tau"] <- paste0("tau_",var)
      names(imp_fit)[names(imp_fit) == "beta"] <- paste0("beta_",var)
      names(imp_fit)[grep("mu\\[",names(imp_fit))] <- paste0("mu_",var,"[",1:36,"]")
      names(imp_fit)[grep("lambda\\[",names(imp_fit))] <- paste0("lambda_",var,"[",1:36,"]")
    }
  }

  if (m == 1){
    samps <- imp_fit
  } else {
    samps <- rbind(samps,imp_fit)
  }
}

saveRDS(samps,"posterior_samples_new_imputations.RDS")
