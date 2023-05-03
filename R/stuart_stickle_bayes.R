rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)

dat <- read.csv("../Data/stuart_stickle_imp_100_raw.csv")
dt <- dat %>%
  group_by(imputation_number,time,gender) %>%
  summarise(mu = mean(stl),
            n = length(stl)) %>%
  ungroup %>%
  complete(imputation_number,time,gender,fill = list(mu = 0,n = 0))

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
  real<lower=0,upper=1> kappa;
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
  sigma ~ normal(0,2);
  tau ~ normal(0,2);
}
"
sm <- stan_model(model_code = mod, verbose = TRUE)
M <- max(dat$imputation_number)

for (m in 1:M){
  tmp <- dat %>% filter(imputation_number == m)
  y <- tmp$stl
  gender <- ifelse(tmp$gender == "Female",1,0)
  time <- tmp$time

  N <- length(y)
  T <- max(time)

  fit <- sampling(sm,
                  list(N=N,T=T,y=y,gender=gender,time=time),
                  iter=10000,chains=1,warmup=5000,thin=5,
                  pars = c("theta_f","theta_m","mu_f","mu_m"))

  if (m == 1){
    samps <- as.data.frame(fit)
  } else {
    samps <- rbind(samps,as.data.frame(fit))
  }
}


samps %>%
  rename(Female = theta_f) %>%
  rename(Male = theta_m) %>%
  select(Female,Male) %>%
  gather("Gender","Value") %>%
  ggplot() +
    geom_density(aes(x=Value,color=Gender)) +
    theme_bw() +
    xlab(expression(theta)) +
    ylab(expression(paste("p(",theta,"|",x,")")))

length(which(samps$theta_m > samps$theta_f))/nrow(samps)

mu_f <- samps %>% select(contains("mu_f"))
mu_m <- samps %>% select(contains("mu_m"))

mu_diff <- mu_m - mu_f

data.frame(time = 1:18,
           med = apply(mu_diff,2,quantile,probs=0.5),
           low = apply(mu_diff,2,quantile,probs=0.05),
           upp = apply(mu_diff,2,quantile,probs=0.95)) %>%
ggplot(aes(x=time,y=med)) +
  geom_point() +
  geom_errorbar(aes(ymin=low,ymax=upp),color = "red") +
  theme_bw() +
  xlab(expression(t)) +
  ylab(expression(paste(mu[m],"-",mu[f])))

names(mu_diff) <- c(1:18)
as.data.frame(mu_diff) %>%
  gather("Time","Value") %>%
  mutate(Time = as.numeric(Time)) %>%
  ggplot(aes(x=Value,y=Time,group=Time)) +
  geom_density_ridges(fill="blue",alpha=0.5) +
  xlim(c(-15,20)) +
  scale_y_reverse() +
  geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
  xlab(expression(paste(mu[mt],"-",mu[ft]))) +
  ylab(expression(t)) +
  theme_bw()

apply(mu_diff,2,function(x){length(which(x > 0))})/nrow(mu_diff)

