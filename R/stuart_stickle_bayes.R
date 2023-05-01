library(tidyverse)

dat <- read.csv("stuart_stickle_imp_100_raw.csv")
# M <- max(dat$imputation_number)
# tmp <- dat %>% filter(imputation_number == 1)
T <- max(tmp$time)
N <- nrow(tmp)

new_dat <- dat %>%
  group_by(gender,index) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from="gender",values_from = "n") %>%
  replace(is.na(.),0) %>%
  mutate(prop = Female/M) %>%
  select(index,prop) %>%
  right_join(dat %>% filter(imputation_number == 1)) %>%
  arrange(index) %>%
  select(-c(imputation_number,gender))

cat("model{

    # Define the observation likelihood function
    for (i in 1:N){
      y[i] ~ dnorm(mu[i],omega)
      mu[i] <- gender[i]*(theta_f + u_f[time[i]]) +
               (1 - gender[i])*(theta_m + u_m[time[i]])
      gender[i] ~ dbern(p[i])
    }

    # Define the gender-specific state likelihood
    for (t in 2:T){
      u_f[t] ~ dnorm(kappa*u_f[t-1],tau2)
      u_m[t] ~ dnorm(kappa*u_m[t-1],tau2)
    }
    u_f[1] ~ dnorm(0,(1-kappa^2)*tau2)
    u_m[1] ~ dnorm(0,(1-kappa^2)*tau2)

    # Define the prior distributions
    theta_f ~ dnorm(54,1.0E-4)
    theta_m ~ dnorm(54,1.0E-4)
    kappa ~ dnorm(0.5,1)
    omega ~ dgamma(1,1)
    tau2 ~ dgamma(1,1)
}", file = "jags_model.txt")

inputs <- list(N=N,T=T,y=new_dat$stl,time=new_dat$time,p=new_dat$prop)

jags.m <- jags.model(file = "jags_model.txt", data=inputs, n.chains=1, n.adapt=10000)

samps <- coda.samples(jags.m,
                      variable.names = c("theta_f","theta_m","kappa","omega","tau2"),
                      n.iter=10000)

theta_f <- samps[[1]][,4]
theta_m <- samps[[1]][,5]
length(which(theta_m > theta_f))/10000
plot(theta_f)
plot(theta_m)

length(which(theta_m > theta_f))
