rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)

DIC_cont_ind <- function(y,X,mu,sigma){
  dnorm(y,as.vector(X %*% mu),sigma,log=TRUE) %>% sum
}

DIC_cont_dep <- function(y,X,mu,stl,gamma,sigma){
  dnorm(y,as.vector(X %*% mu) + gamma*stl,sigma,log=TRUE) %>% sum
}

DIC_disc_ind <- function(y,X,mu){
  dpois(y,exp(as.vector(X %*% mu)),log=TRUE) %>% sum
}

DIC_disc_dep <- function(y,X,mu,stl,gamma){
  dpois(y,exp(as.vector(X %*% mu) + gamma*stl),log=TRUE) %>% sum
}

dat <- read.csv("Data/updated_stickle_imputations_100.csv")
times <- read.csv("RawData_ReadOnly/KSampleMeanTimes.csv")
times$time <- c(1:18)
dat <- dat %>% left_join(times %>% select(time,Inverted.Year))
dat$Inverted.Year <- dat$Inverted.Year/1000
dat$group <- factor(paste(dat$gender,letters[dat$time]))

M <- max(dat$imp)
vars <- c("stl","lps","ect","tpg",
          "cle","pmx","ds1","ds2",
          "ds3","lpt","mdf","mav",
          "maf","mcv","mds","mpt")
stl <- matrix(dat$stl,ncol=M)

DICs <- NULL

apply(samps,)

for (var in vars){
  y <- matrix(dat[,var],ncol=M)
  for (mod in 0:3){
    samps <- readRDS(paste0("Figures/",var,"/posterior_samples_",var,"_",mod,".RDS"))

    N <- nrow(samps)

    mu <- samps %>% select(contains("mu"))
    if (substr(var,1,1) != "m"){
      sigma <- samps %>% select(sigma)
    }
    gamma <- samps %>% select(gamma)

    if (substr(var,1,1) != "m"){
      params <- cbind(mu,gamma,sigma)
    } else {
      params <- cbind(mu,gamma)
    }

    apply(params,1,function(pars){
      sapply(1:M,function(m){
        X <- model.matrix(~0+group,dat %>% filter(imp==m))
        if(var == "mav" | var == "mcv"){
          DIC_disc_dep(y[,m],X,pars[1:36],stl[,m],pars[37])
        } else if (substr(var,1,1) == "m") {
          DIC_disc_ind(y[,m],X,pars[1:36])
        } else if (var=="stl") {
          DIC_cont_ind(y[,m],X,pars[1:36],pars[38])
        } else {
          DIC_cont_dep(y[,m],X,pars[1:36],stl[,m],pars[37],pars[38])
        }
      })
    })

    D_theta <- 0
    for (i in 1:N){
      print(paste0("Running DIC for ",var," model ",mod," for iteration ",i))
      probs <- rep(NA,M)
      for (m in 1:M){
        X <- model.matrix(~0+group,dat %>% filter(imp==m))
        if(var == "mav" | var == "mcv"){
          probs[m] <- DIC_disc_dep(y[,m],X,t(mu[i,]),stl[,m],gamma[i,1])
        } else if (substr(var,1,1) == "m") {
          probs[m] <- DIC_disc_ind(y[,m],X,t(mu[i,]))
        } else if (var=="stl") {
          probs[m] <- DIC_cont_ind(y[,m],X,t(mu[i,]),sigma[i,1])
        } else {
          probs[m] <- DIC_cont_dep(y[,m],X,t(mu[i,]),stl[,m],gamma[i,1],sigma[i,1])
        }
      }
      D_theta <- D_theta + 1/N*(log(mean(exp(probs - max(probs)))) + max(probs))
    }
    for (m in 1:M){
      X <- model.matrix(~0+group,dat %>% filter(imp==m))
      if(var == "mav" | var == "mcv"){
        probs[m] <- DIC_disc_dep(y[,m],X,apply(mu,2,mean),stl[,m],mean(gamma[,1]))
      } else if (substr(var,1,1) == "m") {
        probs[m] <- DIC_disc_ind(y[,m],X,apply(mu,2,mean))
      } else if (var=="stl") {
        probs[m] <- DIC_cont_ind(y[,m],X,apply(mu,2,mean),mean(sigma[,1]))
      } else {
        probs[m] <- DIC_cont_dep(y[,m],X,apply(mu,2,mean),stl[,m],mean(gamma[,1]),mean(sigma[,1]))
      }
      probs[m] <- DIC_cont_ind(y[,m],X,apply(mu,2,mean),mean(sigma[,1]))
    }
    DICs <- c(DICs,-4*D_theta + 2*(log(mean(exp(probs - max(probs)))) + max(probs)))
  }
}
