rm(list = ls())
library(tidyverse)
library(rstan)
library(ggridges)
library(bayesplot)

DIC_cont_ind <- function(y,X,mu,sigma){
  (dnorm(y,as.vector(X %*% mu),sigma,log=TRUE) %>% sum)*-2
}

DIC_cont_dep <- function(y,X,mu,stl,gamma,sigma){
  (dnorm(y,as.vector(X %*% mu) + gamma*stl,sigma,log=TRUE) %>% sum)*-2
}

DIC_disc_ind <- function(y,X,mu){
  (dpois(y,exp(as.vector(X %*% mu)),log=TRUE) %>% sum)*-2
}

DIC_disc_dep <- function(y,X,mu,stl,gamma){
  (dpois(y,exp(as.vector(X %*% mu) + gamma*stl),log=TRUE) %>% sum)*-2
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
cont_w_zeros <- c("lps","tpg","ds1","ds2","ds3")
stl <- matrix(dat$stl,ncol=M)
N <- nrow(stl)
invlogit <- function(x){exp(x)/(1 + exp(x))}

for (var in vars){
  y <- matrix(dat[,var],ncol=M)
  if (substr(var,1,1) == "m"){
    specific_x_values <- as.vector(y) %>% unique() %>% sort()
  } else {
    specific_x_values <- seq(min(density(y)$x),max(density(y)$x),length.out=100)
  }
  num_dens <- 100
  dens <- rep(0,length(specific_x_values))
  for (m in 1:M){
    if(substr(var,1,1) == "m"){
      dens <- dens + 1/M*sapply(specific_x_values,function(x){length(which(y[,1] == x))})/N
    } else {
      kde <- density(y[,m])
      dens <- dens + 1/M*approx(kde$x,kde$y,xout = specific_x_values,method = "linear")$y
    }
  }
  dens <- data.frame(x=specific_x_values,y=dens)

  for (mod in 0:3){
    probs <- NULL
    samps <- readRDS(paste0("Figures/",var,"/posterior_samples_",var,"_",mod,".RDS"))

    mu <- samps %>% select(contains("mu"))
    if (substr(var,1,1) != "m"){
      sigma <- samps %>% select(sigma)
    }
    gamma <- samps %>% select(gamma)

    print(paste0("Calculating DIC for ",var," for model ",mod))

    for (m in 1:M){
      X <- model.matrix(~0+group,dat %>% filter(imp==m))
      for (i in (100*(m-1) + 1):(100*m)){
        if(var == "mav" | var == "mcv"){
          probs <- c(probs,DIC_disc_dep(y[,m],X,t(mu[i,]),stl[,m],gamma[i,1]))
        } else if (substr(var,1,1) == "m") {
          probs <- c(probs,DIC_disc_ind(y[,m],X,t(mu[i,])))
        } else if (var=="stl") {
          probs <- c(probs,DIC_cont_ind(y[,m],X,t(mu[i,]),sigma[i,1]))
        } else {
          probs <- c(probs,DIC_cont_dep(y[,m],X,t(mu[i,]),stl[,m],gamma[i,1],sigma[i,1]))
        }
      }
    }
    DIC <- mean(probs) + 0.5*var(probs)
  }

  DICList <- as.data.frame(matrix(DICs,nrow=1,byrow=TRUE))
  colnames(DICList) <- c("OU-Trend","No OU-Trend","OU-No Trend","No OU-No Trend")
  DICList$var <- vars
  DICList <- DICList %>% pivot_longer(cols=-var,names_to="model",values_to="DIC")

  for (mod in 0:3){
    tmp <- data.frame(x=specific_x_values)
    samps <- readRDS(paste0("Figures/",var,"/posterior_samples_",var,"_",mod,".RDS"))

    n_samps <- nrow(samps)

    mu <- samps %>% select(contains("mu"))
    if (substr(var,1,1) != "m"){
      sigma <- samps %>% select(sigma)
    }
    gamma <- samps %>% select(gamma)

    for (i in 1:num_dens){
      idx <- sample(1:n_samps,1)
      m <- (idx - 1) %/% M + 1
      X <- model.matrix(~0+group,dat %>% filter(imp==m))
      #for (j in 1:N){
      if(var == "mav" | var == "mcv"){
        y_sim <- rpois(N,exp(X %*% t(mu[idx,]) + gamma[idx,1] * stl[,m]))
      } else if (substr(var,1,1) == "m") {
        y_sim <- rpois(N,exp(X %*% t(mu[idx,])))
      } else if (var=="stl") {
        y_sim <- rnorm(N,X %*% t(mu[idx,]),sigma[idx,1])
      } else {
        y_sim <- rnorm(N,X %*% t(mu[idx,]) + gamma[idx,1] * stl[,m],sigma[idx,1])
      }
      if(substr(var,1,1) == "m"){
        tmp <- cbind(tmp,sapply(specific_x_values,function(x){length(which(y_sim == x))})/N)
      } else {
        kde <- density(y_sim)
        tmp <- cbind(tmp,dens=approx(kde$x,kde$y,xout = specific_x_values,method = "linear")$y)
      }
    }
    names(tmp) <- c("x",paste0("dens",c(1:num_dens)))
    tmp <- tmp %>%
      pivot_longer(cols=-c("x"),
                   names_to="Simulation",
                   values_to="Value") %>%
      mutate(model=case_when(mod==0 ~ "OU-Trend",
                             mod==1 ~ "No OU-Trend",
                             mod==2 ~ "OU-No Trend",
                             mod==3 ~ "No OU-No Trend"))
    assign(paste0("tmp",mod),tmp)
  }
  plot_dat <- dens %>%
    left_join(rbind(tmp0,tmp1,tmp2,tmp3)) %>%
    left_join(DICList %>% filter(var==vars[j])) %>%
    mutate(DIC = paste("DIC =",round(DIC,2)))
  if (var == "mav" | var == "mds"){
    plot_dat %>%
      ggplot() +
      geom_segment(aes(x=x,xend=x,y=0,yend=Value,group=Simulation),colour="grey",alpha=0.4,lineend = "round") +
      geom_segment(aes(x=x,xend=x,y=0,yend=y),color="red",lineend = "round") +
      geom_point(aes(x=x,y=Value,group=Simulation)) +
      geom_point(aes(x=x,y=y),color="red") +
      facet_wrap(model ~ .) +
      geom_text(aes(
        x=quantile(x,0.7),
        y=max(quantile(Value,0.975,na.rm=TRUE),
              quantile(y,0.975,na.rm=TRUE)),
        label=DIC)) +
      theme_bw() +
      xlab(var) +
      ylab(expression(paste("p(",tilde(y),"|",y,")")))
  } else if (substr(var,1,1) == "m"){
    plot_dat %>%
      ggplot() +
      geom_segment(aes(x=x,xend=x,y=0,yend=Value,group=Simulation),colour="grey",alpha=0.4,lineend = "round") +
      geom_segment(aes(x=x,xend=x,y=0,yend=y),color="red",lineend = "round") +
      geom_point(aes(x=x,y=Value,group=Simulation)) +
      geom_point(aes(x=x,y=y),color="red") +
      facet_wrap(model ~ .) +
      geom_text(aes(
        x=quantile(x,0.7),
        y=max(quantile(Value,0.975,na.rm=TRUE),
              quantile(y,0.975,na.rm=TRUE)),
        label=DIC)) +
      theme_bw() +
      xlab(var) +
      ylab(expression(paste("p(",tilde(y),"|",y,")")))
  } else {
    plot_dat %>%
      ggplot() +
      geom_line(aes(x=x,y=Value,group=Simulation),colour="grey",alpha=0.4) +
      geom_line(aes(x=x,y=y),color="red") +
      facet_wrap(model ~ .) +
      geom_text(aes(
        x=quantile(x,0.8),
        y=max(quantile(Value,0.975,na.rm=TRUE),
              quantile(y,0.975,na.rm=TRUE)),
        label=DIC)) +theme_bw() +
      xlab(var) +
      ylab(expression(paste("p(",tilde(y),"|",y,")")))
  }
  ggsave(paste0("Figures/",var,"/post_pred_check.pdf"),width=7,height=10)
}


