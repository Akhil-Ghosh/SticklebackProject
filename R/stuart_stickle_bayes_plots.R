library(tidyverse)
library(ggridges)

vars <- c("stl","lps","ect","tpg",
          "cle","pmx","ds1","ds2",
          "ds3","lpt","mdf","mav",
          "maf","mcv","mds","mpt")

samps <- readRDS("posterior_samples_new_imputations.RDS")

samps <- samps %>% select(!contains("lp_"))

for (i in 1:16){
  phenotype = vars[i]

  thetalab <- if (i <= 10){
    xlab(expression(theta[g]))
  } else {
    xlab(expression(gamma[g]))
  }

  samps %>%
    rename(Female = paste0("theta_f_",phenotype)) %>%
    rename(Male = paste0("theta_m_",phenotype)) %>%
    select(Female,Male) %>%
    gather("Gender","Value") %>%
    ggplot() +
    geom_density(aes(x=Value,color=Gender)) +
    scale_color_manual(values = c("red","blue"))+
    theme_bw() +
    thetalab +
    ylab(expression(paste("p(",theta[g],"|",bold(X),")")))
  if (i <= 10){
    ggsave(paste0("Figures/",phenotype,"/theta.pdf"),width=5,height=8)
  } else {
    ggsave(paste0("Figures/",phenotype,"/gamma.pdf"),width=5,height=8)
  }

  if (i <= 10){
    mu <- samps %>% select(contains(paste0("mu_",phenotype)))
  } else {
    mu <- samps %>% select(contains(paste0("lambda_",phenotype)))
  }
  mu_f <- mu[,1:T]
  mu_m <- mu[,(T+1):(2*T)]

  mu_diff <- mu_m - mu_f

  mulab <- if (i <= 10){
    xlab(expression(paste(mu[mt],"-",mu[ft])))
  } else {
    xlab(expression(paste(lambda[mt],"-",lambda[ft])))
  }

  mu_flab <- if (i <= 10){
    xlab(expression(mu[ft]))
  } else {
    xlab(expression(lambda[ft]))
  }

  mu_mlab <- if (i <= 10){
    xlab(expression(mu[mt]))
  } else {
    xlab(expression(lambda[mt]))
  }

  names(mu_diff) <- c(1:18)
  as.data.frame(mu_diff) %>%
    gather("Time","Value") %>%
    mutate(Time = as.numeric(Time)) %>%
    ggplot(aes(x=Value,y=Time,group=Time)) +
    geom_density_ridges(fill="blue",alpha=0.5) +
    scale_y_reverse() +
    geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
    mulab +
    ylab(expression(t)) +
    theme_bw()
  if (i <= 10){
    ggsave(paste0("Figures/",phenotype,"/mu_diff.pdf"),width=5,height=8)
  } else {
    ggsave(paste0("Figures/",phenotype,"/lambda_diff.pdf"),width=5,height=8)
  }

  names(mu_f) <- c(1:18)
  as.data.frame(mu_f) %>%
    gather("Time","Value") %>%
    mutate(Time = as.numeric(Time)) %>%
    ggplot(aes(x=Value,y=Time,group=Time)) +
    geom_density_ridges(fill="red",alpha=0.5) +
    scale_y_reverse() +
    mu_flab +
    ylab(expression(t)) +
    theme_bw()
  if (i <= 10){
    ggsave(paste0("Figures/",phenotype,"/mu_f.pdf"),width=5,height=8)
  } else {
    ggsave(paste0("Figures/",phenotype,"/lambda_f.pdf"),width=5,height=8)
  }

  names(mu_m) <- c(1:18)
  as.data.frame(mu_m) %>%
    gather("Time","Value") %>%
    mutate(Time = as.numeric(Time)) %>%
    ggplot(aes(x=Value,y=Time,group=Time)) +
    geom_density_ridges(fill="blue",alpha=0.5) +
    scale_y_reverse() +
    mu_mlab +
    ylab(expression(t)) +
    theme_bw()
  if (i <= 10){
    ggsave(paste0("Figures/",phenotype,"/mu_m.pdf"),width=5,height=8)
  } else {
    ggsave(paste0("Figures/",phenotype,"/lambda_m.pdf"),width=5,height=8)
  }

  if (i >=2 & i <= 10){
    samps %>%
      rename(beta = paste0("beta_",phenotype))%>%
      select(beta) %>%
      ggplot() +
      geom_density(aes(x=beta)) +
      geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
      theme_bw() +
      xlab(expression(beta[g])) +
      ylab(expression(paste("p(",beta[g],"|",bold(X),")")))
    ggsave(paste0("Figures/",phenotype,"/beta.pdf"),width=5,height=8)
  }
  if (phenotype == "mav"){
    samps %>%
      rename(beta = paste0("beta_",phenotype))%>%
      select(beta) %>%
      ggplot() +
      geom_density(aes(x=beta)) +
      geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
      theme_bw() +
      xlab(expression(beta[g])) +
      ylab(expression(paste("p(",beta[g],"|",bold(X),")")))
    ggsave(paste0("Figures/",phenotype,"/beta.pdf"),width=5,height=8)
  }
  if (phenotype == "mcv"){
    samps %>%
      rename(beta = paste0("beta_",phenotype))%>%
      select(beta) %>%
      ggplot() +
      geom_density(aes(x=beta)) +
      geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
      theme_bw() +
      xlab(expression(beta[g])) +
      ylab(expression(paste("p(",beta[g],"|",bold(X),")")))
    ggsave(paste0("Figures/",phenotype,"/beta.pdf"),width=5,height=8)
  }
}

armor <- c("ds1","ds2","ds3","lps","lpt","mpt","tpg")
nonarmor <- vars[c(!vars %in% armor)]

for (i in 1:length(armor)){
  if (substr(armor[i],1,1) == "m"){
    mu <- samps %>% select(contains(paste0("lambda_",armor[i])))
  } else {
    mu <- samps %>% select(contains(paste0("mu_",armor[i])))
  }
  mu_f <- mu[,1:T]
  mu_m <- mu[,(T+1):(2*T)]

  mu_diff <- mu_m - mu_f

  if (i == 1){
    probs <- apply(mu_diff,2,function(x){length(which(x > 0))/length(x)})
  } else {
    probs <- cbind(probs,
             apply(mu_diff,2,function(x){length(which(x > 0))/length(x)}))
  }
}

colnames(probs) <- armor
as.data.frame(probs) %>%
  gather("Phenotype","Value") %>%
  mutate(Time = rep(1:T,length(armor))) %>%
  ggplot() +
  geom_line(aes(x=Time,y=Value,linetype=Phenotype,colour=Phenotype)) +
  theme_bw() +
  ylab(expression(paste("P(",mu[mt],"-",mu[ft],"|",bold(y),")")))
ggsave(paste0("Figures/post_probs_armor.pdf"),width=5,height=8)

for (i in 1:length(nonarmor)){
  if (substr(nonarmor[i],1,1) == "m"){
    mu <- samps %>% select(contains(paste0("lambda_",nonarmor[i])))
  } else {
    mu <- samps %>% select(contains(paste0("mu_",nonarmor[i])))
  }
  mu_f <- mu[,1:T]
  mu_m <- mu[,(T+1):(2*T)]

  mu_diff <- mu_m - mu_f

  if (i == 1){
    probs <- apply(mu_diff,2,function(x){length(which(x > 0))/length(x)})
  } else {
    probs <- cbind(probs,
                   apply(mu_diff,2,function(x){length(which(x > 0))/length(x)}))
  }
}

colnames(probs) <- nonarmor
as.data.frame(probs) %>%
  gather("Phenotype","Value") %>%
  mutate(Time = rep(1:T,length(nonarmor))) %>%
  ggplot() +
  geom_line(aes(x=Time,y=Value,linetype=Phenotype,colour=Phenotype)) +
  theme_bw() +
  ylab(expression(paste("P(",mu[mt],"-",mu[ft],"|",bold(y),")")))
ggsave(paste0("Figures/post_probs_nonarmor.pdf"),width=5,height=8)

