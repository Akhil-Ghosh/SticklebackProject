library(tidyverse)
library(ggridges)

vars <- c("stl","cle","pmx")

for (x in vars){
  samps <- readRDS(paste0("posterior_samples_",x,".RDS"))

  xlim <- if (x == "stl"){
    xlim(30,80)
  } else if (x == "cle"){
    xlim(5,12)
  } else {
    xlim(25,85)
  }

  samps %>%
    rename(Female = theta_f) %>%
    rename(Male = theta_m) %>%
    select(Female,Male) %>%
    gather("Gender","Value") %>%
    ggplot() +
    geom_density(aes(x=Value,color=Gender)) +
    scale_color_manual(values = c("red","blue"))+
    theme_bw() +
    xlim +
    xlab(expression(theta[g])) +
    ylab(expression(paste("p(",theta[g],"|",bold(X),")")))
  ggsave(paste0("Figures/",x,"/theta.pdf"),width=5,height=8)

  length(which(samps$theta_m > samps$theta_f))/nrow(samps)

  mu_f <- samps %>% select(contains("mu_f"))
  mu_m <- samps %>% select(contains("mu_m"))

  mu_diff <- mu_m - mu_f

  xlim <- if (x == "stl"){
    xlim(-15,25)
  } else if (x == "cle"){
    xlim(-1.75,3.75)
  } else {
    xlim(-15,25)
  }

  names(mu_diff) <- c(1:18)
  as.data.frame(mu_diff) %>%
    gather("Time","Value") %>%
    mutate(Time = as.numeric(Time)) %>%
    ggplot(aes(x=Value,y=Time,group=Time)) +
    geom_density_ridges(fill="blue",alpha=0.5) +
    # xlim(c(-15,20)) +
    scale_y_reverse() +
    geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
    xlab(expression(paste(mu[mt],"-",mu[ft]))) +
    ylab(expression(t)) +
    theme_bw()
  ggsave(paste0("Figures/",x,"/mu_diff.pdf"),width=5,height=8)

  names(mu_f) <- c(1:18)
  as.data.frame(mu_f) %>%
    gather("Time","Value") %>%
    mutate(Time = as.numeric(Time)) %>%
    ggplot(aes(x=Value,y=Time,group=Time)) +
    geom_density_ridges(fill="red",alpha=0.5) +
    scale_y_reverse() +
    xlab(expression(paste(mu[ft]))) +
    ylab(expression(t)) +
    theme_bw()
  ggsave(paste0("Figures/",x,"/mu_f.pdf"),width=5,height=8)

  names(mu_m) <- c(1:18)
  as.data.frame(mu_m) %>%
    gather("Time","Value") %>%
    mutate(Time = as.numeric(Time)) %>%
    ggplot(aes(x=Value,y=Time,group=Time)) +
    geom_density_ridges(fill="blue",alpha=0.5) +
    scale_y_reverse() +
    xlab(expression(paste(mu[mt]))) +
    ylab(expression(t)) +
    theme_bw()
  ggsave(paste0("Figures/",x,"/mu_m.pdf"),width=5,height=8)
}


