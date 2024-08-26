rm(list=ls())
library(tidyverse)
library(ggridges)
library(LaplacesDemon)

vars <- c("stl","lps","ect","tpg",
          "cle","pmx","ds1","ds2",
          "ds3","lpt","mdf","mav",
          "maf","mcv","mds","mpt")

cont_w_zeros <- c("lps","tpg","ds1","ds2","ds3")

invlogit <- function(x){
  exp(x)/(1 + exp(x))
}

softmax <- function(x){
  exp(x)/sum(exp(x))
}

dat <- read.csv("../Data/updated_stickle_imputations_100.csv")
times <- read.csv("../RawData_ReadOnly/KSampleMeanTimes.csv")
years <- times$Inverted.Year/1000
T <- max(dat$time)

R <- 1000

armor <- c("ds1","ds2","ds3","lps","lpt","mpt","tpg")
nonarmor <- vars[c(!vars %in% armor)]
post_probs_armor <- NULL
post_probs_nonarmor <- NULL

slopes <- NULL
pvals <- NULL

# for (var in vars){
#   file.remove(paste0("Figures/",var,"/mu_diff.pdf"))
#   file.remove(paste0("Figures/",var,"/lambda_diff.pdf"))
#   file.remove(paste0("Figures/",var,"/mu_m.pdf"))
#   file.remove(paste0("Figures/",var,"/lambda_m.pdf"))
#   file.remove(paste0("Figures/",var,"/mu_f.pdf"))
#   file.remove(paste0("Figures/",var,"/lambda_f.pdf"))
#   file.remove(paste0("Figures/",var,"/beta.pdf"))
#   file.remove(paste0("Figures/",var,"/gamma.pdf"))
#   file.remove(paste0("Figures/",var,"/theta.pdf"))
# }

for (i in 1:length(vars)){
  phenotype = vars[i]

  if (phenotype == "tpg"){
    bestDIC = 1
  } else {
    bestDIC = 0
  }

  samps <- readRDS(paste0("../Figures/",phenotype,"/posterior_samples_",phenotype,"_",bestDIC,".RDS"))

  if (i != 1){
    samps_stl <- readRDS(paste0("../Figures/stl/posterior_samples_stl_0.RDS"))
  }

  if (phenotype == "tpg"){
    mu_hi <- samps %>% select(contains("mu_hi"))
    mu_lo <- samps %>% select(contains("mu_lo"))
    mu_prob <- samps %>% select(contains("mu_prob"))
  } else if (phenotype == "mds" | phenotype == "mpt"){
    beta0 <- samps %>% select(matches("beta0.*"))
    beta1 <- samps %>% select(matches("beta1.*"))
    u <- samps %>% select(matches("^u\\[.*"))
    beta0_f <- samps %>% select(matches("beta0\\[1,.*"))
    beta1_f <- samps %>% select(matches("beta1\\[1,.*"))
    u_f <- u %>% select(!matches(".*([1][9]|[2-3][0-9]),"))
    beta0_m <- samps %>% select(matches("beta0\\[2,.*"))
    beta1_m <- samps %>% select(matches("beta1\\[2,.*"))
    u_m <- u %>% select(matches(".*([1][9]|[2-3][0-9]),"))

    if (phenotype == "mds"){
      for (j in 1:4){
        tmp <- (beta0_f[,j] + rep(years,each = 10000) * beta1_f[,j]) %>% matrix(nrow=10000) +
          u_f[,1:18 + (j-1)*18]
        assign(paste0("softmax_f_",j),tmp)
        tmp <- (beta0_m[,j] + rep(years,each = 10000) * beta1_m[,j]) %>% matrix(nrow=10000) +
          u_m[,1:18 + (j-1)*18]
        assign(paste0("softmax_m_",j),tmp)
      }

      prob_f_1 <- exp(softmax_f_1)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4))
      prob_f_2 <- exp(softmax_f_2)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4))
      prob_f_3 <- exp(softmax_f_3)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4))
      prob_f_4 <- exp(softmax_f_4)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4))
      prob_m_1 <- exp(softmax_m_1)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4))
      prob_m_2 <- exp(softmax_m_2)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4))
      prob_m_3 <- exp(softmax_m_3)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4))
      prob_m_4 <- exp(softmax_m_4)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4))

      mu_f <- prob_f_2 + 2*prob_f_3 + 3*prob_f_4
      mu_m <- prob_m_2 + 2*prob_m_3 + 3*prob_f_4
    } else if (phenotype == "mpt"){
      for (j in 1:5){
        tmp <- (beta0_f[,j] + rep(years,each = 10000) * beta1_f[,j]) %>% matrix(nrow=10000) +
          u_f[,1:18 + (j-1)*18]
        assign(paste0("softmax_f_",j),tmp)
        tmp <- (beta0_m[,j] + rep(years,each = 10000) * beta1_m[,j]) %>% matrix(nrow=10000) +
          u_m[,1:18 + (j-1)*18]
        assign(paste0("softmax_m_",j),tmp)
      }

      prob_f_1 <- exp(softmax_f_1)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4) + exp(softmax_f_5))
      prob_f_2 <- exp(softmax_f_2)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4) + exp(softmax_f_5))
      prob_f_3 <- exp(softmax_f_3)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4) + exp(softmax_f_5))
      prob_f_4 <- exp(softmax_f_4)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4) + exp(softmax_f_5))
      prob_f_5 <- exp(softmax_f_5)/(exp(softmax_f_1) + exp(softmax_f_2) + exp(softmax_f_3) + exp(softmax_f_4) + exp(softmax_f_5))
      prob_m_1 <- exp(softmax_m_1)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4) + exp(softmax_m_5))
      prob_m_2 <- exp(softmax_m_2)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4) + exp(softmax_m_5))
      prob_m_3 <- exp(softmax_m_3)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4) + exp(softmax_m_5))
      prob_m_4 <- exp(softmax_m_4)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4) + exp(softmax_m_5))
      prob_m_5 <- exp(softmax_m_5)/(exp(softmax_m_1) + exp(softmax_m_2) + exp(softmax_m_3) + exp(softmax_m_4) + exp(softmax_m_5))

      mu_f <- 3*prob_f_1 + 4*prob_f_2 + 5*prob_f_3 + 6*prob_f_4 + 7*prob_f_5
      mu_m <- 3*prob_m_1 + 4*prob_m_2 + 5*prob_m_3 + 6*prob_m_4 + 7*prob_m_5
    }
  } else {
    mu <- samps %>% select(contains("mu"))
  }
  if (i != 1){
    mu_stl <- samps_stl %>% select(contains("mu"))
    sigma_stl <- samps_stl %>% select(sigma)
    if (phenotype == "tpg"){
      gamma_hi <- samps %>% select(gamma_hi)
      gamma_lo <- samps %>% select(gamma_lo)
    } else if (phenotype != "mds" & phenotype != "mpt"){
      gamma <- samps %>% select(gamma)
    }
  }

  if (phenotype == "mav" | phenotype == "mcv"){
    mu_f <- exp(mu[,1:T] + as.data.frame(rep(gamma,18)) * mu_stl[,1:T] + as.data.frame(rep((gamma*sigma_stl)^2,18))/2)
    mu_m <- exp(mu[,(T+1):(2*T)] + as.data.frame(rep(gamma,18)) * mu_stl[,(T+1):(2*T)] + as.data.frame(rep((gamma*sigma_stl)^2,18))/2)
  } else if (phenotype == "mdf" | phenotype == "maf"){
    mu_f <- exp(mu[,1:T])
    mu_m <- exp(mu[,(T+1):(2*T)])
  } else if (phenotype == "tpg"){
    mu_f <- (mu_hi[,1:T] + as.data.frame(rep(gamma_hi,18)) * mu_stl[,1:T]) * invlogit(mu_prob[,1:T]) +
      (mu_lo[,1:T]  + as.data.frame(rep(gamma_lo,18)) * mu_stl[,1:T]) * (1 - invlogit(mu_prob[,1:T]))
    mu_m <- mu_hi[,(T+1):(2*T)] * invlogit(mu_prob[,(T+1):(2*T)]) +
      mu_lo[,(T+1):(2*T)] * (1 - invlogit(mu_prob[,(T+1):(2*T)]))
  } else if (i  %in% c(2:10)){
    mu_f <- mu[,1:T] + as.data.frame(rep(gamma,18)) * mu_stl[,1:T]
    mu_m <- mu[,(T+1):(2*T)] + as.data.frame(rep(gamma,18)) * mu_stl[,(T+1):(2*T)]
  } else if (phenotype == "stl"){
    mu_f <- mu[,1:T]
    mu_m <- mu[,(T+1):(2*T)]
  }

  mu_diff_med <- apply(mu_m - mu_f,2,median)

  # if(i > 10){
  #   mu_diff <- exp(mu_diff)
  # }

  # mulab <- if (i <= 10){
  #   ylab(expression(paste(mu[mt],"-",mu[ft])))
  # } else {
  #   ylab(expression(paste(lambda[mt],"-",lambda[ft])))
  # }

  # mu_flab <- if (i <= 10){
  #   xlab(expression(mu[ft]))
  # } else {
  #   xlab(expression(lambda[ft]))
  # }
  #
  # mu_mlab <- if (i <= 10){
  #   xlab(expression(mu[mt]))
  # } else {
  #   xlab(expression(lambda[mt]))
  # }

  # best_model <- case_when(bestDIC == 0 ~ "OU - Trend",
  #                         bestDIC == 1 ~ "No OU - Trend",
  #                         bestDIC == 2 ~ "OU - No Trend",
  #                         bestDIC == 3 ~ "No OU - No Trend")

  plot_dat <- data.frame(years,
                    mu_diff_med)
  lin_mod <- lm(mu_diff_med ~ years,plot_dat)
  curr_plot <- ggplot(plot_dat,aes(x=years,y=mu_diff_med)) +
    geom_point(aes(x=years,mu_diff_med)) +
    geom_smooth(method="lm",se=FALSE) +
    # geom_hline(aes(yintercept = 0),colour = "red",linetype=2) +
    ylab(expression(paste(mu[mt],"-",mu[ft]))) +
    xlab(expression(t)) +
    ggtitle(paste(phenotype,": ",
                   expression(beta[1])," = ",
                   summary(lin_mod)$coefficients[2,1] %>% signif(4),
                   ", p = ",
                   summary(lin_mod)$coefficients[2,4] %>% signif(4))) +
    theme_bw()

  if (max(plot_dat$mu_diff_med) > 0 & min(plot_dat$mu_diff_med) < 0){
    final_plot <- curr_plot +
      geom_hline(aes(yintercept = 0),colour = "red",linetype=2)
  } else if (max(plot_dat$mu_diff_med) < 0 & min(plot_dat$mu_diff_med) > 0){
    final_plot <- curr_plot +
      geom_hline(aes(yintercept = 0),colour = "red",linetype=2)
  } else {
    final_plot <- curr_plot
  }

  ggsave(paste0("../Figures/SD_Figures/",phenotype,".pdf"),
         final_plot,
         width=5,height=8)

  slopes <- c(slopes,summary(lm(mu_diff_med ~ years,plot_dat))$coefficients[2,1])
  pvals <- c(pvals,summary(lm(mu_diff_med ~ years,plot_dat))$coefficients[2,4])

  # names(mu_f) <- c(1:18)
  # as.data.frame(mu_f) %>%
  #   gather("Time","Value") %>%
  #   mutate(Time = as.numeric(Time)) %>%
  #   ggplot(aes(x=Value,y=Time,group=Time)) +
  #   geom_density_ridges(fill="red",alpha=0.5) +
  #   scale_y_reverse() +
  #   mu_flab +
  #   ylab(expression(t)) +
  #   theme_bw()
  # if (i <= 10){
  #   ggsave(paste0("Figures/",phenotype,"/mu_f.pdf"),width=5,height=8)
  # } else {
  #   ggsave(paste0("Figures/",phenotype,"/lambda_f.pdf"),width=5,height=8)
  # }
  #
  # names(mu_m) <- c(1:18)
  # as.data.frame(mu_m) %>%
  #   gather("Time","Value") %>%
  #   mutate(Time = as.numeric(Time)) %>%
  #   ggplot(aes(x=Value,y=Time,group=Time)) +
  #   geom_density_ridges(fill="blue",alpha=0.5) +
  #   scale_y_reverse() +
  #   mu_mlab +
  #   ylab(expression(t)) +
  #   theme_bw()
  # if (i <= 10){
  #   ggsave(paste0("Figures/",phenotype,"/mu_m.pdf"),width=5,height=8)
  # } else {
  #   ggsave(paste0("Figures/",phenotype,"/lambda_m.pdf"),width=5,height=8)
  # }

#   if (i >=2 & i <= 10){
#     samps %>%
#       rename(beta = paste0("beta_",phenotype))%>%
#       select(beta) %>%
#       ggplot() +
#       geom_density(aes(x=beta)) +
#       geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
#       theme_bw() +
#       xlab(expression(beta[g])) +
#       ylab(expression(paste("p(",beta[g],"|",bold(X),")")))
#     ggsave(paste0("Figures/",phenotype,"/beta.pdf"),width=5,height=8)
#   }
#   if (phenotype == "mav"){
#     samps %>%
#       rename(beta = paste0("beta_",phenotype))%>%
#       select(beta) %>%
#       ggplot() +
#       geom_density(aes(x=beta)) +
#       geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
#       theme_bw() +
#       xlab(expression(beta[g])) +
#       ylab(expression(paste("p(",beta[g],"|",bold(X),")")))
#     ggsave(paste0("Figures/",phenotype,"/beta.pdf"),width=5,height=8)
#   }
#   if (phenotype == "mcv"){
#     samps %>%
#       rename(beta = paste0("beta_",phenotype))%>%
#       select(beta) %>%
#       ggplot() +
#       geom_density(aes(x=beta)) +
#       geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
#       theme_bw() +
#       xlab(expression(beta[g])) +
#       ylab(expression(paste("p(",beta[g],"|",bold(X),")")))
#     ggsave(paste0("Figures/",phenotype,"/beta.pdf"),width=5,height=8)
#   }
# }
#
  # if (phenotype %in% armor){
  #   post_probs_armor <- cbind(post_probs_armor,
  #                         apply(mu_diff,2,function(x){length(which(x > 0))/length(x)}))
  # } else {
  #   post_probs_nonarmor <- cbind(post_probs_nonarmor,
  #                            apply(mu_diff,2,function(x){length(which(x > 0))/length(x)}))
  # }


# for (i in 1:length(armor)){
#   if (substr(armor[i],1,1) == "m"){
#     mu <- samps %>% select(contains(paste0("lambda_",armor[i])))
#   } else {
#     mu <- samps %>% select(contains(paste0("mu_",armor[i])))
#   }
#   mu_f <- mu[,1:T]
#   mu_m <- mu[,(T+1):(2*T)]
#
#   mu_diff <- mu_m - mu_f
#
#   if (i == 1){
#     probs <- apply(mu_diff,2,function(x){length(which(x > 0))/length(x)})
#   } else {
#     probs <- cbind(probs,
#              apply(mu_diff,2,function(x){length(which(x > 0))/length(x)}))
#   }
# }
#
# colnames(probs) <- armor
# as.data.frame(probs) %>%
#   gather("Phenotype","Value") %>%
#   mutate(Time = rep(1:T,length(armor))) %>%
#   ggplot() +
#   geom_line(aes(x=Time,y=Value,linetype=Phenotype,colour=Phenotype)) +
#   theme_bw() +
#   ylab(expression(paste("P(",mu[mt],"-",mu[ft],"|",bold(y),")")))
# ggsave(paste0("Figures/post_probs_armor.pdf"),width=5,height=8)
#
# for (i in 1:length(nonarmor)){
#   if (substr(nonarmor[i],1,1) == "m"){
#     mu <- samps %>% select(contains(paste0("lambda_",nonarmor[i])))
#   } else {
#     mu <- samps %>% select(contains(paste0("mu_",nonarmor[i])))
#   }
#   mu_f <- mu[,1:T]
#   mu_m <- mu[,(T+1):(2*T)]
#
#   mu_diff <- mu_m - mu_f
#
#   if (i == 1){
#     probs <- apply(mu_diff,2,function(x){length(which(x > 0))/length(x)})
#   } else {
#     probs <- cbind(probs,
#                    apply(mu_diff,2,function(x){length(which(x > 0))/length(x)}))
#   }
}
#
# colnames(post_probs_armor) <- armor
# as.data.frame(post_probs_armor) %>%
#   gather("Phenotype","Value") %>%
#   mutate(Time = rep(1:T,length(armor))) %>%
#   ggplot() +
#   geom_line(aes(x=Time,y=Value,linetype=Phenotype,colour=Phenotype)) +
#   theme_bw() +
#   ylab(expression(paste("P(",mu[mt],"-",mu[ft],"|",bold(y),")")))
# ggsave(paste0("Figures/post_probs_armor.pdf"),width=5,height=8)
#
# colnames(post_probs_nonarmor) <- nonarmor
# as.data.frame(post_probs_nonarmor) %>%
#   gather("Phenotype","Value") %>%
#   mutate(Time = rep(1:T,length(nonarmor))) %>%
#   ggplot() +
#   geom_line(aes(x=Time,y=Value,linetype=Phenotype,colour=Phenotype)) +
#   theme_bw() +
#   ylab(expression(paste("P(",mu[mt],"-",mu[ft],"|",bold(y),")")))
# ggsave(paste0("Figures/post_probs_nonarmor.pdf"),width=5,height=8)
#
# for (i in 1:length(cont_w_zeros)){
#   phenotype = cont_w_zeros[i]
#
#   bestDIC = mods_zeros[i]
#
#   samps <- readRDS(paste0("Figures/",phenotype,"/posterior_samples_",phenotype,"prob_",bestDIC,".RDS"))
#   samps_stl <- readRDS(paste0("Figures/stl/posterior_samples_stl_",mods[1],".RDS"))
#
#   N <- nrow(samps)
#
#   mu <- samps %>% select(contains("mu"))
#
#   mu_stl <- samps_stl %>% select(contains("mu"))
#   sigma_stl <- samps_stl %>% select(sigma)
#   gamma <- samps %>% select(gamma)
#
#   p <- sapply(1:N,function(i){
#     stl_sim <- matrix(rnorm(R*2*T,as.vector(unlist(mu_stl[1,])),sigma_stl[1,]),ncol=2*T,byrow=TRUE)
#     apply(invlogit(matrix(unlist(mu[i,]),nrow=R,ncol=2*T,byrow=TRUE) + gamma[i,1] * stl_sim),2,mean)
#   })
#
#   p <- 1 - p
#
#   p_f <- p[1:T,]
#   p_m <- p[(T+1):(2*T),]
#
#   p_diff <- p_m - p_f
#
#   plab <- xlab(expression(paste(p[mt],"-",p[ft])))
#
#   rownames(p_diff) <- c(1:18)
#   rownames(p_m) <- c(1:18)
#   rownames(p_f) <- c(1:18)
#
#   as.data.frame(t(p_diff)) %>%
#     gather("Time","Value") %>%
#     mutate(Time = as.numeric(Time)) %>%
#     ggplot(aes(x=Value,y=Time,group=Time)) +
#     geom_density_ridges(fill="blue",alpha=0.5) +
#     scale_y_reverse() +
#     geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
#     plab +
#     ylab(expression(t)) +
#     ggtitle(paste0(phenotype,": ",best_model)) +
#     theme_bw()
#   ggsave(paste0("Figures/",phenotype,"/p_diff.pdf"),width=5,height=8)
#
#   as.data.frame(t(p_m)) %>%
#     gather("Time","Value") %>%
#     mutate(Time = as.numeric(Time)) %>%
#     ggplot(aes(x=Value,y=Time,group=Time)) +
#     geom_density_ridges(fill="blue",alpha=0.5) +
#     scale_y_reverse() +
#     xlab(expression(paste(p[mt]))) +
#     ylab(expression(t)) +
#     ggtitle(paste0(phenotype,": ",best_model)) +
#     theme_bw()
#   ggsave(paste0("Figures/",phenotype,"/p_m.pdf"),width=5,height=8)
#
#   as.data.frame(t(p_f)) %>%
#     gather("Time","Value") %>%
#     mutate(Time = as.numeric(Time)) %>%
#     ggplot(aes(x=Value,y=Time,group=Time)) +
#     geom_density_ridges(fill="red",alpha=0.5) +
#     scale_y_reverse() +
#     xlab(expression(paste(p[ft]))) +
#     ylab(expression(t)) +
#     ggtitle(paste0(phenotype,": ",best_model)) +
#     theme_bw()
#   ggsave(paste0("Figures/",phenotype,"/p_f.pdf"),width=5,height=8)
# }
for (i in 1:length(cont_w_zeros)){
  phenotype = cont_w_zeros[i]

  samps <- readRDS(paste0("../Figures/",phenotype,"/posterior_samples_",phenotype,"prob_0.RDS"))
  samps_stl <- readRDS(paste0("../Figures/stl/posterior_samples_stl_0.RDS"))

  N <- nrow(samps)

  mu <- samps %>% select(contains("mu"))
  mu_stl <- samps_stl %>% select(contains("mu"))
  sigma_stl <- samps_stl %>% select(sigma)
  gamma <- samps %>% select(gamma)

  p <- sapply(1:N,function(i){
    stl_sim <- matrix(rnorm(R*2*T,as.vector(unlist(mu_stl[1,])),sigma_stl[1,]),ncol=2*T,byrow=TRUE)
    apply(invlogit(matrix(unlist(mu[i,]),nrow=R,ncol=2*T,byrow=TRUE) + gamma[i,1] * stl_sim),2,mean)
  })

  p <- 1 - p

  p_f <- p[1:T,]
  p_m <- p[(T+1):(2*T),]

  p_diff_med <- apply(p_m - p_f,1,median)

  plot_dat <- data.frame(years,
                         p_diff_med)
  lin_mod <- lm(p_diff_med ~ years,plot_dat)
  curr_plot <- ggplot(plot_dat,aes(x=years,y=p_diff_med)) +
    geom_point(aes(x=years,p_diff_med)) +
    geom_smooth(method="lm",se=FALSE) +
    # geom_hline(aes(yintercept = 0),colour = "red",linetype=2) +
    ylab(expression(paste(p[mt],"-",p[ft]))) +
    xlab(expression(t)) +
    ggtitle(paste0(phenotype,": ",
                   expression(beta[1])," = ",
                   summary(lin_mod)$coefficients[2,1] %>% signif(4),
                   ", p = ",
                   summary(lin_mod)$coefficients[2,4] %>% signif(4))) +
    theme_bw()

  if (max(plot_dat$p_diff_med) > 0 & min(plot_dat$p_diff_med) < 0){
    final_plot <- curr_plot +
      geom_hline(aes(yintercept = 0),colour = "red",linetype=2)
  } else if (max(plot_dat$p_diff_med) < 0 & min(plot_dat$p_diff_med) > 0){
    final_plot <- curr_plot +
      geom_hline(aes(yintercept = 0),colour = "red",linetype=2)
  } else {
    final_plot <- curr_plot
  }

  ggsave(paste0("../Figures/SD_Figures/",phenotype,"_prob_zero.pdf"),
         final_plot,
         width=5,height=8)

  slopes <- c(slopes,summary(lm(p_diff_med ~ years,plot_dat))$coefficients[2,1])
  pvals <- c(pvals,summary(lm(p_diff_med ~ years,plot_dat))$coefficients[2,4])
}
