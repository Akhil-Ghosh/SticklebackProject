rm(list=ls())
library(tidyverse)
library(ggridges)
library(LaplacesDemon)
library(officer)
library(flextable)
library(xtable)

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
times_for_join <- data.frame(time = c(1:T),
                             years = times$Inverted.Year/1000)

summ_prop_tbl <- dat %>%
  left_join(times_for_join,by = "time") %>%
  group_by(time,years,imp) %>%
  summarise(samp_size = n(),
            n_male = sum(gender=="Male"),
            n_female = sum(gender=="Female")) %>%
  group_by(time,years) %>%
  summarise(`Sample Size` = min(samp_size),
            `Number of Males (95% CI)` = paste0(round(median(n_male),2),
                                    " (",
                                    round(quantile(n_male,0.025),2),
                                    ",",
                                    round(quantile(n_male,0.975),2),
                                    ")"),
            `Number of Females (95% CI)` = paste0(round(median(n_female),2),
                                                " (",
                                                round(quantile(n_female,0.025),2),
                                                ",",
                                                round(quantile(n_female,0.975),2),
                                                ")")) %>%
  ungroup() %>%
  rename(Time = time,Year = years)
doc <- read_docx()
ft <- flextable(summ_prop_tbl)
doc %>%
  body_add_flextable(ft)
print(doc, target = "summary_prop_males.docx")

tau_dat <- dat %>%
  left_join(times_for_join,by = "time") %>%
  group_by(years,imp) %>%
  summarise(success = sum(gender == "Male"),
            failure = sum(gender == "Female"),
            .groups = "drop") %>%
  uncount(100) %>%
  mutate(iter = row_number() %% 100,
         post_prop = rbeta(n = n(),0.5+success,0.5+failure)) %>%
  group_by(imp,iter) %>%
  summarise(tau = cor(years,post_prop,method = "kendall"),
            .groups = "drop")
p_tau_prob <- length(which(tau_dat$tau > 0))/10000
ggplot(tau_dat) +
  geom_histogram(aes(x=tau)) +
  xlab(expression(tau)) + ylab("") + xlim(-1,1) +
  ggtitle("") +
  geom_vline(aes(xintercept=mean(tau))) +
  geom_vline(aes(xintercept=0),colour="red",linetype=2)+
  annotate("text",
         x = -1,
         y = Inf,
         label = sprintf("P(tau > 0) == %.4f", p_tau_prob),
         parse = TRUE,
         hjust = 0,
         vjust = 1.5) +
  theme_bw(base_size = 12)
ggsave("../Figures/Figure_3_tau.pdf",width=5,height=3)

dat %>%
  left_join(times_for_join,by = "time") %>%
  group_by(years,imp) %>%
  summarise(prop_male = mean(gender == "Male"),
            .groups = "drop") %>%
  group_by(years) %>%
  mutate(median_prop = median(prop_male)) %>%
  ggplot(aes(x = prop_male, y = years, group = years)) +
  stat_density_ridges(quantile_lines = TRUE,quantiles=2,scale = 1.2, rel_min_height = 0.01,alpha=0.3,fill="blue") +
  xlab("Male Proportion") + ylab("Relative Time of Deposition (1000s of years)") +
  coord_flip() +
  theme_bw()
ggsave("../Figures/Figure_3.pdf",width=8,height=7)

R <- 1000

post_probs <- data.frame(Phenotype = NULL,
                         Time = NULL,
                         Value = NULL)
post_probs_zero <- data.frame(Phenotype = NULL,
                         Time = NULL,
                         Value = NULL)
armor <- c("ds1","ds2","ds3","lps","lpt","mpt","tpg")
nonarmor <- vars[c(!vars %in% armor)]

means <- NULL
sds <- NULL
taus <- NULL
tau_prob <- NULL

means_zero <- NULL
sds_zero <- NULL
taus_zero <- NULL
tau_prob_zero <- NULL

tau_means <- NULL
tau_probs <- NULL

for (i in 1:length(vars)){
  print(i)
  phenotype = vars[i]
  var_dat <- dat %>%
    select(!ends_with(".sc")) %>%
    left_join(times_for_join,by = "time") %>%
    select(years,imp,contains(phenotype)) %>%
    rename(trait = all_of(phenotype)) %>%
    group_by(years,imp) %>%
    summarise(size = n(),
              xbar = mean(trait),
              s2 = var(trait),
              .groups = "drop") %>%
    uncount(100) %>%
    mutate(iter = row_number() %% 100,
           post_var = 1/rgamma(n = n(),
                             shape = 1 + size/2,
                             rate = 1 + (size-1)*s2/2 + 0*size/(size + 0)*xbar^2/2))

  var_scat_plot <- var_dat %>%
    group_by(years) %>%
    summarise(med = median(post_var),
              low = quantile(post_var,0.025),
              upp = quantile(post_var,0.975)) %>%
    ggplot() +
    geom_point(aes(x=years,med)) +
    geom_segment(aes(x=years,xend=years,y=low,yend=upp)) +
    ylab(expression(s^2)) +
    xlab("Time") +
    ggtitle(paste("Posterior sample variance for",phenotype)) +
    theme_bw(base_size = 12)

  ggsave(paste0("../Figures/SampVarScatterPlots/",phenotype,".pdf"),
         var_scat_plot,
         width=5,height=8)

  tau_dat <- var_dat %>%
    group_by(imp,iter) %>%
    summarise(tau = cor(years,post_var,method = "kendall"),
              .groups = "drop")
  tau_means <- c(tau_means,mean(tau_dat$tau))
  p_tau_prob <- length(which(tau_dat$tau > 0))/10000
  tau_probs <- c(tau_probs,p_tau_prob)
  tau_plot <- ggplot(tau_dat) +
    geom_histogram(aes(x=tau)) +
    xlab(expression(tau)) + ylab("") + xlim(-1,1) +
    ggtitle(phenotype) +
    geom_vline(aes(xintercept=mean(tau))) +
    geom_vline(aes(xintercept=0),colour="red",linetype=2)+
    annotate("text",
             x = -1,
             y = Inf,
             label = sprintf("P(tau > 0) == %.4f", p_tau_prob),
             parse = TRUE,
             hjust = 0,
             vjust = 1.5) +
    theme_bw(base_size = 12)
  ggsave(paste0("../Figures/TauHistograms/",phenotype,"SampVar.pdf"),
        tau_plot,
        width=5,height=3)
}

doc <- read_docx()
ft <- flextable(cbind(vars,round(tau_means,4),tau_probs) %>%
                  as.data.frame)
doc %>%
  body_add_flextable(ft)
print(doc, target = "summary_table.docx")
for (i in 1:length(vars)){
  phenotype = vars[i]
  samps <- readRDS(paste0("../Figures/posterior_samples_",phenotype,"_0.RDS"))

  if (i != 1){
    samps_stl <- readRDS(paste0("../Figures/posterior_samples_stl_0.RDS"))
  }
  mu <- samps %>% select(contains("mu"))
  # }
  if (i != 1){
    mu_stl <- samps_stl %>% select(contains("mu"))
    sigma_stl <- samps_stl %>% select(sigma)
    if (phenotype != "mds" & phenotype != "mpt"){
      gamma <- samps %>% select(gamma)
    }
  }

  if (phenotype == "mav" | phenotype == "mcv"){
    mu_f <- exp(mu[,1:T] + as.data.frame(rep(gamma,18)) * mu_stl[,1:T] + as.data.frame(rep((gamma*sigma_stl)^2,18))/2)
    mu_m <- exp(mu[,(T+1):(2*T)] + as.data.frame(rep(gamma,18)) * mu_stl[,(T+1):(2*T)] + as.data.frame(rep((gamma*sigma_stl)^2,18))/2)
  } else if (substr(phenotype,1,1) == "m"){
    mu_f <- exp(mu[,1:T])
    mu_m <- exp(mu[,(T+1):(2*T)])
  } else if (i  %in% c(2:10)){
    mu_f <- mu[,1:T] + as.data.frame(rep(gamma,18)) * mu_stl[,1:T]
    mu_m <- mu[,(T+1):(2*T)] + as.data.frame(rep(gamma,18)) * mu_stl[,(T+1):(2*T)]
  } else if (phenotype == "stl"){
    mu_f <- mu[,1:T]
    mu_m <- mu[,(T+1):(2*T)]
  }

  mu_diff_med <- apply(mu_m - mu_f,2,median)
  mu_diff_low <- apply(mu_m - mu_f,2,quantile,0.025)
  mu_diff_upp <- apply(mu_m - mu_f,2,quantile,0.975)

  curr_plot <- ggplot() +
    geom_point(aes(x=years,mu_diff_med)) +
    geom_segment(aes(x=years,xend=years,y=mu_diff_low,yend=mu_diff_upp)) +
    ylab(expression(paste(mu[mt],"-",mu[ft]))) +
    xlab(expression(t)) +
    ggtitle(paste(phenotype)) +
    geom_hline(aes(yintercept=0),colour = "red",linetype=2) +
    theme_bw(base_size = 12)

  if (phenotype == "pmx"){
    dat_pmx_curr <- data.frame(years=years,
                               mu_diff_med=mu_diff_med,
                               mu_diff_low=mu_diff_low,
                               mu_diff_upp=mu_diff_upp)
    fig_4ai <- ggplot(dat_pmx_curr) +
      geom_point(aes(x=years,mu_diff_med)) +
      geom_segment(aes(x=years,xend=years,y=mu_diff_low,yend=mu_diff_upp)) +
      ylab(expression(paste(mu[mt],"-",mu[ft]))) +
      xlab("Relative Time of Deposition (1000s of years)") +
      ggtitle("Premaxilla length (A)") +
      geom_hline(aes(yintercept=0),colour = "red",linetype=2) +
      theme_bw(base_size = 12)
  } else if (phenotype == "tpg"){
    dat_tpg_curr <- data.frame(years=years,
                               mu_diff_med=mu_diff_med,
                               mu_diff_low=mu_diff_low,
                               mu_diff_upp=mu_diff_upp)
    fig_4bi <- ggplot(dat_tpg_curr) +
      geom_point(aes(x=years,mu_diff_med)) +
      geom_segment(aes(x=years,xend=years,y=mu_diff_low,yend=mu_diff_upp)) +
      ylab(expression(paste(mu[mt],"-",mu[ft]))) +
      xlab("Relative Time of Deposition (1000s of years)") +
      ggtitle("Pelvic girdle length (C)") +
      geom_hline(aes(yintercept=0),colour = "red",linetype=2) +
      theme_bw(base_size = 12)
  } else if (phenotype == "ds2"){
    dat_ds2_curr <- data.frame(years=years,
                               mu_diff_med=mu_diff_med,
                               mu_diff_low=mu_diff_low,
                               mu_diff_upp=mu_diff_upp)
    fig_4ci <- ggplot(dat_ds2_curr) +
      geom_point(aes(x=years,mu_diff_med)) +
      geom_segment(aes(x=years,xend=years,y=mu_diff_low,yend=mu_diff_upp)) +
      ylab(expression(paste(mu[mt],"-",mu[ft]))) +
      xlab("Relative Time of Deposition (1000s of years)") +
      ggtitle("Dorsal spine 2 length (E)") +
      geom_hline(aes(yintercept=0),colour = "red",linetype=2) +
      theme_bw(base_size = 12)
  } else if (phenotype == "mcv"){
    dat_mcv_curr <- data.frame(years=years,
                               mu_diff_med=mu_diff_med,
                               mu_diff_low=mu_diff_low,
                               mu_diff_upp=mu_diff_upp)
    fig_4di <- ggplot(dat_mcv_curr) +
      geom_point(aes(x=years,mu_diff_med)) +
      geom_segment(aes(x=years,xend=years,y=mu_diff_low,yend=mu_diff_upp)) +
      ylab(expression(paste(mu[mt],"-",mu[ft]))) +
      xlab("Relative Time of Deposition (1000s of years)") +
      ggtitle("Caudal vert. count (G)") +
      geom_hline(aes(yintercept=0),colour = "red",linetype=2) +
      theme_bw(base_size = 12)
  }

  means <- c(means,mean(apply(mu_m - mu_f,1,mean)))
  sds <- c(sds,var(apply(mu_m - mu_f,1,mean)) + mean(apply(mu_m - mu_f,1,var)))
  tau_inter <- sapply(1:10000,function(r){cor(years,unlist(mu_m[r,] - mu_f[r,]),method="kendall")})
  taus <- c(taus,mean(tau_inter))
  tau_prob <- c(tau_prob,length(which(tau_inter > 0))/10000)
  p_tau <- mean(tau_inter > 0)
    tau_plot <- ggplot() +
      geom_histogram(aes(x=tau_inter)) +
      xlab(expression(tau)) + ylab("") + xlim(-1,1) +
      ggtitle(phenotype) +
      geom_vline(aes(xintercept=mean(tau_inter))) +
      geom_vline(aes(xintercept=0),colour="red",linetype=2) +
      geom_text(
        aes(
          x = -1,
          y = Inf,
          label = sprintf("P(tau > 0) == %.4f", p_tau)
        ),
        parse = TRUE,
        hjust = 0,
        vjust = 1.5
      ) +
      theme_bw(base_size = 12)

    if (phenotype == "pmx"){
      dat_pmx_tau <-  data.frame(tau_inter=tau_inter)
      p_tau_pmx <- p_tau
      fig_4aii <- ggplot(dat_pmx_tau) +
        geom_histogram(aes(x=tau_inter)) +
        xlab(expression(tau)) + ylab("") + xlim(-1,1) +
        ggtitle("Premaxilla length (B)") +
        geom_vline(aes(xintercept=mean(tau_inter))) +
        geom_vline(aes(xintercept=0),colour="red",linetype=2)+
        annotate("text",
                 x = -1,
                 y = Inf,
                 label = sprintf("P(tau > 0) == %.4f", p_tau_pmx),
                 parse = TRUE,
                 hjust = 0,
                 vjust = 1.5) +
        theme_bw(base_size = 12)
    } else if (phenotype == "tpg"){
      dat_tpg_tau <-  data.frame(tau_inter=tau_inter)
      p_tau_tpg <- p_tau
      fig_4bii <- ggplot(dat_tpg_tau) +
        geom_histogram(aes(x=tau_inter)) +
        xlab(expression(tau)) + ylab("") + xlim(-1,1) +
        ggtitle("Pelvic girdle length (D)") +
        geom_vline(aes(xintercept=mean(tau_inter))) +
        geom_vline(aes(xintercept=0),colour="red",linetype=2)+
        annotate("text",
                 x = -1,
                 y = Inf,
                 label = sprintf("P(tau > 0) == %.4f", p_tau_tpg),
                 parse = TRUE,
                 hjust = 0,
                 vjust = 1.5) + theme_bw(base_size = 12)
    } else if (phenotype == "ds2"){
      dat_ds2_tau <-  data.frame(tau_inter=tau_inter)
      p_tau_ds2 <- p_tau
      fig_4cii <- ggplot(dat_ds2_tau) +
        geom_histogram(aes(x=tau_inter)) +
        xlab(expression(tau)) + ylab("") + xlim(-1,1) +
        ggtitle("Dorsal spine 2 length (F)") +
        geom_vline(aes(xintercept=mean(tau_inter))) +
        geom_vline(aes(xintercept=0),colour="red",linetype=2)+
        annotate("text",
                 x = -1,
                 y = Inf,
                 label = sprintf("P(tau > 0) == %.4f", p_tau_ds2),
                 parse = TRUE,
                 hjust = 0,
                 vjust = 1.5) + theme_bw(base_size = 12)
    } else if (phenotype == "mcv"){
      dat_mcv_tau <-  data.frame(tau_inter=tau_inter)
      p_tau_mcv <- p_tau
      fig_4dii <- ggplot(dat_mcv_tau) +
        geom_histogram(aes(x=tau_inter)) +
        xlab(expression(tau)) + ylab("") + xlim(-1,1) +
        ggtitle("Caudal vert. count (H)") +
        geom_vline(aes(xintercept=mean(tau_inter))) +
        geom_vline(aes(xintercept=0),colour="red",linetype=2)+
        annotate("text",
                 x = -1,
                 y = Inf,
                 label = sprintf("P(tau > 0) == %.4f", p_tau_mcv),
                 parse = TRUE,
                 hjust = 0,
                 vjust = 1.5) + theme_bw(base_size = 12)
    }

    gridExtra::grid.arrange(curr_plot,tau_plot,nrow=2) %>%
    ggsave(paste0("../Figures/",phenotype,"_supplemental.pdf"),
             .,
             width=8.5,height=11)
}
gridExtra::grid.arrange(fig_4ai,fig_4aii,fig_4bi,fig_4bii,fig_4ci,fig_4cii,fig_4di,fig_4dii,nrow=4) %>%
  ggsave(filename="../Figures/Figure_4.pdf",width=11,height=14)

for (i in 1:length(cont_w_zeros)){
  phenotype = cont_w_zeros[i]

  samps <- readRDS(paste0("../Figures/posterior_samples_",phenotype,"prob_0.RDS"))
  samps_stl <- readRDS(paste0("../Figures/posterior_samples_stl_0.RDS"))

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
  p_diff_low <- apply(p_m - p_f,1,quantile,0.025)
  p_diff_upp <- apply(p_m - p_f,1,quantile,0.975)

  post_probs_zero <- rbind(post_probs_zero,
                      data.frame(Phenotype = rep(phenotype,18),
                                 Time = years,
                                 Value = apply(p_m - p_f,1,function(x){length(which(x > 0))/length(x)})))


  curr_plot <- ggplot() +
    geom_point(aes(x=years,p_diff_med)) +
    geom_segment(aes(x=years,xend=years,y=p_diff_low,yend=p_diff_upp)) +
    ylab(expression(paste(p[mt],"-",p[ft]))) +
    xlab(expression(t)) +
    ggtitle(paste(phenotype)) +
    geom_hline(aes(yintercept=0),colour = "red",linetype=2) +
    theme_bw()

  p_m_diff_med <- apply(p_m,1,median)
  p_m_diff_low <- apply(p_m,1,quantile,0.025)
  p_m_diff_upp <- apply(p_m,1,quantile,0.975)

  p_f_diff_med <- apply(p_f,1,median)
  p_f_diff_low <- apply(p_f,1,quantile,0.025)
  p_f_diff_upp <- apply(p_f,1,quantile,0.975)

  p_dat <- data.frame(years = rep(years,2),
                      sex = factor(rep(c("Males","Females"),each=T),
                                   levels = c("Males","Females")),
                      p_med = c(p_m_diff_med,p_f_diff_med),
                      p_low = c(p_m_diff_low,p_f_diff_low),
                      p_upp = c(p_m_diff_upp,p_f_diff_upp))
  curr_sep_plot <- ggplot(p_dat) +
    geom_point(aes(x=years,p_med)) +
    geom_segment(aes(x=years,xend=years,y=p_low,yend=p_upp)) +
    facet_wrap(~sex) +
    ylab(expression(p)) +
    xlab(expression(t)) +
    ggtitle(paste(phenotype)) +
    theme_bw()

  curr_sep_plot %>%
    ggsave(paste0("../Figures/",phenotype,"_supplemental_zero_mvf.pdf"),
           .,
           width=8.5,height=11)

  means_zero <- c(means_zero,mean(apply(p_m - p_f,1,mean)))
  sds_zero <- c(sds_zero,var(apply(p_m - p_f,1,mean)) + mean(apply(p_m - p_f,1,var)))
  tau_inter <- sapply(1:10000,function(r){cor(years,unlist(p_m[,r] - p_f[,r]),method="kendall")})
  taus_zero <- c(taus_zero,mean(tau_inter))
  tau_prob_zero <- c(tau_prob_zero,length(which(tau_inter > 0))/10000)
  p_tau <- mean(tau_inter > 0)
  tau_plot <- ggplot() +
    geom_histogram(aes(x=tau_inter)) +
    xlab(expression(tau)) + ylab("") + xlim(-1,1) +
    ggtitle(phenotype) +
    geom_vline(aes(xintercept=mean(tau_inter))) +
    geom_vline(aes(xintercept=0),colour="red",linetype=2)+
    geom_text(
      aes(
        x = -1,
        y = Inf,
        label = sprintf("P(tau > 0) == %.4f", p_tau)
      ),
      parse = TRUE,
      hjust = 0,
      vjust = 1.5
    ) + theme_bw()

  gridExtra::grid.arrange(curr_plot,tau_plot,nrow=2) %>%
    ggsave(paste0("../Figures/",phenotype,"_supplemental_zero.pdf"),
           .,
           width=8.5,height=11)


}
