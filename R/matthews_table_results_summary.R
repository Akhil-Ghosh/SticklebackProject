library(tidyverse)
library(ROCR)
model_list <- list()
#Indicate which model was used for each trai
model_list[["cle"]] <- 3
model_list[["ds1"]] <- 3
model_list[["ds2"]] <- 1
model_list[["ds3"]] <- 1
model_list[["ect"]] <- 1
model_list[["lps"]] <- 3
model_list[["lpt"]] <- 2
model_list[["maf"]] <- 1
model_list[["mav"]] <- 2
model_list[["mcv"]] <- 0
model_list[["mdf"]] <- 1
#model_list[["mds"]] <- 3
#model_list[["mpt"]] <- 0
model_list[["pmx"]] <- 3
model_list[["stl"]] <- 2
#model_list[["tpg"]] <- 1

library(here)
samps_stl <- readRDS(here("Figures/stl","posterior_samples_stl_2.RDS"))
mu_stl <- samps_stl %>% select(contains("mu"))
maxtime <- 18
results_prob_mu_diff_gt_mu_diff0 <- results_prob_mu_diff_gt_0 <- results_auc_mu_diff <-  data.frame()
for (phenotype in names(model_list)){
  samps <- readRDS(here(paste0("Figures/",phenotype),paste0("posterior_samples_",phenotype,"_",model_list[[phenotype]],".RDS")))
  gamma <- samps %>% select(gamma)
  mu <- samps %>% select(contains("mu"))
  if (phenotype == "stl") {
    mu_f <-
      mu[, 1:maxtime]
    mu_m <-
      mu[, (maxtime + 1):(2 * maxtime)]
  } else {

    if (phenotype %in% c("mdf","maf","mav","mcv")){
      mu_f <-exp(mu[, 1:maxtime])
      mu_m <- exp(mu[, (maxtime + 1):(2 * maxtime)])

    } else {
      mu_f <-
        mu[, 1:maxtime] + as.data.frame(rep(gamma, 18)) * mu_stl[, 1:maxtime]
      mu_m <-
        mu[, (maxtime + 1):(2 * maxtime)] + as.data.frame(rep(gamma, 18)) * mu_stl[, (maxtime +
                                                                                        1):(2 * maxtime)]
    }



    }

  mu_diff <- mu_m - mu_f

  #Probability males are greater than females
  results_prob_mu_diff_gt_0 <- rbind(results_prob_mu_diff_gt_0,data.frame(
    phenotype = phenotype,
    time = 1:18,
    prob_mu_diff_gt_0 = apply(mu_diff, 2, function(x) {
      mean(x > 0)
    })
  ))

  #Differences in the distributions from time point 1 to time point 18
  auc_mu_diff <- c()
  for (j in 1:18){
  pred <- prediction(c(mu_diff[,1], mu_diff[,j]), rep(c(1,0), each = nrow(mu_diff)))
  auc_mu_diff[j] <- perf <- performance(pred, "auc")@y.values[[1]]
  }


  results_auc_mu_diff <- rbind(results_auc_mu_diff,data.frame(
    phenotype = phenotype,
    time = 1:18,
    auc_mu_diff = auc_mu_diff
    ))



  # plot(perf,
  #      avg= "threshold",
  #      colorize=TRUE,
  #      lwd= 3,
  #      main= "With ROCR you can produce standard plots\nlike ROC curves ...")
  # plot(perf,
  #      lty=3,
  #      col="black",
  #      add=TRUE)
  #


}



# ggplot(aes(x = time, y = auc_mu_diff), data = results_auc_mu_diff) + facet_wrap(~phenotype) + geom_path() +  geom_path() + geom_hline(yintercept = 0.5, lty = 3) + ylim(0,1)
# ggplot(aes(x = time, y = prob_mu_diff_gt_mu_diff0), data = results_prob_mu_diff_gt_mu_diff0) + facet_wrap(~phenotype) + geom_path() +  geom_path() + geom_hline(yintercept = 0.5, lty = 3)
#ggplot(aes(x = time, y = prob_mu_diff_gt_0), data = results_prob_mu_diff_gt_0) + facet_wrap(~phenotype) + geom_path() + geom_hline(yintercept = 0.5, lty = 3) + ylim(0,1)


#Now make some tables.
#Rows are traits.




#Which model was used
mod_df <- data.frame(phenotype = names(unlist(model_list)), mod = unlist(model_list))
mod_code <- data.frame(mod = 0:3, mod_name = c("OU-Trend","No OU-Trend",
                                               "OU-No Trend","NO OU-No Trend"))
#Columsn are times.
#Posterior probability that male mean is greater than female mean.
#rows are traits, coliumns are times
#results_prob_mu_diff_gt_0_table <- results_prob_mu_diff_gt_0 %>% mutate(prob_mu_diff_gt_0 = round(prob_mu_diff_gt_0,3)) %>% pivot_wider(names_from = "time", values_from = "prob_mu_diff_gt_0") %>% left_join(mod_df, by = "phenotype") %>% left_join(mod_code, by = "mod") %>%  relocate(mod_name, .after = phenotype)

#cols are traits, rows are times
results_prob_mu_diff_gt_0 %>%  mutate(prob_mu_diff_gt_0 = round(prob_mu_diff_gt_0,3)) %>% pivot_wider(names_from = "phenotype", values_from = "prob_mu_diff_gt_0")
#Model code:
#0: OU-Trend
#1: No OU- Trend
#2: OU - No Trend
#3: No OU - No Trend

