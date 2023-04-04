library(tidyverse)
library(mixgb)
library(mice)
library(janitor)
library(missMethods)
library(ROCR)
library(shiny)
library(missForest)

stickle <- read_csv('stickle_revamped.csv')

stick <- stickle %>%
  mutate(gender = as.factor(case_when(
    grepl(".M.",univID,fixed=TRUE) ~ "Male",
    grepl(".F.",univID,fixed=TRUE) ~ "Female",
    grepl("Fos",univID,fixed=TRUE) ~ NA_character_
  )), ips = as.factor(ips)
  ) %>%
  select(gender, everything())

stick_whole <- stick %>% select(-c(univID, ips))

extant_stickle <- stick %>% filter(gender != "Fossil") %>% select(!univID)
extant_stickle <- extant_stickle %>% select(!ips)

univ <- stick %>%
  mutate(time = as.integer(sub(".*([0-9]{2})\\..*", "\\1",univID)))

times <- univ %>% drop_na(time) %>% pull(time)



times_vector <- c(rep(1,43)) %>%
  append(21-times)


imputed_stick_whole <- mice(stick_whole,m=5)

completed_stick_lists <- list()
M=5

for(m in 1:M){
  completed_stick_lists[[m]] <- complete(imputed_stick_whole, m)
}



imputed_means <- function(M,target_variable, completed_stick_lists){
  target_variable <- enquo(target_variable)

  temp_comp_lists <- list()
  est_list <- list()
  for(m in 1:M){

    temp_comp_lists[[m]] <- completed_stick_lists[[m]][-c(1:367),]
    temp_comp_lists[[m]]$time <- times_vector


    est_list[[m]] <- temp_comp_lists[[m]] %>% group_by(gender, time) %>% summarize(xbar = mean(!! target_variable), se = sqrt(var(!! target_variable)/n()), n())
    est_list[[m]]$imp = m
  }


  # stacking the datasets
  combined_est <- do.call(rbind,est_list)
  return(combined_est)
}


combined_se_function <- function(M, variable, completed_stick_lists) {
  variable <- enquo(variable)
  combined_estimate <- imputed_means(5, !!variable, completed_stick_lists)

  combined_est_errors <- combined_estimate %>% group_by(gender, time) %>% summarize(x_bar_bar  = mean(xbar), B_m = var(xbar), W_m = mean(se^2), se_xbar_bar = sqrt(W_m+(1+1/5)*(B_m)), r = (1+1/5)*B_m/W_m, U = (5-1)*(1+(1/r)^2), t_val = qt(0.05/2, df = U, lower.tail = FALSE), uci = x_bar_bar + t_val * se_xbar_bar, lci = x_bar_bar - t_val * se_xbar_bar, t_se = t_val * se_xbar_bar)

  return(combined_est_errors)
}

stickle_imputation_plotting <- function(M, var, var2, completed_stick_lists) {
  var = enquo(var)
  combined_errors <- combined_se_function(5, !!var, completed_stick_lists)

  bb <- ggplot(data= combined_errors, aes(x = time, y = x_bar_bar, col = gender)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = lci, ymax = uci, group = gender),
                                                                                                                               width = 0.2) + labs(x = "Time", y = paste('Mean of imputed', var2))

  return(bb)
}



shinyServer(function(input, output) {

  imputed_stick_whole_mice <- reactive({
    mice(stick_whole, m = 5)
  })

  imputed_stick_whole_mixgb <- reactive({
    mixgb(stick_whole, m = 5, nrounds = 100)
  })

  completed_stick_lists_mice <- reactive({
    lapply(seq_len(5), function(m) complete(imputed_stick_whole_mice(), m))
  })

  completed_stick_lists_mixgb <- reactive({
    imputed_stick_whole_mixgb()
  })

  # Create the plot for MICE imputation
  output$plot_mice <- renderPlot({
    stickle_imputation_plotting(5, get(input$variable), input$variable, completed_stick_lists_mice())
  })

  # Create the plot for mixgb imputation
  output$plot_mixgb <- renderPlot({
    stickle_imputation_plotting(5, get(input$variable), input$variable, completed_stick_lists_mixgb())
  })

  output$accuracy_plot <- renderPlot({
    M = 5
    k = 5
    mice(extant_stickle, method = 'rf', rfpackage = "randomForest")

    stick_matrix <- matrix(NA, nrow=nrow(extant_stickle), ncol = M)
    set.seed(2023)
    fold <- sample(1:k, size = nrow(extant_stickle), replace=TRUE)
    for(i in 1:k){
      extant_stickle_temp <- extant_stickle
      extant_stickle_temp$gender[fold==i] <- NA
      extant_stickle_temp_df <- as.data.frame(extant_stickle_temp)
      doParallel::registerDoParallel(cores = 8) # set based on number of CPU cores
      doRNG::registerDoRNG(seed = 2023)

      imputed_stick <- missForest(extant_stickle_temp_df, parallelize = 'variables')$ximp



      for(n in 1:M){

        imputed_stick$gender[fold==i]
        stick_matrix[fold==i, n] <- as.character(imputed_stick$gender[fold==i])
      }

    }

    pct_male <- apply(stick_matrix,1,function(x){
      mean(x == "Male")
    })

    imptest <- data_frame(truth = as.factor(extant_stickle$gender), pct_male)
    roc_pred <- prediction(imptest$pct_male, imptest$truth)

    auc_rocr <- performance(roc_pred, measure="auc")

    auc_missForest <<- auc_rocr@y.values[[1]]
    plot(performance(roc_pred, measure="tpr", x.measure="fpr"))
  })

  output$accuracy_plot_mixgb <- renderPlot({
    M = 5
    k = 5
    stick_matrix <- matrix(NA, nrow=nrow(extant_stickle), ncol = M)

    set.seed(2023) # seed used for 86% was 12

    fold <- sample(1:k, size = nrow(extant_stickle), replace=TRUE)
    for(i in 1:k){
      extant_stickle_temp <- extant_stickle
      extant_stickle_temp$gender[fold==i] <- NA

      #imputed_stick <- mixgb(extant_stickle_temp,nrounds = 100, m = 5, xgb.params = list(max_depth = 12, eta = 0.014311225124087, gamma = 0,min_child_weight = 1,subsample = 1, colsample_bytree = 1, colsample_bylevel = 1, colsample_bynode = 1,tree_method = "auto", gpu_id = 0, predictor = "auto"))

      imputed_stick <- mixgb(extant_stickle_temp, nrounds = 100, m=5)
      #all that tuning for it to be a worse imputer....
      #whatever


      for(n in 1:M){

        complete_stick <- imputed_stick[[n]]

        complete_stick$gender[fold==i]
        stick_matrix[fold==i, n] <- as.character(complete_stick$gender[fold==i])
      }

    }

    pct_male <- apply(stick_matrix,1,function(x){
      mean(x == "Male")
    })

    imptest <- data_frame(truth = as.factor(extant_stickle$gender), pct_male)
    roc_pred <- prediction(imptest$pct_male, imptest$truth)

    auc_rocr <- performance(roc_pred, measure="auc")
    auc_mixgb <<- auc_rocr@y.values[[1]]
    plot(performance(roc_pred, measure="tpr", x.measure="fpr"))
    })

  output$accuracy_plot_mice <- renderPlot({
    M = 5
    k = 5
    stick_matrix <- matrix(NA, nrow=nrow(extant_stickle), ncol = M)
    set.seed(12)
    fold <- sample(1:k, size = nrow(extant_stickle), replace=TRUE)
    for(i in 1:k){
      extant_stickle_temp <- extant_stickle
      extant_stickle_temp$gender[fold==i] <- NA
      imputed_stick <- mice(extant_stickle_temp,m=M)
      for(n in 1:M){

        complete_stick <- complete(imputed_stick, n)

        complete_stick$gender[fold==i]
        stick_matrix[fold==i, n] <- as.character(complete_stick$gender[fold==i])
      }

    }
    pct_male <- apply(stick_matrix,1,function(x){
      mean(x == "Male")
    })

    imptest <- data_frame(truth = as.factor(extant_stickle$gender), pct_male)
    roc_pred <- prediction(imptest$pct_male, imptest$truth)

    auc_rocr <- performance(roc_pred, measure="auc")

    auc_mice <<- auc_rocr@y.values[[1]]
    plot(performance(roc_pred, measure="tpr", x.measure="fpr"))
  })

output$auc_mice <- renderText({auc_mice})
output$auc_mixgb<- renderText({auc_mixgb})
output$auc_missForest <- renderText({auc_missForest})




})
