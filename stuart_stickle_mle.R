library(tidyverse)

dat <- read.csv("stuart_stickle_imp_100_raw.csv")
M <- max(dat$imputation_number)
tmp <- dat %>% filter(imputation_number == 1)
T <- max(tmp$time)

new_dat <- dat %>%
  group_by(gender,index) %>%
  tally() %>%
  pivot_wider(names_from="gender",values_from = "n") %>%
  mutate(prop = Female/M) %>%
  select(index,prop) %>%
  right_join(dat %>% filter(imputation_number == 1)) %>%
  select(-c(imputation_number,gender))
