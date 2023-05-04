samps <- readRDS("posterior_samples.RDS")

samps %>%
  rename(Female = theta_f) %>%
  rename(Male = theta_m) %>%
  select(Female,Male) %>%
  gather("Gender","Value") %>%
  ggplot() +
  geom_density(aes(x=Value,color=Gender)) +
  theme_bw() +
  xlim(c(30,80)) +
  xlab(expression(theta)) +
  ylab(expression(paste("p(",theta,"|",x,")")))

length(which(samps$theta_m > samps$theta_f))/nrow(samps)

mu_f <- samps %>% select(contains("mu_f"))
mu_m <- samps %>% select(contains("mu_m"))

mu_diff <- mu_m - mu_f

data.frame(time = 1:18,
           med = apply(mu_diff,2,quantile,probs=0.5),
           low = apply(mu_diff,2,quantile,probs=0.025),
           upp = apply(mu_diff,2,quantile,probs=0.975)) %>%
  ggplot(aes(x=time,y=med)) +
  geom_point() +
  geom_errorbar(aes(ymin=low,ymax=upp),color = "red") +
  theme_bw() +
  xlab(expression(t)) +
  ylab(expression(paste(mu[m],"-",mu[f])))

names(mu_diff) <- c(1:18)
as.data.frame(mu_diff) %>%
  gather("Time","Value") %>%
  mutate(Time = as.numeric(Time)) %>%
  ggplot(aes(x=Value,y=Time,group=Time)) +
  geom_density_ridges(fill="blue",alpha=0.5) +
  xlim(c(-15,20)) +
  scale_y_reverse() +
  geom_vline(aes(xintercept = 0),colour = "red",linetype=2) +
  xlab(expression(paste(mu[mt],"-",mu[ft]))) +
  ylab(expression(t)) +
  theme_bw()

post_probs <- apply(mu_diff,2,function(x){length(which(x > 0))})/nrow(mu_diff)

data.frame(time=1:18,probs=post_probs) %>%
  ggplot() +
  geom_line(aes(x=time,y=probs)) +
  xlab(expression(t)) +
  ylab(expression(paste(P,"(",mu[mt] > mu[ft],"|",X,")"))) +
  theme_bw()


