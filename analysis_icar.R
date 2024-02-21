library(brms)
library(posterior)
person_wise_adj <- as_adj(sim_graph)

sim_patients <- sim_patients %>%
  mutate(infected_bin = ifelse(infected,1,0))

sim_patients_obs <- sim_patients %>%
  filter(detected)

#simple random effect for ward
fit_popn <- brm(infected_bin ~ (1|ward),
           data = sim_patients,
           family = bernoulli())

fit_obs <- brm(infected_bin ~ (1|ward),
                data = sim_patients_obs,
                family = bernoulli())

newdata = data.frame(ward = factor(1:10))

ward_preds <- posterior_linpred(fit_popn, newdata = newdata,
                                transform = TRUE)


obs_ward <- sim_patients_obs %>%
  group_by(ward)%>%
  summarise(observed = TRUE,
            num_obs = n())%>%
  ungroup()%>%
  mutate(ward = factor(paste0("x",ward)))

true_fulldata <- sim_patients %>%
  group_by(ward)%>%
  summarise(popn_truth = mean(infected))%>%
  ungroup()%>%
  mutate(ward = factor(paste0("x",ward))) %>%
  left_join(obs_ward)%>%
  mutate(observed = ifelse(is.na(observed),FALSE,observed))


ward_ests_full <- janitor::clean_names(as_draws_df(ward_preds)) %>%
  pivot_longer(x1:x10, names_to = "ward")%>%
  group_by(ward)%>%
  summarise(estimate = median(value),
            low = quantile(value, .025),
            up = quantile(value, .975)) %>%
  ungroup()%>%
  mutate(data = "full",
         model = "re")%>%
  full_join(true_fulldata)

ward_preds_obs <- posterior_linpred(fit_obs, newdata = newdata,
                                    allow_new_levels = TRUE,
                                    sample_new_levels = "gaussian",
                                    transform = TRUE)

ward_ests_obs <- janitor::clean_names(as_draws_df(ward_preds_obs)) %>%
  pivot_longer(x1:x10, names_to = "ward")%>%
  group_by(ward)%>%
  summarise(estimate = median(value),
            low = quantile(value, .025),
            up = quantile(value, .975)) %>%
  ungroup()%>%
  mutate(data = "obs",
         model = "re")%>%
  full_join(true_fulldata)


model_data = rbind(ward_ests_obs, 
                   ward_ests_full)

ggplot(model_data, aes(x = popn_truth, y = estimate, colour = observed,
                       ymin = low, ymax = up))+
  geom_point()+
  geom_abline()+
  geom_errorbar()+
  facet_grid(.~data)+
  ylab("Ward estimate")+
  xlab("Ward truth")
