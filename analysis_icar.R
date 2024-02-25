library(brms)
library(posterior)
person_wise_adj <- as_adj(sim_graph)
ward_wise_adj <- ward_adj

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

#icar models
ward_infected <- sim_patients %>%
  group_by(ward)%>%
  summarise(n_pos = sum(infected_bin), num_tested = n())%>%
  ungroup()

ward_infected_obs <- sim_patients_obs %>%
  group_by(ward)%>%
  summarise(n_pos = sum(infected_bin), num_tested = n())%>%
  ungroup() %>%
  rbind(data.frame(ward = c(3, 8),
                   n_pos = c(0,0),
                   num_tested = c(0,0)))

fit_popn_icar <- brm(n_pos|trials(num_tested) ~ car(W, gr = ward, type = "icar") ,
                data = ward_infected, data2 = list(W = ward_wise_adj), 
                family = binomial())

fit_obs_icar <- brm(n_pos|trials(num_tested) ~ car(W, gr = ward, type = "icar"),
               data = ward_infected_obs, data2 = list(W = ward_wise_adj), 
               family = binomial())

### icar with strength 

fit_popn_icar_strength <- brm(n_pos|trials(num_tested) ~ car(W, gr = ward, type = "icar"),
                     data = ward_infected, data2 = list(W = ward_adj_strength), 
                     family = binomial())

fit_obs_icar_strength <- brm(n_pos|trials(num_tested) ~ car(W, gr = ward, type = "icar"),
                    data = ward_infected_obs, data2 = list(W = ward_adj_strength), 
                    family = binomial())


newdata = data.frame(ward = factor(1:45), 
                     num_tested = rep(1,45))


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


ward_preds <- posterior_linpred(fit_popn, newdata = newdata,
                                    allow_new_levels = TRUE,
                                    sample_new_levels = "gaussian",
                                    transform = TRUE)


ward_ests_full <- janitor::clean_names(as_draws_df(ward_preds)) %>%
  pivot_longer(x1:x45, names_to = "ward")%>%
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
  pivot_longer(x1:x45, names_to = "ward")%>%
  group_by(ward)%>%
  summarise(estimate = median(value),
            low = quantile(value, .025),
            up = quantile(value, .975)) %>%
  ungroup()%>%
  mutate(data = "obs",
         model = "re")%>%
  full_join(true_fulldata)


#### icar model predictions
ward_preds_full_icar <- posterior_linpred(fit_popn_icar, newdata = newdata,
                                    allow_new_levels = TRUE,
                                    sample_new_levels = "gaussian",
                                    transform = TRUE)

ward_ests_full_icar <- janitor::clean_names(as_draws_df(ward_preds_full_icar)) %>%
  pivot_longer(x1:x45, names_to = "ward")%>%
  group_by(ward)%>%
  summarise(estimate = median(value),
            low = quantile(value, .025),
            up = quantile(value, .975)) %>%
  ungroup()%>%
  mutate(data = "full",
         model = "icar")%>%
  full_join(true_fulldata)

ward_preds_obs_icar <- posterior_linpred(fit_obs_icar, newdata = newdata,
                                    allow_new_levels = TRUE,
                                    sample_new_levels = "gaussian",
                                    transform = TRUE)

ward_ests_obs_icar <- janitor::clean_names(as_draws_df(ward_preds_obs_icar)) %>%
  pivot_longer(x1:x45, names_to = "ward")%>%
  group_by(ward)%>%
  summarise(estimate = median(value),
            low = quantile(value, .025),
            up = quantile(value, .975)) %>%
  ungroup()%>%
  mutate(data = "obs",
         model = "icar")%>%
  full_join(true_fulldata)

#### icar model predictions with strength of prediction
ward_preds_full_icar_strength <- posterior_linpred(fit_popn_icar_strength, newdata = newdata,
                                          allow_new_levels = TRUE,
                                          sample_new_levels = "gaussian",
                                          transform = TRUE)

ward_ests_full_icar_strength <- janitor::clean_names(as_draws_df(ward_preds_full_icar_strength)) %>%
  pivot_longer(x1:x45, names_to = "ward")%>%
  group_by(ward)%>%
  summarise(estimate = median(value),
            low = quantile(value, .025),
            up = quantile(value, .975)) %>%
  ungroup()%>%
  mutate(data = "full",
         model = "icar - strength")%>%
  full_join(true_fulldata)

ward_preds_obs_icar_strength <- posterior_linpred(fit_obs_icar_strength, newdata = newdata,
                                         allow_new_levels = TRUE,
                                         sample_new_levels = "gaussian",
                                         transform = TRUE)

ward_ests_obs_icar_strength <- janitor::clean_names(as_draws_df(ward_preds_obs_icar_strength)) %>%
  pivot_longer(x1:x45, names_to = "ward")%>%
  group_by(ward)%>%
  summarise(estimate = median(value),
            low = quantile(value, .025),
            up = quantile(value, .975)) %>%
  ungroup()%>%
  mutate(data = "obs",
         model = "icar - strength")%>%
  full_join(true_fulldata)

disagg_obs <- sim_patients_obs %>%
  group_by(ward)%>%
  summarise(estimate = mean(infected))%>%
  ungroup()%>%
  mutate(low = NA, up = NA, data = "obs",
         model = "raw_obs",
         ward = factor(paste0("x",ward)))%>%
  left_join(true_fulldata)

model_data = rbind(ward_ests_obs, 
                   ward_ests_full,
                   ward_ests_full_icar,
                   ward_ests_obs_icar,
                   ward_ests_full_icar_strength,
                   ward_ests_obs_icar_strength,
                   disagg_obs)

ggplot(model_data, aes(x = popn_truth, y = estimate, colour = observed,
                       ymin = low, ymax = up))+
  geom_point()+
  geom_abline()+
  geom_errorbar()+
  facet_grid(model~data)+
  ylab("Ward estimate")+
  xlab("Ward truth")

ggplot(model_data, aes(x = popn_truth, y = estimate, colour = model,
                       ymin = low, ymax = up))+
  geom_point()+
  geom_abline()+
  geom_errorbar()+
  facet_grid(model~data)+
  ylab("Ward estimate")+
  xlab("Ward truth")

