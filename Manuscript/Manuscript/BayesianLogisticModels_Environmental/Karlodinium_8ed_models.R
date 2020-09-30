mydataset <- asv.count.HAB10.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  dplyr::rename("presence" = "x") %>%
  filter(pH > 7.5) %>%
  filter(Taxon == "Karlodinium_8ed") %>%
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity))


arm.fit1 <- stan_glmer(presence ~ (1 + pHStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit2 <- stan_glmer(presence ~ (1 + TempStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit3 <- stan_glmer(presence ~ (1 + SalinityStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit4 <- stan_glmer(presence ~ TempStd + (1 + pHStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit5 <- stan_glmer(presence ~ pHStd + (1 + TempStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit6 <- stan_glmer(presence ~ pHStd + SalinityStd + (1 + TempStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit7 <- stan_glm(presence ~ pHStd + TempStd,
                     data = mydataset,
                     family = "binomial",
                     prior_intercept = normal(0, 1),
                     prior = normal(0,1),
                     iter = 1000,
                     chains = 4)
arm.fit8 <- stan_glm(presence ~ pHStd * TempStd,
                     data = mydataset,
                     family = "binomial",
                     prior_intercept = normal(0, 1),
                     prior = normal(0,1),
                     iter = 1000,
                     chains = 4)
arm.fit9 <- stan_glm(presence ~ pHStd + SalinityStd,
                     data = mydataset,
                     family = "binomial",
                     prior_intercept = normal(0, 1),
                     prior = normal(0,1),
                     iter = 1000,
                     chains = 4)
arm.fit10 <- stan_glmer(presence ~ pHStd + SalinityStd + (1 | Season),
                        data = mydataset,
                        family = "binomial",
                        prior_intercept = normal(0, 1),
                        prior = normal(0,1),
                        iter = 1000,
                        chains = 4)
arm.fit11 <- stan_glmer(presence ~ pHStd + (1 + SalinityStd | Season),
                        data = mydataset,
                        family = "binomial",
                        prior_intercept = normal(0, 1),
                        prior = normal(0,1),
                        iter = 1000,
                        chains = 4)
arm.fit12 <- stan_glmer(presence ~ SalinityStd + (1 + pHStd | Season),
                        data = mydataset,
                        family = "binomial",
                        prior_intercept = normal(0, 1),
                        prior = normal(0,1),
                        iter = 1000,
                        chains = 4)
arm.fit13 <- stan_glmer(presence ~ pHStd + (0 + SalinityStd | Season),
                        data = mydataset,
                        family = "binomial",
                        prior_intercept = normal(0, 1),
                        prior = normal(0,1),
                        iter = 1000,
                        chains = 4)
arm.fit14 <- stan_glmer(presence ~ SalinityStd + (0 +  pHStd | Season),
                        data = mydataset,
                        family = "binomial",
                        prior_intercept = normal(0, 1),
                        prior = normal(0,1),
                        iter = 1000,
                        chains = 4)
arm.fit15 <- stan_glmer(presence ~ 0 + (1 + TempStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit16 <- stan_glmer(presence ~ 1 + (0 + TempStd | Season),
                        data = mydataset,
                        family = "binomial",
                        prior_intercept = normal(0, 1),
                        prior = normal(0,1),
                        iter = 1000,
                        chains = 4)

arm.fit17 <- stan_glmer(presence ~ 0 + (0 + TempStd | Season),
                        data = mydataset,
                        family = "binomial",
                        prior_intercept = normal(0, 1),
                        prior = normal(0,1),
                        iter = 1000,
                        chains = 4)
#model compare
modList <- paste0("arm.fit", 1:17)
names(modList) <- paste0("arm.fit", 1:17)
modList %>% 
  map(as.symbol) %>% 
  map(eval) %>% 
  map(waic) %>% 
  loo_compare()

#best model is seasonal, depending only on temperature; look at params
plot(arm.fit17)

mydataset %>% 
  add_fitted_draws(arm.fit17, n = 1000) %>% 
  ggplot(aes(x = Temperature, y = presence, color = Season)) +
  geom_point() +
  facet_grid(~ Season, scales = "free_x") +
  stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
  scale_fill_brewer()

PREDERROR(arm.fit17, mydataset, "presence")

saveRDS(arm.fit17, file = "BayesianLogisticModels_Environmental/Karlodinium_8ed_BestModel.RDS")

