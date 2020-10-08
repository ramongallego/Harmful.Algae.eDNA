#changed the initial readin to match other model testing files, now ONLY d36! needs to be run again to generate BestModel without including 4e5.
mydataset <- asv.count.HAB10.enviro %>%
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  dplyr::rename("presence" = "x") %>%
  filter(pH > 7.5) %>%
  filter(Taxon == "Pseudonitzschia_d36") %>%
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity)) 

mydataset %>% 
  ggplot(aes(x = SalinityStd, y = presence)) +
    geom_point() +
    facet_grid(~Season, scales = "free_x")


arm.fit1 <- stan_glmer(presence ~ pHStd + TempStd + SalinityStd + (1 | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0, 2),
                       iter = 1000,
                       chains = 4)
arm.fit2 <- stan_glmer(presence ~ pHStd + SalinityStd + (1 | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0, 2),
                       iter = 1000,
                       chains = 4)
arm.fit3 <- stan_glmer(presence ~ 0 + pHStd + SalinityStd + (1 | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0, 2),
                       iter = 1000,
                       chains = 4)
arm.fit4 <- stan_glmer(presence ~ 0 + SalinityStd + (1 | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0, 2),
                       iter = 1000,
                       chains = 4)
arm.fit5 <- stan_glm(presence ~ TempStd + SalinityStd + pHStd,
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0, 2),
                       iter = 1000,
                       chains = 4)
arm.fit6 <- stan_glmer(presence ~ 0 + (1 + SalinityStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0, 2),
                       iter = 1000,
                       chains = 4)
arm.fit7 <- stan_glmer(presence ~ 0 + TempStd + (1 + SalinityStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0, 2),
                       iter = 1000,
                       chains = 4)

arm.fit8 <- stan_glmer(presence ~ 0 + TempStd + pHStd + (1 + SalinityStd | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0, 2),
                       iter = 1000,
                       chains = 4)

loo_compare(waic(arm.fit1),
            waic(arm.fit2),
            waic(arm.fit3),
            waic(arm.fit4),
            waic(arm.fit5),
            waic(arm.fit6),
            waic(arm.fit7),
            waic(arm.fit8))

PREDERROR(arm.fit3, mydataset, "presence")
saveRDS(arm.fit3, file = "BayesianLogisticModels_Environmental/Pseudonitzschia_d36_BestModel.RDS")
