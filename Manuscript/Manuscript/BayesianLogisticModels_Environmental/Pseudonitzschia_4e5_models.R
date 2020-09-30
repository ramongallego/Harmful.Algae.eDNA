mydataset <- asv.count.HAB10.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  dplyr::rename("presence" = "x") %>%
  filter(pH > 7.5) %>%
  filter(Taxon == "Pseudonitzschia_4e5") %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity))

mydataset %>% 
  ggplot(aes(x = pHStd, y = presence)) +
    geom_point() +
    facet_grid(~Season)


arm.fit1 <- stan_glmer(presence ~ 0 + TempStd + SalinityStd + (1 | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0,2),
                       iter = 1000,
                       chains = 4)
arm.fit2 <- stan_glmer(presence ~ TempStd + SalinityStd + (1 | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0,2),
                       iter = 1000,
                       chains = 4)
arm.fit3 <- stan_glm(presence ~ TempStd + SalinityStd,
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0,2),
                       iter = 1000,
                       chains = 4)
arm.fit4 <- stan_glm(presence ~ TempStd + SalinityStd + pHStd,
                     data = mydataset,
                     family = "binomial",
                     prior_intercept = normal(0, 2),
                     prior = normal(0,2),
                     iter = 1000,
                     chains = 4)
arm.fit5 <- stan_glmer(presence ~ pHStd + TempStd + SalinityStd + (1 | Season),
                       data = mydataset,
                       family = "binomial",
                       prior_intercept = normal(0, 2),
                       prior = normal(0,2),
                       iter = 1000,
                       chains = 4)


loo_compare(waic(arm.fit1),
            waic(arm.fit2),
            waic(arm.fit3),
            waic(arm.fit4),
            waic(arm.fit5)
            )


PREDERROR(arm.fit4, mydataset, "presence")
saveRDS(arm.fit4, file = "BayesianLogisticModels_Environmental/Pseudonitzschia_4e5_BestModel.RDS")
