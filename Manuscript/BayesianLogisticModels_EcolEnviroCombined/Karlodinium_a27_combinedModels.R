tempdataset <- asv.count.all.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  mutate(presence = ifelse(nReads > 0, 1, 0)) %>%
  filter(pH > 7.5) %>%
  filter(Taxon %in% c("Karlodinium_a27", "Balanus_2eb", "Hincksia_8a1", "Skeletonema_0a6")) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity)) %>% 
  pivot_wider(id_cols = c(TempStd, pHStd, SalinityStd, Season),
              names_from = Taxon,
              values_from = presence)


#best model from previous attempts
Karlodinium_a27_envir_model <- stan_glmer(Karlodinium_a27 ~ pHStd + (0 + SalinityStd | Season),
                                          data = tempdataset,
                                          family = "binomial",
                                          prior_intercept = normal(0, 1),
                                          prior = normal(0,1),
                                          iter = 1000,
                                          chains = 4)

Karlodinium_a27_combined_model1 <- stan_glmer(Karlodinium_a27 ~ pHStd + (0 + SalinityStd | Season) + (1 | Balanus_2eb),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Karlodinium_a27_combined_model2 <- stan_glmer(Karlodinium_a27 ~ pHStd + (0 + SalinityStd | Season) + (1 | Balanus_2eb) + (1 | Hincksia_8a1),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Karlodinium_a27_combined_model3 <- stan_glmer(Karlodinium_a27 ~ pHStd + (0 + SalinityStd | Season) + (1 | Balanus_2eb) + (1 | Hincksia_8a1) + (1 | Skeletonema_0a6),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Karlodinium_a27_ecol <- stan_glmer(Karlodinium_a27 ~ (1 | Balanus_2eb) + (1 | Hincksia_8a1) + (1 | Skeletonema_0a6),
                              data = tempdataset,
                              family = "binomial",
                              prior_intercept = normal(0, 1),
                              prior = normal(0,1),
                              iter = 1000,
                              chains = 4)


#adding these observations improves the model a lot!
loo_compare(waic(Karlodinium_a27_envir_model),
            waic(Karlodinium_a27_combined_model1),
            waic(Karlodinium_a27_combined_model2),
            waic(Karlodinium_a27_combined_model3),
            waic(Karlodinium_a27_ecol))  #just other species associations (w no enviro data) is actually a better model than most of the others.

#adding biology improves the pres/abs predictions
PREDERROR(Karlodinium_a27_envir_model, tempdataset, outcomename = "Karlodinium_a27")
# PREDERROR(Karlodinium_a27_combined_model1, tempdataset, outcomename = "Karlodinium_a27")
# PREDERROR(Karlodinium_a27_combined_model2, tempdataset, outcomename = "Karlodinium_a27")
PREDERROR(Karlodinium_a27_combined_model3, tempdataset, outcomename = "Karlodinium_a27")
PREDERROR(Karlodinium_a27_ecol, tempdataset, outcomename = "Karlodinium_a27")

#xtabs(~ Karlodinium_a27 + Balanus_2eb, data = tempdataset)

saveRDS(Karlodinium_a27_envir_model, "BayesianLogisticModels_EcolEnviroCombined/Karlodinium_a27_envir_model.RDS")
saveRDS(Karlodinium_a27_combined_model1, "BayesianLogisticModels_EcolEnviroCombined/Karlodinium_a27_combined_model.RDS")
saveRDS(tempdataset, "BayesianLogisticModels_EcolEnviroCombined/Karlodinium_a27_tempdataset.RDS")

