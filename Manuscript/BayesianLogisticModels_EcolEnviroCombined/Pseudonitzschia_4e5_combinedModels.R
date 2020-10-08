tempdataset <- asv.count.all.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  mutate(presence = ifelse(nReads > 0, 1, 0)) %>%
  filter(pH > 7.5) %>%
  filter(Taxon %in% c("Pseudonitzschia_4e5", "Calanus_79b", "Pseudonitzschia_d36", "Rhizosolenia_e51")) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity)) %>% 
  pivot_wider(id_cols = c(TempStd, pHStd, SalinityStd, Season),
              names_from = Taxon,
              values_from = presence)


Pseudonitzschia_4e5_envir_model <- stan_glm(Pseudonitzschia_4e5 ~ TempStd + SalinityStd + pHStd,
           data = tempdataset,
           family = "binomial",
           prior_intercept = normal(0, 2),
           prior = normal(0, 2),
           iter = 1000,
           chains = 4)

Pseudonitzschia_4e5_combined_model1 <- stan_glmer(Pseudonitzschia_4e5 ~ TempStd + SalinityStd + pHStd + (1 | Calanus_79b),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 2),
                                         prior = normal(0, 2),
                                         iter = 1000,
                                         chains = 4)

Pseudonitzschia_4e5_combined_model2 <- stan_glmer(Pseudonitzschia_4e5 ~ 0 + TempStd + SalinityStd + pHStd + (1 | Calanus_79b),
                                                  data = tempdataset,
                                                  family = "binomial",
                                                  prior_intercept = normal(0, 2),
                                                  prior = normal(0, 2),
                                                  iter = 1000,
                                                  chains = 4)

#adding these observations improves the model a lot!
loo_compare(waic(Pseudonitzschia_4e5_envir_model),
            waic(Pseudonitzschia_4e5_combined_model1),
            waic(Pseudonitzschia_4e5_combined_model2)
            )

#adding biology improves the pres/abs predictions
PREDERROR(Pseudonitzschia_4e5_envir_model, tempdataset, outcomename = "Pseudonitzschia_4e5")
PREDERROR(Pseudonitzschia_4e5_combined_model2, tempdataset, outcomename = "Pseudonitzschia_4e5")

saveRDS(Pseudonitzschia_4e5_envir_model, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_4e5_envir_model.RDS")
saveRDS(Pseudonitzschia_4e5_combined_model2, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_4e5_combined_model.RDS")
saveRDS(tempdataset, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_4e5_tempdataset.RDS")

