tempdataset <- asv.count.all.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  mutate(presence = ifelse(nReads > 0, 1, 0)) %>%
  filter(pH > 7.5) %>%
  filter(Taxon %in% c("Pseudonitzschia_d40", "Ditylum_ce4", "Synchaeta_06c", "Calanus_79b")) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity)) %>% 
  pivot_wider(id_cols = c(TempStd, pHStd, SalinityStd, Season),
              names_from = Taxon,
              values_from = presence)


Pseudonitzschia_d40_envir_model <- stan_glm(Pseudonitzschia_d40 ~ SalinityStd,
           data = tempdataset,
           family = "binomial",
           prior_intercept = normal(0, 1),
           prior = normal(0, 1),
           iter = 1000,
           chains = 4)

Pseudonitzschia_d40_envir_model <- stan_glm(Pseudonitzschia_d40 ~ SalinityStd,
                                            data = tempdataset,
                                            family = "binomial",
                                            prior_intercept = normal(0, 1),
                                            prior = normal(0, 1),
                                            iter = 1000,
                                            chains = 4)



Pseudonitzschia_d40_combined_model1 <- stan_glmer(Pseudonitzschia_d40 ~ 0 + SalinityStd + (1 | Ditylum_ce4),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 2),
                                         prior = normal(0, 2),
                                         iter = 1000,
                                         chains = 4)

Pseudonitzschia_d40_combined_model2 <- stan_glmer(Pseudonitzschia_d40 ~ SalinityStd + (1 | Ditylum_ce4),
                                                  data = tempdataset,
                                                  family = "binomial",
                                                  prior_intercept = normal(0, 2),
                                                  prior = normal(0, 2),
                                                  iter = 1000,
                                                  chains = 4)

#adding these observations improves the model a lot!
loo_compare(waic(Pseudonitzschia_d40_envir_model),
            waic(Pseudonitzschia_d40_combined_model1),
            waic(Pseudonitzschia_d40_combined_model2)
            )

#adding biology improves the pres/abs predictions
PREDERROR(Pseudonitzschia_d40_envir_model, tempdataset, outcomename = "Pseudonitzschia_d40")
PREDERROR(Pseudonitzschia_d40_combined_model1, tempdataset, outcomename = "Pseudonitzschia_d40")

saveRDS(Pseudonitzschia_d40_envir_model, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_d40_envir_model.RDS")
saveRDS(Pseudonitzschia_d40_combined_model1, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_d40_combined_model.RDS")
saveRDS(tempdataset, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_d40_tempdataset.RDS")

