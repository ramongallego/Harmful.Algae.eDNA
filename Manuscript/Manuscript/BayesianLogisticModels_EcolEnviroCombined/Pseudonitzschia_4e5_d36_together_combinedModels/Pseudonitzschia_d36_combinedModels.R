#first, combine P-N haplotypes 
tempdataset <- asv.count.HAB10.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  filter(Taxon %in% c("Pseudonitzschia_d36", "Pseudonitzschia_4e5")) %>% 
  group_by(Site, Month, Year) %>% 
  summarise(nReads = sum(nReads),
            Temperature = unique(Temperature),
            Salinity = unique(Salinity),
            pH = unique(pH),
            Season = unique(Season),
            Taxon = "Pseudonitzschia_d36")
  
  #then merge in predictors
tempdataset <- asv.count.all.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  filter(Taxon %in% c("Calanus_79b", "Clytia_8b4")) %>% 
  full_join(tempdataset) %>% 
  group_by(Taxon) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity),
         presence = ifelse(nReads > 0, 1, 0)) %>%
  pivot_wider(id_cols = c(TempStd, pHStd, SalinityStd, Season),
              names_from = Taxon,
              values_from = presence)

Pseudonitzschia_d36_envir_model <- stan_glm(Pseudonitzschia_d36 ~ TempStd + SalinityStd + pHStd,
                                            data = tempdataset,
                                            family = "binomial",
                                            prior_intercept = normal(0, 2),
                                            prior = normal(0, 2),
                                            iter = 1000,
                                            chains = 4)

Pseudonitzschia_d36_combined_model1 <- stan_glmer(Pseudonitzschia_d36 ~ 0 + TempStd + SalinityStd + pHStd + (1 | Calanus_79b),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 2),
                                         prior = normal(0, 2),
                                         iter = 1000,
                                         chains = 4)

Pseudonitzschia_d36_combined_model2 <- stan_glmer(Pseudonitzschia_d36 ~ 0 + TempStd + SalinityStd + (1 | Calanus_79b),
                                                  data = tempdataset,
                                                  family = "binomial",
                                                  prior_intercept = normal(0, 2),
                                                  prior = normal(0, 2),
                                                  iter = 1000,
                                                  chains = 4)

Pseudonitzschia_d36_combined_model3 <- stan_glmer(Pseudonitzschia_d36 ~ 0 + (1 | Season) + (1 | Calanus_79b),
                                                  data = tempdataset,
                                                  family = "binomial",
                                                  prior_intercept = normal(0, 2),
                                                  prior = normal(0, 2),
                                                  iter = 1000,
                                                  chains = 4)


#adding these observations improves the model a lot!
loo_compare(waic(Pseudonitzschia_d36_envir_model),
            waic(Pseudonitzschia_d36_combined_model1),
            waic(Pseudonitzschia_d36_combined_model2),
            waic(Pseudonitzschia_d36_combined_model3)
            )

#adding biology improves the pres/abs predictions
PREDERROR(Pseudonitzschia_d36_envir_model, tempdataset, outcomename = "Pseudonitzschia_d36")
PREDERROR(Pseudonitzschia_d36_combined_model1, tempdataset, outcomename = "Pseudonitzschia_d36")

saveRDS(Pseudonitzschia_d36_envir_model, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_d36_envir_model.RDS")
saveRDS(Pseudonitzschia_d36_combined_model1, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_d36_combined_model.RDS")
saveRDS(tempdataset, "BayesianLogisticModels_EcolEnviroCombined/Pseudonitzschia_d36_tempdataset.RDS")

