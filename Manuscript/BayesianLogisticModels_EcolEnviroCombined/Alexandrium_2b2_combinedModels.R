tempdataset <- asv.count.all.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  mutate(presence = ifelse(nReads > 0, 1, 0)) %>%
  filter(pH > 7.5) %>%
  filter(Taxon %in% c("Alexandrium_2b2", "Ditylum_ba4", "Ditylum_a31", "Thalassiosira_47a")) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity)) %>% 
  pivot_wider(id_cols = c(TempStd, pHStd, SalinityStd, Season),
              names_from = Taxon,
              values_from = presence)


#best model from previous attempts (modified to have no universal intercept)
Alexandrium_2b2_envir_model <- stan_glmer(Alexandrium_2b2 ~ 0 + pHStd + (1 + TempStd | Season),
           data = tempdataset,
           family = "binomial",
           prior_intercept = normal(0, 1),
           prior = normal(0,1),
           iter = 1000,
           chains = 4)



Alexandrium_2b2_combined_model1 <- stan_glmer(Alexandrium_2b2 ~ 0 + pHStd + TempStd + (1 | Ditylum_ba4),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Alexandrium_2b2_combined_model2 <- stan_glmer(Alexandrium_2b2 ~ pHStd + TempStd + (1 | Ditylum_ba4) + (1 | Ditylum_a31),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Alexandrium_2b2_combined_model3 <- stan_glmer(Alexandrium_2b2 ~ pHStd + TempStd + (1 | Ditylum_ba4) + (1 | Ditylum_a31) + (1 | Thalassiosira_47a),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Alexandrium_2b2_combined_model4 <- stan_glmer(Alexandrium_2b2 ~ 0 + pHStd + (1 + TempStd | Season) + (1 | Ditylum_ba4),
                                            data = tempdataset,
                                            family = "binomial",
                                            prior_intercept = normal(0, 1),
                                            prior = normal(0,1),
                                            iter = 1000,
                                            chains = 4)
Alexandrium_2b2_combined_model5 <- stan_glmer(Alexandrium_2b2 ~ 0 + pHStd + (1 + TempStd | Season) + (1 | Ditylum_ba4) + (1 | Ditylum_a31),
                                            data = tempdataset,
                                            family = "binomial",
                                            prior_intercept = normal(0, 1),
                                            prior = normal(0,1),
                                            iter = 1000,
                                            chains = 4)
Alexandrium_2b2_combined_model6 <- stan_glmer(Alexandrium_2b2 ~ 0 + pHStd + (1 + TempStd | Season) + (1 | Ditylum_ba4) + (1 | Ditylum_a31) + (1 | Thalassiosira_47a),
                                            data = tempdataset,
                                            family = "binomial",
                                            prior_intercept = normal(0, 1),
                                            prior = normal(0,1),
                                            iter = 1000,
                                            chains = 4)
Alexandrium_2b2_ecol <- stan_glmer(Alexandrium_2b2 ~ (1 | Ditylum_ba4) + (1 | Ditylum_a31) + (1 | Thalassiosira_47a),
                              data = tempdataset,
                              family = "binomial",
                              prior_intercept = normal(0, 1),
                              prior = normal(0,1),
                              iter = 1000,
                              chains = 4)


#adding these observations improves the model a lot!
loo_compare(waic(Alexandrium_2b2_envir_model),
            waic(Alexandrium_2b2_combined_model1),
            waic(Alexandrium_2b2_combined_model2),
            waic(Alexandrium_2b2_combined_model3),
            waic(Alexandrium_2b2_combined_model4),
            waic(Alexandrium_2b2_combined_model5),
            waic(Alexandrium_2b2_combined_model6),
            waic(Alexandrium_2b2_ecol))  #just other species associations (w no enviro data) is actually a better model than most of the others.

#adding biology improves the pres/abs predictions
PREDERROR(Alexandrium_2b2_envir_model, tempdataset, outcomename = "Alexandrium_2b2")
PREDERROR(Alexandrium_2b2_combined_model6, tempdataset, outcomename = "Alexandrium_2b2")
# PREDERROR(Alexandrium_2b2_combined_model5, tempdataset, outcomename = "Alexandrium_2b2")
# PREDERROR(Alexandrium_2b2_combined_model3, tempdataset, outcomename = "Alexandrium_2b2")


saveRDS(Alexandrium_2b2_envir_model, "BayesianLogisticModels_EcolEnviroCombined/Alexandrium_2b2_envir_model.RDS")
saveRDS(Alexandrium_2b2_combined_model1, "BayesianLogisticModels_EcolEnviroCombined/Alexandrium_2b2_combined_model.RDS")
saveRDS(tempdataset, "BayesianLogisticModels_EcolEnviroCombined/Alexandrium_2b2_tempdataset.RDS")


# tempdataset %>% 
#   add_fitted_draws(Hematodinium_combined_model1, n = 1000) %>% 
#   ggplot(aes(x = pHStd, y = Alexandrium_2b2)) +
#   geom_point() +
#   facet_grid(TempStd > 0 ~ Ditylum_ba4, scales = "free_x") +
#   stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
#   scale_fill_brewer()
