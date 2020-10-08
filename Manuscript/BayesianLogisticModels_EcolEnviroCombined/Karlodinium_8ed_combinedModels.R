tempdataset <- asv.count.all.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  mutate(presence = ifelse(nReads > 0, 1, 0)) %>%
  filter(pH > 7.5) %>%
  filter(Taxon %in% c("Karlodinium_8ed", "Ditylum_a31", "Thalassiosira_6b3", "Poteriospumella_b57")) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity)) %>% 
  pivot_wider(id_cols = c(TempStd, pHStd, SalinityStd, Season),
              names_from = Taxon,
              values_from = presence)


#MODEIFIED model from previous attempts
Karlodinium8ed_envir_model <- stan_glmer(Karlodinium_8ed ~  0 + (0 + TempStd | Season),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0, 1),
                                         iter = 1000,
                                         chains = 4)

Karlodinium8ed_combined_model1 <- stan_glmer(Karlodinium_8ed ~  0 + (0 + TempStd | Season) + (1 | Ditylum_a31),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
# Karlodinium8ed_combined_model2 <- stan_glmer(Karlodinium_8ed ~  0 + (0 + TempStd | Season) + (1 | Ditylum_a31) + (1 | Thalassiosira_6b3),
#                                              data = tempdataset,
#                                              family = "binomial",
#                                              prior_intercept = normal(0, 1),
#                                              prior = normal(0,1),
#                                              iter = 1000,
#                                              chains = 4)
# Karlodinium8ed_combined_model3 <- stan_glmer(Karlodinium_8ed ~  0 + (0 + TempStd | Season) + (1 | Ditylum_a31) + (1 | Thalassiosira_6b3) + (1 | Poteriospumella_b57),
#                                              data = tempdataset,
#                                              family = "binomial",
#                                              prior_intercept = normal(0, 1),
#                                              prior = normal(0,1),
#                                              iter = 1000,
#                                              chains = 4)
# Karlodinium8ed_ecol <- stan_glmer(Karlodinium_8ed ~ (1 | Ditylum_a31) + (1 | Thalassiosira_6b3) + (1 | Poteriospumella_b57),
#                                              data = tempdataset,
#                                              family = "binomial",
#                                              prior_intercept = normal(0, 1),
#                                              prior = normal(0,1),
#                                              iter = 1000,
#                                              chains = 4)
# Karlodinium8ed_min <- stan_glmer(Karlodinium_8ed ~  (1 | Ditylum_a31) + (1 | Thalassiosira_6b3),
#                                 data = tempdataset,
#                                 family = "binomial",
#                                 prior_intercept = normal(0, 1),
#                                 prior = normal(0,1),
#                                 iter = 1000,
#                                 chains = 4)

#adding these observations improves the model a lot!
loo_compare(waic(Karlodinium8ed_envir_model),
            waic(Karlodinium8ed_combined_model1)
            # waic(Karlodinium8ed_combined_model2),
            # waic(Karlodinium8ed_combined_model3),
            # waic(Karlodinium8ed_ecol),
            # waic(Karlodinium8ed_min)
            )  #just other species associations (w no enviro data) is actually a better model than most of the others.


#adding biology improves the pres/abs predictions
PREDERROR(Karlodinium8ed_envir_model, tempdataset, outcomename = "Karlodinium_8ed")
# PREDERROR(Karlodinium8ed_combined_model1, tempdataset, outcomename = "Karlodinium_8ed")
#PREDERROR(Karlodinium8ed_combined_model2, tempdataset, outcomename = "Karlodinium_8ed")
PREDERROR(Karlodinium8ed_combined_model1, tempdataset, outcomename = "Karlodinium_8ed")
# PREDERROR(Karlodinium8ed_ecol, tempdataset, outcomename = "Karlodinium_8ed")

#plot(Karlodinium8ed_combined_model3)  


# xtabs(~ Karlodinium_8ed + Thalassiosira_6b3, data = tempdataset)
# 
# 
# tempdataset %>% 
#   add_fitted_draws(Karlodinium8ed_combined_model3, n = 1000) %>% 
#   ggplot(aes(x = TempStd, y = Karlodinium_8ed)) +
#   geom_point() +
#   facet_grid(Season ~ Ditylum_a31, scales = "free_x") +
#   stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
#   scale_fill_brewer()

saveRDS(Karlodinium8ed_envir_model, "../BayesianLogisticModels_EcolEnviroCombined/Karlodinium_8ed_envir_model.RDS")
saveRDS(Karlodinium8ed_combined_model1, "../BayesianLogisticModels_EcolEnviroCombined/Karlodinium_8ed_combined_model.RDS")
saveRDS(tempdataset, "../BayesianLogisticModels_EcolEnviroCombined/Karlodinium_8ed_tempdataset.RDS")

