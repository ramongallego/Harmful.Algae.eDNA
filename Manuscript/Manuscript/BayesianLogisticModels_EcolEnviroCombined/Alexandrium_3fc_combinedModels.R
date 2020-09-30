tempdataset <- asv.count.all.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  mutate(presence = ifelse(nReads > 0, 1, 0)) %>%
  filter(pH > 7.5) %>%
  filter(Taxon %in% c("Alexandrium_3fc", "Prasinoderma_6ac", "Stephanodiscus_0c0", "Chaetoceros_dcd")) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity)) %>% 
  pivot_wider(id_cols = c(TempStd, pHStd, SalinityStd, Season),
              names_from = Taxon,
              values_from = presence)


#best model from previous attempts (slightly modified to remove universal intercept, which wasn't doing anything)
Alexandrium_3fc_envir_model <- stan_glmer(Alexandrium_3fc ~ 0 + (1 + TempStd | Season),
                                          data = tempdataset,
                                          family = "binomial",
                                          prior_intercept = normal(0, 1),
                                          prior = normal(0,1),
                                          iter = 1000,
                                          chains = 4)

Alexandrium_3fc_combined_model1 <- stan_glmer(Alexandrium_3fc ~ 0 + (1 + TempStd | Season) + (1 | Prasinoderma_6ac),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Alexandrium_3fc_combined_model2 <- stan_glmer(Alexandrium_3fc ~ 0 + (1 + TempStd | Season) + (1 | Prasinoderma_6ac) + (1 | Stephanodiscus_0c0),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Alexandrium_3fc_combined_model3 <- stan_glmer(Alexandrium_3fc ~ 0 + (1 + TempStd | Season) + (1 | Prasinoderma_6ac) + (1 | Stephanodiscus_0c0) + (1 | Chaetoceros_dcd),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Alexandrium_3fc_ecol <- stan_glmer(Alexandrium_3fc ~ (1 | Prasinoderma_6ac) + (1 | Stephanodiscus_0c0) + (1 | Chaetoceros_dcd),
                              data = tempdataset,
                              family = "binomial",
                              prior_intercept = normal(0, 1),
                              prior = normal(0,1),
                              iter = 1000,
                              chains = 4)


#adding these observations improves the model a lot!
loo_compare(waic(Alexandrium_3fc_envir_model),
            waic(Alexandrium_3fc_combined_model1),
            waic(Alexandrium_3fc_combined_model2),
            waic(Alexandrium_3fc_combined_model3),
            waic(Alexandrium_3fc_ecol))  #just other species associations (w no enviro data) is actually a better model than most of the others.

#adding biology improves the pres/abs predictions
PREDERROR(Alexandrium_3fc_envir_model, tempdataset, outcomename = "Alexandrium_3fc")
# PREDERROR(Alexandrium_3fc_combined_model1, tempdataset, outcomename = "Alexandrium_3fc")
# PREDERROR(Alexandrium_3fc_combined_model2, tempdataset, outcomename = "Alexandrium_3fc")
PREDERROR(Alexandrium_3fc_combined_model3, tempdataset, outcomename = "Alexandrium_3fc")
PREDERROR(Alexandrium_3fc_ecol, tempdataset, outcomename = "Alexandrium_3fc")


saveRDS(Alexandrium_3fc_envir_model, "BayesianLogisticModels_EcolEnviroCombined/Alexandrium_3fc_envir_model.RDS")
saveRDS(Alexandrium_3fc_combined_model1, "BayesianLogisticModels_EcolEnviroCombined/Alexandrium_3fc_combined_model.RDS")
saveRDS(tempdataset, "BayesianLogisticModels_EcolEnviroCombined/Alexandrium_3fc_tempdataset.RDS")




#plot(Hematodinium_combined_model3)  #to see that all 3 of these predictor species are individually helpful, and each better than the existing enviro predictors.

#but it still feels like there's a risk of overfitting here. waic/loo guard against this, but let's figure out how much to trust that. 

# invLogit <- function(x){exp(x)/(1 + exp(x))}
# #at mean temperature and pH, the chance of seeing Hematodinium without Saxidomus is 
# Hematodinium_combined_model1$coefficients[1] %>% invLogit()
# #and, given Saxidomus, the chance of seeing Hematodinium is
# sum(Hematodinium_combined_model1$coefficients[c(1,4)]) %>% invLogit()
# #confirm this modeled conclusion by looking at the crosstabs of observed data
# xtabs(~ Alexandrium_3fc + Prasinoderma_6ac, data = tempdataset)
# 
# 
# tempdataset %>% 
#   add_fitted_draws(Hematodinium_combined_model1, n = 1000) %>% 
#   ggplot(aes(x = pHStd, y = Alexandrium_3fc)) +
#   geom_point() +
#   facet_grid(TempStd > 0 ~ Prasinoderma_6ac, scales = "free_x") +
#   stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
#   scale_fill_brewer()
