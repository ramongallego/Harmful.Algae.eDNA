
#as an exercise, condition our model on the presence of one or more predictor taxa:
tempdataset <- asv.count.all.enviro %>% 
  separate(sample, into = c("Site", "Month", "Year", "z")) %>%
  mutate(Season = ifelse(Month %in% c("10","11","12","1","2","3"), "Winter", "Summer")) %>%
  mutate(presence = ifelse(nReads > 0, 1, 0)) %>%
  filter(pH > 7.5) %>%
  filter(Taxon %in% c("Hematodinium_449", "Saxidomus_33e", "Chrysochromulina_7aa", "Cylindrotheca_799")) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity)) %>% 
  pivot_wider(id_cols = c(TempStd, pHStd, SalinityStd, Season),
              names_from = Taxon,
              values_from = presence)


#MODIFIED model from previous attempts
Hematodinium_449_envir_model <- stan_glm(Hematodinium_449 ~ pHStd + TempStd,
                                       data = tempdataset,
                                       family = "binomial",
                                       prior_intercept = normal(0, 1),
                                       prior = normal(0,1),
                                       iter = 1000,
                                       chains = 4)

Hematodinium_449_combined_model1 <- stan_glmer(Hematodinium_449 ~ pHStd + TempStd + (1 | Saxidomus_33e),
                                         data = tempdataset,
                                         family = "binomial",
                                         prior_intercept = normal(0, 1),
                                         prior = normal(0,1),
                                         iter = 1000,
                                         chains = 4)
Hematodinium_449_combined_model2 <- stan_glmer(Hematodinium_449 ~ pHStd + TempStd + (1 | Saxidomus_33e) + (1 | Chrysochromulina_7aa),
                                             data = tempdataset,
                                             family = "binomial",
                                             prior_intercept = normal(0, 1),
                                             prior = normal(0,1),
                                             iter = 1000,
                                             chains = 4)
Hematodinium_449_combined_model3 <- stan_glmer(Hematodinium_449 ~ pHStd + TempStd + (1 | Saxidomus_33e) + (1 | Chrysochromulina_7aa) + (1 | Cylindrotheca_799),
                                             data = tempdataset,
                                             family = "binomial",
                                             prior_intercept = normal(0, 1),
                                             prior = normal(0,1),
                                             iter = 1000,
                                             chains = 4)
Hematodinium_449_ecol <- stan_glmer(Hematodinium_449 ~ (1 | Saxidomus_33e) + (1 | Chrysochromulina_7aa) + (1 | Cylindrotheca_799),
                                             data = tempdataset,
                                             family = "binomial",
                                             prior_intercept = normal(0, 1),
                                             prior = normal(0,1),
                                             iter = 1000,
                                             chains = 4)
Hematodinium_449_min <- stan_glmer(Hematodinium_449 ~  (1 | Saxidomus_33e) + (1 | Chrysochromulina_7aa),
                                data = tempdataset,
                                family = "binomial",
                                prior_intercept = normal(0, 1),
                                prior = normal(0,1),
                                iter = 1000,
                                chains = 4)

#adding these observations improves the model a lot!
loo_compare(waic(Hematodinium_449_envir_model),
            waic(Hematodinium_449_combined_model1),
            waic(Hematodinium_449_combined_model2),
            waic(Hematodinium_449_combined_model3),
            waic(Hematodinium_449_ecol),
            waic(Hematodinium_449_min))  #just other species associations (w no enviro data) is actually a better model than most of the others.


#adding biology improves the pres/abs predictions
PREDERROR(Hematodinium_449_envir_model, tempdataset, outcomename = "Hematodinium_449")
# PREDERROR(Hematodinium_449_combined_model1, tempdataset, outcomename = "Hematodinium_449")
#PREDERROR(Hematodinium_449_combined_model2, tempdataset, outcomename = "Hematodinium_449")
PREDERROR(Hematodinium_449_combined_model3, tempdataset, outcomename = "Hematodinium_449")
# PREDERROR(Hematodinium_449_ecol, tempdataset, outcomename = "Hematodinium_449")

#plot(Hematodinium_449_combined_model3)  

#but it still feels like there's a risk of overfitting here. waic/loo guard against this, but let's figure out how much to trust that. 

tempdataset %>% 
  add_fitted_draws(Hematodinium_449_combined_model3, n = 1000) %>% 
  ggplot(aes(x = TempStd, y = Hematodinium_449)) +
  geom_point() +
  facet_grid(Season ~ Saxidomus_33e, scales = "free_x") +
  stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
  scale_fill_brewer()

saveRDS(Hematodinium_449_envir_model, "BayesianLogisticModels_EcolEnviroCombined/Hematodinium_449_envir_model.RDS")
saveRDS(Hematodinium_449_combined_model1, "BayesianLogisticModels_EcolEnviroCombined/Hematodinium_449_combined_model.RDS")
saveRDS(tempdataset, "BayesianLogisticModels_EcolEnviroCombined/Hematodinium_449_tempdataset.RDS")



