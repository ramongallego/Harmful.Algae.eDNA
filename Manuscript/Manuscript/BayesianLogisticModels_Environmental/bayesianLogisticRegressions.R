#RPK bayesian models for selected HAB taxa, environmental variables

library(tidyverse)
library(bayesplot)
library(tidybayes)
library(rstanarm)
library(loo)
options(mc.cores = parallel::detectCores())

Alex2b2 <- read.csv("tempAlex.csv", row.names = 1) %>% 
  separate(sample, into = c("Site","Month","Year", "u"), remove = F) %>% 
  select(-u) %>% 
  mutate(Month = as.numeric(Month),
         SeasonBinary = ifelse(Season == "Summer", 1 , 2)) %>% 
  filter(Taxon == "Alexandrium_2b2",
         pH > 7) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity))


Alex3fc <- read.csv("tempAlex.csv", row.names = 1) %>% 
  separate(sample, into = c("Site","Month","Year", "u"), remove = F) %>% 
  select(-u) %>% 
  mutate(Month = as.numeric(Month),
         SeasonBinary = ifelse(Season == "Summer", 1 , 2)) %>% 
  filter(Taxon == "Alexandrium_3fc",
         pH > 7) %>% 
  mutate(TempStd = (Temperature - mean(Temperature))/sd(Temperature),
         pHStd = (pH - mean(pH))/sd(pH),
         SalinityStd = (Salinity - mean(Salinity))/sd(Salinity))

##############################################
#Both Alex species are seasonal. 
#Fit models in which pH and temperature and salinity might influence their occurrence, w 
#parameters varying by season. Then test model fits. 
#Standardizing the variables helps with fitting and w parameter interpretation, since all params will be on the same scale. 

#Alex2b2

arm.fit1 <- stan_glmer(presence ~ (1 + pHStd | Season),
                       data = Alex2b2,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit2 <- stan_glmer(presence ~ (1 + TempStd | Season),
                       data = Alex2b2,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit3 <- stan_glmer(presence ~ (1 + SalinityStd | Season),
                       data = Alex2b2,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit4 <- stan_glmer(presence ~ TempStd + (1 + pHStd | Season),
                       data = Alex2b2,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit5 <- stan_glmer(presence ~ pHStd + (1 + TempStd | Season),
                       data = Alex2b2,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit6 <- stan_glmer(presence ~ pHStd + SalinityStd + (1 + TempStd | Season),
                       data = Alex2b2,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)

#without season; just a glm, not a hierarchical model
arm.fit7 <- stan_glm(presence ~ pHStd + TempStd,
                       data = Alex2b2,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit8 <- stan_glm(presence ~ pHStd * TempStd,
                     data = Alex2b2,
                     family = "binomial",
                     prior_intercept = normal(0, 1),
                     prior = normal(0,1),
                     iter = 1000,
                     chains = 4)

#threshold question: 
#does salinity matter at all?
plot(arm.fit3)  #maybe a bit, so keep it for now.

#which model fits are best?
loo_compare(waic(arm.fit1),
            waic(arm.fit2),
            waic(arm.fit3),
            waic(arm.fit4),
            waic(arm.fit5),
            waic(arm.fit6),
            waic(arm.fit7), #model likes the hierarchical structure w season
            waic(arm.fit8)  #model likes the hierarchical structure w season
            )
#model 5 wins overall: an overall effect of pH, with a season-specific baseline probability of occurrence (intercept) and effect of temperature
#interestingly, model 1 is close behind, considering only pH.

#moving forward w model5:

#see https://mjskay.github.io/tidybayes/index.html
Alex2b2 %>% 
  add_fitted_draws(arm.fit5, n = 1000) %>% 
  ggplot(aes(x = Temperature, y = presence, color = Season)) +
  geom_point() +
  facet_grid(~ Season, scales = "free_x") +
  stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
  scale_fill_brewer()


Alex2b2 %>% 
  add_fitted_draws(arm.fit5, n = 1000) %>% 
  ggplot(aes(x = pH, y = presence, color = Season)) +
  geom_point() +
  facet_grid(~ Season, scales = "free_x") +
  stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
  scale_fill_brewer()


##########################################################################################
# same analyses w Alex3fc:

arm.fit11 <- stan_glmer(presence ~ (1 + pHStd | Season),
                       data = Alex3fc,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit12 <- stan_glmer(presence ~ (1 + TempStd | Season),
                       data = Alex3fc,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit13 <- stan_glmer(presence ~ (1 + SalinityStd | Season),
                       data = Alex3fc,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit14 <- stan_glmer(presence ~ TempStd + (1 + pHStd | Season),
                       data = Alex3fc,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit15 <- stan_glmer(presence ~ pHStd + (1 + TempStd | Season),
                       data = Alex3fc,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
arm.fit16 <- stan_glmer(presence ~ pHStd + SalinityStd + (1 + TempStd | Season),
                       data = Alex3fc,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1),
                       iter = 1000,
                       chains = 4)
#non-hier
arm.fit17 <- stan_glm(presence ~ TempStd,
                      data = Alex3fc,
                      family = "binomial",
                      prior_intercept = normal(0, 1),
                      prior = normal(0,1),
                      iter = 1000,
                      chains = 4)
arm.fit18 <- stan_glm(presence ~ pHStd + TempStd,
                     data = Alex3fc,
                     family = "binomial",
                     prior_intercept = normal(0, 1),
                     prior = normal(0,1),
                     iter = 1000,
                     chains = 4)
arm.fit19 <- stan_glm(presence ~ pHStd * TempStd,
                     data = Alex3fc,
                     family = "binomial",
                     prior_intercept = normal(0, 1),
                     prior = normal(0,1),
                     iter = 1000,
                     chains = 4)



#threshold questions: 
#does salinity matter at all?
plot(arm.fit13)  #no

#which model fits are best?
loo_compare(waic(arm.fit11),
            waic(arm.fit12),
            waic(arm.fit13),
            waic(arm.fit14),
            waic(arm.fit15),
            waic(arm.fit16),
            waic(arm.fit17),
            waic(arm.fit18),
            waic(arm.fit19)
)
#model 12 wins, w a seasonal difference in intercept (baseline probability of presence) and effect of temperature
#model 15 is next best, which is the same model as for Alex2b2, and includes pH

#since temperature is the only relevant variable here, we can plot it in one plot
Alex3fc %>% 
  add_fitted_draws(arm.fit12, n = 1000) %>% 
  ggplot(aes(x = Temperature, y = presence, color = Season)) +
  geom_point() +
  facet_grid(~ Season, scales = "free_x") +
  stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
  scale_fill_brewer()





