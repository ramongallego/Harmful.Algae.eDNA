library(tidyverse)
library(bayesplot)
library(tidybayes)
library(rethinking)
library(rstanarm)
library(loo)
library(lme4)
options(mc.cores = parallel::detectCores())

tempAlex <- read.csv("tempAlex.csv", row.names = 1)

tempAlex <- tempAlex %>% 
  separate(sample, into = c("Site","Month","Year", "u"), remove = F) %>% 
  select(-u) %>% 
  mutate(Month = as.numeric(Month),
         MonthsFromJune = abs(Month -6),
         SeasonBinary = ifelse(Season == "Summer", 1 , 2)) %>% 
  filter(Taxon == "Alexandrium_2b2",
         pH > 7) %>% 
  mutate(TempCentered = Temperature - mean(Temperature),
         pHCentered = pH - mean(pH),
         SalinityCentered = Salinity - mean(Salinity))

mod1 <- glm(presence ~ TempCentered + pH + Season, 
            data = tempAlex[tempAlex$Taxon == "Alexandrium_2b2",],
            family = "binomial")

mod2 <- glm(presence ~ TempCentered + MonthsFromJune, 
            data = tempAlex[tempAlex$Taxon == "Alexandrium_2b2",],
            family = "binomial")

mod3 <- glm(presence ~ TempCentered + Season, 
            data = tempAlex[tempAlex$Taxon == "Alexandrium_2b2",],
            family = "binomial")

mod4 <- glm(presence ~ pH + MonthsFromJune, 
            data = tempAlex[tempAlex$Taxon == "Alexandrium_2b2",],
            family = "binomial")

mod5 <- glm(presence ~ pH + TempCentered, 
            data = tempAlex[tempAlex$Taxon == "Alexandrium_2b2",],
            family = "binomial")


AIC(mod1, mod2, mod3, mod4, mod5)


pairs(tempAlex[,c("MonthsFromJune", "pHCentered", "TempCentered")], pch = 19)


#lme4 can't fit this model; needs prior
# mod4 <- glmer(presence ~ (1 + TempCentered | Season), 
#             data = tempAlex[tempAlex$Taxon == "Alexandrium_2b2",],
#             family = "binomial") 

arm.fit1 <- stan_glmer(presence ~ pHCentered + (1 + TempCentered | Season),
                      data = tempAlex,
                      family = "binomial",
                      prior_intercept = normal(0, 1),
                      prior = normal(0,1))

arm.fit2 <- stan_glmer(presence ~ pHCentered + (0 + TempCentered | Season),
                       data = tempAlex,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1))

arm.fit3 <- stan_glmer(presence ~ (1 + TempCentered | Season),
                       data = tempAlex,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1))

arm.fit4 <- stan_glmer(presence ~ (1 + TempCentered | Season) + (0 + pHCentered | Season),
                       data = tempAlex,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1))

arm.fit5 <- stan_glmer(presence ~ TempCentered + (1 + pHCentered | Season),
                       data = tempAlex,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1))

arm.fit6 <- stan_glmer(presence ~ pHCentered + (1 + TempCentered | Season),
                       data = tempAlex,
                       family = "binomial",
                       prior_intercept = normal(0, 1),
                       prior = normal(0,1))



#pH is really giving us some information here.
#and that one positive detection in September (Summer) is pretty influential
loo_compare(waic(arm.fit5),
            waic(arm.fit6)
)

#it's clear that one positive detection in September is odd, and the model has a hard time w it
arm.loo <-
  loo(arm.fit1, cores = 12)
plot(arm.loo)  
tempAlex[39,]

arm.loo2 <-
  loo(arm.fit6, cores = 12)


#see https://mjskay.github.io/tidybayes/index.html
tempAlex %>% 
  add_fitted_draws(arm.fit6, n = 1000) %>% 
  ggplot(aes(x = Temperature, y = presence, color = Season)) +
    geom_point() +
    facet_grid(~ Season, scales = "free_x") +
    stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
  scale_fill_brewer()
  

tempAlex %>% 
  add_fitted_draws(arm.fit6, n = 1000) %>% 
  ggplot(aes(x = pH, y = presence, color = Season)) +
  geom_point() +
  facet_grid(~ Season, scales = "free_x") +
  stat_lineribbon(aes(y = .value), .width = c(.95, .5)) +
  scale_fill_brewer()
  
#Conclusion:
#We have only observed this dino between September and March. This is the colder part of the year, and so it appears to occur in colder water.
#But within that time period, Alex2b2 occurs at seasonal-average temperatures (+/- about 2 degrees C of mean), and tolerates a range of salinities. However, we have only observed it when pH was at or below the seasonal mean. 
#Thus, we are most likely to see Alex2b2 between September and March when pH is below about 8.18

a <- tempAlex %>% 
  filter(Month %in% c(9:12, 1:3)) %>% 
  mutate(TempCentered = Temperature - mean(Temperature),
         pHCentered = pH - mean(pH),
         SalinityCentered = Salinity - mean(Salinity))

b <- data.frame(
  presence = a$presence,
  pH_below_seasonal_mean = as.numeric(a$pHCentered < 0),
  T_above_seasonal_mean = as.numeric(a$TempCentered > 0)
)

xtabs(presence ~ pH_below_seasonal_mean + T_above_seasonal_mean, data =b)




ggplot(a, aes(x = TempCentered, y = pHCentered, color = as.factor(presence))) +
  geom_point()



tempAlex$armfit <- posterior_predict(arm.fit) %>% 
  colMeans()



#creating multilevel model here for rethinking/stan
f <- alist(
  presence ~ dbinom( 1 , p ),
  logit(p) <- b1*pHCentered + b2*MonthsFromJune,
  #a ~ dnorm(0, 1),
  b1 ~ normal(0, 1),
  b2 ~ normal(0, 1)
)

#fitting model w mcmc
ulam.fit <- ulam(f, 
                 data = tempAlex,
                 chains = 4)
precis(ulam.fit)

#eval priors
#(there's got to be a better way of doing this)
prior.ulam <- extract.prior(ulam.fit, n = 1000)

priorDraws <- link(ulam.fit , 
     post=prior.ulam , 
     data=tempAlex) %>% 
  as.data.frame()
colnames(priorDraws) <- paste0("datum", 1:ncol(priorDraws))
row.names(priorDraws) <- paste0("drawNo", 1:nrow(priorDraws))

priorDraws %>% 
  pivot_longer()



# cv_quap(fit, 
#         lno = 3,
#         cores = 4) #cross-validation


tempAlex$ulam_predicted <- link(ulam.fit) %>% 
    colMeans()

tempAlex$glm_predicted <- mod4$fitted.values


ggplot(tempAlex, aes(x = TempCentered, y = presence, color = Season)) +
  geom_point() +
  facet_grid(~Season, scales = "free_x") +
  geom_point(aes(x = TempCentered, y = glm_predicted), color = "grey") +
  geom_point(aes(x = TempCentered, y = armfit), color = "black")

