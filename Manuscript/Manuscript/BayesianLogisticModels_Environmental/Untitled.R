
mydataset$InRange <- ifelse(mydataset$TempStd < 0 & mydataset$TempStd > -1.5, 1, 0)
tempdataset$InRange <- ifelse(tempdataset$TempStd < 0 & tempdataset$TempStd > -1.5, 1, 0)


arm.fit <- stan_glmer(Alexandrium_2b2 ~ 0 + (1 | InRange) + (1 | Ditylum_ba4),
                        data = tempdataset,
                        family = "binomial",
                        prior_intercept = normal(0, 1),
                        prior = normal(0,1),
                        iter = 1000,
                        chains = 4)
mydataset$pred <- posterior_predict(arm.fit) %>% 
  colMeans()

plot(presence ~ TempStd, data = mydataset)
points(pred ~ TempStd, data = mydataset, col= "red")

plot(Alexandrium_2b2 ~ TempStd, data = newdataset)


stan_nlm(presence ~ exp(-((TempStd - b)^2) / (2*c^2)),
           data = mydataset,
           family = "binomial",
           prior_intercept = normal(0, 1),
           prior = normal(0,1),
           iter = 1000,
           chains = 4)



