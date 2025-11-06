#########
#model general cover categories
#########
mod1 <- glmer(cov_prop ~ Survey*Label_General +
                (1|SiteName),
              weights=npoints,
              data=covgen,family="binomial")
#########
#models species specific cover
#########
mod3 <- glmer(cov_prop ~ Survey*Species +
                (1|SiteName),
              weights=npoints,
              data=sto2,family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

#########
#model change in coral count
#zero inflated negative binomial using GLMMadaptive zero inflated "hurdle"
#########
zibinom <- mixed_model(fixed = total ~ survey * species, random = ~1 | location_name, 
                       data = demo_df, family = zi.negative.binomial(), 
                       zi_fixed = ~ species, zi_random = NULL,
                       max_coef_value = 30,
                       control = list(iter_EM = 0))

##########
#model disease prevalence
#beta binomial regression in GLMMadaptative that incorporates counts of healthy and tl
##########
betabinom <- mixed_model(
  fixed = cbind(n_tl, n_healthy) ~ survey * species,
  random = ~ 1 | location_name,
  data = tldf_sub,
  family = beta.binomial(),
  max_coef_value = 30,
  control = list(iter_EM = 0))
