library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)
library(kableExtra)
library(webshot2)

# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/GitHub/SCTLD_demographics/")

pd <- position_dodge(width = 0.93)
se<-function(x)sqrt(var(x)/length(x))

# Establish time levels
TimeLevels <- c("October19","November19","December19",
                "January20", "February20", "March20" ,"April20","May20", "June20","July20",
                "August20" ,"September20", "October20","November20","December20",
                "January21","February21","March21" ,"April21","May21", "June21",
                "July21","August21" ,"September21", "October21","November21","December21",
                "January22","February22","March22" ,"April22","May22", "June22","July22",
                "August22" ,"September22", "October22","November22","December22")


sitecolors = c('CBC Central'='red3','CBC Lagoon'='darkorchid','CBC30N'='gold1',
               'Curlew Patch' = 'blue', 'House Reef' ='seagreen4', 
               'South Reef Central' = 'orange', 'SR30N' = 'pink')


isSig <- function(p) {
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}
############################
#DEMO
############################
demo <- read.csv("CBC_demo_curated.csv")
cond <- read.csv("CBC_conditions_curated.csv")

demolong <- demo %>% 
  group_by(location_name, time_point, survey, event, scientific_name) %>%
  summarize(total = sum(total), total_lg = sum(total_lg),
            total_small = sum(total_small), adult = sum(adult)) %>%
  ungroup() %>%
  complete(scientific_name, nesting(location_name, time_point, survey, event), 
           fill = list(total = 0, total_lg = 0,
                       total_small = 0, adult = 0)) %>%
  rename("species" = scientific_name) %>%
  mutate_at(.vars = vars("time_point"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


demolong <- demolong %>% 
  mutate(survey = case_when(
    survey == "October19" | survey == "January20" ~ "November19",
    TRUE ~ as.factor(survey))) %>%
  unite("merge_event", c(survey,location_name,species))

tldf_long <- cond %>%
  complete(species, nesting(location_name, time_point, survey), 
           fill = list(n_tl = 0, prev_tl = 0,
                       total = 0, n_healthy = 0)) %>%
  select(-(event)) %>%
  unite("merge_event", c(survey,location_name,species), remove = FALSE)



demo_tl <- demolong %>% left_join(tldf_long, by = "merge_event") %>%
  dplyr::select(-(total.y)) %>%
  dplyr::select(-(time_point.y))

colnames(demo_tl) <- gsub("\\.x","",colnames(demo_tl))

#check thresholds ----

thresh <- demo_tl %>% group_by(species, survey) %>%
  summarize(count = sum(adult > 0), sum = sum(adult)) %>%
  subset(survey == "November19") %>%
  arrange(count)

library(lme4)
library(detectseparation)

############################
#Species-specific count models----
############################

denssum <- demo_tl %>% group_by(species) %>% summarize(total_all = sum(total)) %>%
  arrange(total_all)

demo_tl <- demo_tl %>% left_join(denssum, by = "species") %>% subset(species != "Agaricia spp.") %>% droplevels()

#drop rarest species ----
exclude <- demo_tl %>% subset(total_all < 8 |
                                species == "Helioseris cucullata"|
                                species == "Madracis decactis") %>%
  group_by(species, survey) %>%
  summarize("total count" = sum(total)) %>%
  pivot_wider(id_cols = c(species), names_from = survey, 
              values_from = "total count") %>%
  select(species, November19, May22, December22) %>%
  rename("Oct 2019/Jan 2020" = November19, "May 2022" = May22, "Dec 2022" = December22)

exclude_table <- exclude %>%
  kbl(caption = "<span style='color: black;'> <b>Table S2.</b> Total counts by survey timepoint of 
      species dropped from models (fewer than 8 total observations for all
      except H. cuculatta and M. decactis)<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(exclude_table, file = "Tables/Table S2.pdf")

demo_df <- demo_tl %>% subset(total_all > 7) %>% mutate_at(.vars = vars("time_point"),
                                                           .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) 


levels(as.factor(demo_df$species))

#determine thresholds for inclusion
test <- demo_df %>% 
  group_by(species, survey) %>% summarize(n0 = sum(total == 0)) 

totals <- demo_df %>% group_by(species, total_all) %>% summarize(n = n())

demo_df <- demo_df %>% subset(species != "Helioseris cucullata" &
                                species != "Madracis decactis") %>% droplevels()


#P-A Model

# pa_df <- demo_df %>% mutate(total = ifelse(total > 0, 1, total))
# 
# binommod <- mixed_model(fixed = total ~ survey * species, random = ~1 | location_name, 
#                          data = pa_df, family = zi.binomial(), 
#                         zi_fixed = ~ species, zi_random = NULL,
#                         max_coef_value = 30,
#                         control = list(iter_EM = 0))
# hist(resid(nbinommod))


hist(demo_df$total)


library(lme4)
library(GLMMadaptive)
#remotes::install_github("glmmTMB/glmmTMB/glmmTMB")

densmod <- glmer(total ~ survey*species + (1|location_name), family = "poisson", 
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 data = demo_df)

hist(resid(densmod)) 
Anova(densmod)
summary(densmod)
fixed_effects_coefficients <- summary(densmod)$coefficients
fixed_effects_coefficients

demo_df <- demo_df %>%
  mutate(survey_numeric = as.numeric(factor(survey, levels = c("November19", "May22", "December22")))) %>%
  mutate(survey_scaled = as.numeric(scale(as.numeric(survey_numeric))))

nbinommod <- mixed_model(fixed = total ~ survey + species, random = ~1 | location_name, 
                         data = demo_df, family = negative.binomial(), 
                         max_coef_value = 30,
                         control = list(iter_EM = 0))
hist(resid(nbinommod))

demo_df <- as.data.frame(demo_df)

zibinom <- mixed_model(fixed = total ~ survey * species, random = ~1 | location_name, 
                       data = demo_df, family = zi.negative.binomial(), 
                       zi_fixed = ~ species, zi_random = ~1|location_name,
                       max_coef_value = 30,
                       control = list(iter_EM = 0, optimizer = "optim", optim_method = "BFGS"))
hist(resid(zibinom))
summary(zibinom)

zibinom2 <- mixed_model(fixed = total ~ survey * species, random = ~1 | location_name, 
                        data = demo_df, family = zi.negative.binomial(), 
                        zi_fixed = ~ species, zi_random = NULL,
                        max_coef_value = 30,
                        control = list(iter_EM = 0, optimizer = "optim", optim_method = "BFGS"))
hist(resid(zibinom2))
summary(zibinom2)

anova(zibinom, zibinom2)


#Wald Tests, Demo ----
#note we did Wald tests because they are default for GLMMadaptive 
#model structures, not the most conservative though
library(car)

# Get coefficient names first
coef_names <- names(fixef(zibinom))
print(coef_names)  # Check what you actually have
n_coef <- length(coef_names)

# Find interaction positions (contains ":")
interaction_positions <- grep(":", coef_names)
print(interaction_positions)  # Verify these exist

# Create L matrices
# Survey main effects
survey_idx <- grep("survey", coef_names)
survey_main_idx <- survey_idx[!grepl(":", coef_names[survey_idx])]

L_survey <- matrix(0, nrow = length(survey_main_idx), ncol = n_coef)
for(i in seq_along(survey_main_idx)) {
  L_survey[i, survey_main_idx[i]] <- 1
}

# Species main effects
species_idx <- grep("species", coef_names)
species_main_idx <- species_idx[!grepl(":", coef_names[species_idx])]

L_species <- matrix(0, nrow = length(species_main_idx), ncol = n_coef)
for(i in seq_along(species_main_idx)) {
  L_species[i, species_main_idx[i]] <- 1
}

# Interaction effects
L_interaction <- matrix(0, nrow = length(interaction_positions), ncol = n_coef)
for(i in seq_along(interaction_positions)) {
  L_interaction[i, interaction_positions[i]] <- 1
}

# Run tests
wald_survey <- anova(zibinom, L = L_survey)
wald_species <- anova(zibinom, L = L_species)
wald_interaction <- anova(zibinom, L = L_interaction)

#extract results
demowald1 <- as.data.frame(wald_survey$aovTab.L) %>%
  mutate(Effect = "Survey")
demowald2 <- as.data.frame(wald_species$aovTab.L) %>%
  mutate(Effect = "Species")
demowald3 <- as.data.frame(wald_interaction$aovTab.L) %>%
  mutate(Effect = "Survey x Species")
demowald <- rbind(demowald1,demowald2,demowald3) %>%
  select(Effect, Chisq, df, "Pr(>|Chi|)") %>%
  mutate(`Pr(>|Chi|)` =round(`Pr(>|Chi|)`, digits = 3)) %>%
  mutate(`Pr(>|Chi|)` = ifelse(`Pr(>|Chi|)` < 0.001 ,
                               "<0.001", `Pr(>|Chi|)` )) %>%
  remove_rownames() %>%
  mutate(sig = isSig(`Pr(>|Chi|)`))


demo_table <- demowald %>%
  kbl(caption = "<span style='color: black;'> <b>Table 3.</b> Results of
      zero inflated negative binomial model testing effects of survey timepoint, 
      species, and their interaction on total coral count in transect surveys.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(demo_table, file = "Tables/Table 3.pdf")


# Fit a reduced model for comparison
zibinom_reduced <- mixed_model(
  fixed = total ~ survey + species,  # No interaction
  random = ~1 | location_name, 
  data = demo_df, 
  family = zi.negative.binomial(), 
  zi_fixed = ~ species, 
  zi_random = NULL,
  max_coef_value = 30,
  control = list(iter_EM = 0)
)


# Compare models

AIC(zibinom)
AIC(nbinommod)
AIC(zibinom_reduced)

anova(zibinom_reduced, zibinom)


sd(demo_df$total)

library(emmeans)
library(rcompanion)



emm <- emmeans(zibinom, ~ survey*species)
simple <- pairs(emm, simple = "survey")
pairwise <- as.data.frame(pairs(emm, simple = "survey")) 
pairwisemerge <- pairwise %>%
  separate(contrast, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>%
  mutate(contrast_labeled = paste(time1, species, "-", time2, species)) %>%
  select(-time1, -time2, -df) %>% rename("SE_est" = SE)

eff <- as.data.frame(eff_size(emm, sigma = 26, edf = Inf)) 

est_eff <- pairwisemerge %>% left_join(eff, by = c("contrast_labeled" = "contrast")) %>%
  select(-SE, -contrast_labeled)
est_eff$sig <- isSig(est_eff$p.value)


df <- data.frame()

Spec <- levels(as.factor(pairwise$species))

for(current_Spec in Spec) {
  
  Spec_df <- pairwise %>% subset(species == current_Spec) 
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = Spec_df,
                                   threshold = 0.05)) %>%
    mutate("Species" = current_Spec)
  
  df <- df %>%
    bind_rows(letters)
}

denslet <- df %>% unite("event", c("Species", "Group"))

dem <- demo_df %>% unite("event", c("species", "survey"), remove = FALSE)

dem2 <- dem %>% left_join(denslet, by = "event")

demsum <- dem2 %>% group_by(event) %>%
  summarize(max = max(total), N = sum(total),
            max_tl = max(prev_tl))

dem3 <- dem2 %>% left_join(demsum, by = "event")

# pairwise <- pairwise %>%
#   mutate(`p.value` =round(`p.value`, digits = 3)) %>%
#   mutate(`p.value` = ifelse(`p.value` < 0.001 ,
#                             "<0.001", `p.value` ))

est_eff <- est_eff %>%
  mutate(`p.value` =round(`p.value`, digits = 3)) %>%
  mutate(`p.value` = ifelse(`p.value` < 0.001 ,
                            "<0.001", `p.value` ))


demo_emmefftable <- est_eff %>%
  kbl(caption = "<span style='color: black;'> <b>Table S7.</b> Pairwise contrasts - colony counts over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019. SE_est is standard error 
      of the estimate<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(demo_emmefftable, file = "Tables/Table S7&8.pdf")
