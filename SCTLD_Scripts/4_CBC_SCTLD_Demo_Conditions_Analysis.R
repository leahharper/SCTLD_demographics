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
exclude <- demo_tl %>% subset(total_all < 8) %>%
  group_by(species, survey) %>%
  summarize("total count" = sum(total)) %>%
  pivot_wider(id_cols = c(species), names_from = survey, 
              values_from = "total count") %>%
  select(species, November19, May22, December22) %>%
  rename("Oct 2019/Jan 2020" = November19, "May 2022" = May22, "Dec 2022" = December22)
  
exclude_table <- exclude %>%
  kbl(caption = "<span style='color: black;'> <b>Table S2.</b> Total counts by survey timepoint of 
      species dropped from models (fewer than 8 total observations <span>",
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
  kbl(caption = "<span style='color: black;'> <b>Table 1.</b> Results of
      zero inflated negative binomial model testing effects of survey timepoint, 
      species, and their interaction on total coral count in transect surveys.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(demo_table, file = "Tables/Table 1.pdf")


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


isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}

emm <- emmeans(zibinom, ~ survey*species)
simple <- pairs(emm, simple = "survey")
pairwise <- as.data.frame(pairs(emm, simple = "survey")) 
eff <- as.data.frame(eff_size(emm, sigma = 26, edf = Inf)) 


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

pairwise <- pairwise %>%
  mutate(`p.value` =round(`p.value`, digits = 3)) %>%
  mutate(`p.value` = ifelse(`p.value` < 0.001 ,
                            "<0.001", `p.value` ))


demo_emmtable <- pairwise %>%
  kbl(caption = "<span style='color: black;'> <b>Table S3.</b> Pairwise contrasts - colony counts over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(demo_emmtable, file = "Tables/Table S3.pdf")

demo_efftable <- eff %>%
  kbl(caption = "<span style='color: black;'> <b>Table S4.</b> Effect sizes - colony counts over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(demo_efftable, file = "Tables/Table S4.pdf")

#####################################
#Species-specific condition models----
#####################################
levels(as.factor(dem3$species))

zero_tl <- dem3 %>% group_by(species) %>% summarize(total_tl = sum(n_tl)) %>% 
   subset(total_tl < 1)

tldf_sub <- dem3 %>% subset(!(species %in% zero_tl$species))


##create binomial dataframe where each colony is 0 (healthy) or 1 (diseased)----
sum_demo_long <- tldf_sub %>% dplyr::select(event, survey, location_name, species, 
                                            prev_tl, n_healthy) %>%
  uncount(n_healthy) %>% mutate("tl_val"  = 0)

sum_tl_long <- tldf_sub %>% dplyr::select(event, survey, location_name,
                                          species, prev_tl, n_tl) %>%
  uncount(n_tl) %>% mutate("tl_val" = 1)

tl_df_long <- rbind(sum_demo_long, sum_tl_long) %>% arrange(prev_tl) %>% 
  droplevels()

levels(as.factor(tl_df_long$species))

library(GLMMadaptive)


# zibinom <- mixed_model(fixed = tl_val ~ survey + species, random = ~1 | location_name, 
#                        data = tl_df_long, family = zi.binomial(), 
#                        zi_fixed = ~species, zi_random = ~1|location_name,
#                        max_coef_value = 30,
#                        control = list(iter_EM = 0))
# 
# AIC(zibinom)
# 
# hist(resid(zibinom))
# summary(zibinom)

summary(tldf_sub$prev_tl)
table(tldf_sub$prev_tl == 0)  # All healthy
table(tldf_sub$prev_tl == 1)  #All disease

#tldf_sub <- tldf_sub %>%
  #group_by(survey, species) %>%
  #subset(prev_tl < 1 & prev_tl > 0) %>% ungroup()

betabinom <- mixed_model(
  fixed = cbind(n_tl, n_healthy) ~ survey + species,
  random = ~ 1 | location_name,
  data = tldf_sub,
  family = beta.binomial(),
  max_coef_value = 30,
  control = list(
    iter_EM = 20)
)

hist(resid(betabinom))
AIC(betabinom)
summary(betabinom)

#attempt specifying starting vals
# betabinom_simple <- mixed_model(
#   fixed = cbind(n_tl, n_healthy) ~ survey + species,
#   random = ~ 1 | location_name,
#   data = tldf_sub,
#   family = beta.binomial()
# )
# 
# betabinom <- mixed_model(
#   fixed = cbind(n_tl, n_healthy) ~ survey * species,
#   random = ~ 1 | location_name,
#   data = tldf_sub,
#   family = beta.binomial(),
#   max_coef_value = 50,
#   initial_values = list(
#     betas = coef(simple_glm),
#     D = matrix(0.1)  # Small positive value for random effect variance
#   )
# )
# 
# AIC(betabinom)
# hist(resid(betabinom))
# summary(betabinom)

#Wald Tests, Condition ----

# Get coefficient names first
coef_names <- names(fixef(betabinom))
print(coef_names)  # Check what you actually have
n_coef <- length(coef_names)

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

# Run tests
wald_survey <- anova(betabinom, L = L_survey)
wald_species <- anova(betabinom, L = L_species)

#extract results
condwald1 <- as.data.frame(wald_survey$aovTab.L) %>%
  mutate(Effect = "Survey")
condwald2 <- as.data.frame(wald_species$aovTab.L) %>%
  mutate(Effect = "Species")
condwald <- rbind(condwald1,condwald2) %>%
  select(Effect, Chisq, df, "Pr(>|Chi|)") %>%  
  mutate(`Pr(>|Chi|)` = ifelse(`Pr(>|Chi|)` < 0.001 ,
              "<0.001", `Pr(>|Chi|)` )) %>%
  mutate(sig = isSig(`Pr(>|Chi|)`)) %>%
  remove_rownames()


cond_table <- condwald %>%
  kbl(caption = "<span style='color: black;'> <b>Table 3.</b> Results of
      beta binomial model testing effects of survey timepoint and
      species on disease prevalence in transect surveys.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(cond_table, file = "Tables/Table 3.pdf")


library(emmeans)
library(rcompanion)


isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}

#pairwise surveys
emm_surv <- emmeans(betabinom, ~ survey)
simple_surv <- pairs(emm_surv, simple = "survey")
pairwise_surv <- as.data.frame(pairs(emm_surv, simple = "survey")) %>%
  mutate(`p.value` =round(`p.value`, digits = 3)) %>%
  mutate(`p.value` = ifelse(`p.value` < 0.001 ,
                            "<0.001", `p.value` )) %>%
  mutate(sig = isSig(`p.value`))

eff_surv <- as.data.frame(eff_size(simple_surv, sigma = 26, edf = Inf)) 


surv_cond_emmtable <- pairwise_surv %>%
  kbl(caption = "<span style='color: black;'> <b>Table S5.</b> Pairwise contrasts - tissue loss prevalence over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
      digits = 3,
  format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(surv_cond_emmtable, file = "Tables/Table S5.pdf")


surv_cond_efftable <- eff_surv %>%
  kbl(caption = "<span style='color: black;'> <b>Table S6.</b> Effect sizes - tissue loss prevalence over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(surv_cond_efftable, file = "Tables/Table S6.pdf")

#pairwise species
emm_spp <- emmeans(betabinom, ~ species)
simple_spp <- pairs(emm_spp, simple = "species")
pairwise_spp <- as.data.frame(pairs(emm_spp, simple = "species")) %>%
  mutate(`p.value` =round(`p.value`, digits = 3)) %>%
  mutate(`p.value` = ifelse(`p.value` < 0.001 ,
                            "<0.001", `p.value` )) %>%
  mutate(sig = isSig(`p.value`))


eff_spp <- as.data.frame(eff_size(simple_spp, sigma = 26, edf = Inf))


spp_cond_emmtable <- pairwise_spp %>%
  kbl(caption = "<span style='color: black;'> <b>Table S7.</b> Pairwise contrasts - difference in tissue loss 
  prevalence between pairs of species across all timepoints.",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(spp_cond_emmtable, file = "Tables/Table S7.pdf")


spp_cond_efftable <- eff_spp %>%
  kbl(caption = "<span style='color: black;'> <b>Table S8.</b> Effect sizes - tissue loss prevalence among species.
      <span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(spp_cond_efftable, file = "Tables/Table S8.pdf")

desired_order <- c("November19", "May22", "December22")

lesionsum <- demo_tl %>% group_by(survey, species) %>%
  summarize(n_tl = sum(n_tl), n_healthy = sum(n_healthy)) %>%
  mutate(survey = factor(survey, levels = desired_order)) %>%
  arrange(survey) %>%
  pivot_wider(names_from = survey,
              values_from = c(n_tl, n_healthy)) %>%
  select(species, n_tl_November19, n_healthy_November19,
         n_tl_May22, n_healthy_May22, n_tl_December22, n_healthy_December22) %>%
  rename("Oct 2019/Jan 2020 TL" = n_tl_November19, "Oct 2019/Jan 2020 Healthy" = n_healthy_November19,
         "May 2022 TL" = n_tl_May22, "May 2022 Healthy" = n_healthy_May22,
         "Dec 2022 TL" = n_tl_December22, "Dec 2022 Healthy" = n_healthy_December22)

lesion_summary_table <- lesionsum %>%
  kbl(caption = "<span style='color: black;'> <b>Table 2.</b> Total counts of colonies with tissue loss and
  healthy colonies (without tissue loss) by survey timepoint and species.
      <span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(lesion_summary_table, file = "Tables/Table 2.pdf")

#Generates sig letters for disease, cant do without interaction
# 
# df <- data.frame()
# 
# Spec <- levels(as.factor(pairwise_spp$species))
# 
# for(current_Spec in Spec) {
#   
#   Spec_df <- pairwise_spp %>% subset(species == current_Spec) 
#   
#   pairwise_spp$sig <- isSig(pairwise_spp$p.value)
#   
#   letters <- as.data.frame(cldList(p.value ~ contrast,
#                                    data = Spec_df,
#                                    threshold = 0.05)) %>%
#     mutate("Species" = current_Spec)
#   
#   df <- df %>%
#     bind_rows(letters)
# }
# 
# prevlet <- df %>% unite("event", c("Species", "Group")) %>%
#   rename("prev_letter" = Letter)
# 
# dem3 <- dem3 %>% left_join(prevlet, by = "event") %>% mutate(prev_letter = ifelse(species != "Agaricia tenuifolia" & species != "Siderastrea siderea" ,
#                                                                                   "", prev_letter))
# 
# 
 dem3 <- dem3 %>% mutate(Letter = ifelse(species == "Agaricia tenuifolia" | species == "Eusmilia fastigiata" |
                                             species == "Porites astreoides" | species == "Meandrina meandrites" |
                                             species == "Diploria labyrinthiformis"| species == "Pseudodiploria strigosa"|
                                             species == "Stephanocoenia intersepta",
                                           "", Letter)) %>%
   mutate_at(.vars = vars("time_point"),
             .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
   mutate_at(.vars = vars("survey"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

dem3$species <- recode(dem3$species,
                       "Porites porities" = "Porites porites")

dem3$survey <- as.factor(dem3$survey)
############################################
#plot count with bubble scaled to prevalence----
############################################
high <- dem3 %>% subset(species == "Agaricia agaricites" | species == "Agaricia tenuifolia"|
                           species == "Porites astreoides")

highdensp <- ggplot() +
  geom_jitter(data = high, aes(x = time_point, y = total, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  geom_text(data = high, aes(x = survey, y = max+20, label = Letter)) +
  geom_boxplot(data = high, aes(x = survey, y = total), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = high, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("# Colonies", limits = c(0,205), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_size_continuous("TL Prevalence", range = c(2,7)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

highdensp

med <- dem3 %>% subset(species == "Orbicella spp." | species == "Porites porites"|
                          species == "Siderastrea siderea")

meddensp <- ggplot() +
  geom_jitter(data = med, aes(x = time_point, y = total, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  geom_text(data = med, aes(x = survey, y = max+15, label = Letter)) +
  geom_boxplot(data = med, aes(x = survey, y = total), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = med, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("# Colonies", limits = c(0,135), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_size_continuous("TL Prevalence", range = c(2,7)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

meddensp


low <- dem3 %>% subset(species == "Montastraea cavernosa"|
                          species == "Pseudodiploria strigosa" | species == "Stephanocoenia intersepta")


lowdensp <- ggplot() +
  geom_jitter(data = low, aes(x = time_point, y = total, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  geom_text(data = low, aes(x = survey, y = max+2, label = Letter)) +
  geom_boxplot(data = low, aes(x = survey, y = total), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = low, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("# Colonies", limits = c(0,14), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_size_continuous("TL Prevalence", range = c(2,7)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

lowdensp

vlow <- dem3 %>% subset(species == "Diploria labyrinthiformis" | 
                         species == "Meandrina meandrites" | species == "Eusmilia fastigiata" |
                         species == "Dichocoenia stokesii")

vlowdensp <- ggplot() +
  geom_jitter(data = vlow, aes(x = time_point, y = total, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  geom_text(data = vlow, aes(x = survey, y = max+0.9, label = Letter)) +
  geom_boxplot(data = vlow, aes(x = survey, y = total), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = vlow, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("# Colonies", limits = c(0,6), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5.5))) + 
  scale_size_continuous("TL Prevalence", range = c(2,7)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

vlowdensp

legendp <- ggplot() +
  geom_jitter(data = vlow, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  scale_size_continuous("TL Prevalence", range = c(2,7)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

library(gridExtra)
library(cowplot)

legend <- cowplot::get_legend(legendp)

demplot <- grid.arrange(highdensp, meddensp, lowdensp, ncol=1, nrow =3)

demplot2 <- cowplot::plot_grid(demplot, legend, rel_widths = c(3/4, 1/4), axis = 't', align = "v")

demplot2

demplot3 <- cowplot::plot_grid(demplot2, vlowdensp, rel_heights = c(7/10,3/10),
                               ncol = 1, axis = 't', align = "v")

demplot3

#tiff("Figures/demo_tl_plot.tif",width = 9, height = 8, units = "in", res = 400)
#demplot3
#dev.off()

png("Figures/current/Fig 4.png", width = 9, height = 8, units = "in", res = 400)
demplot3
dev.off()

############################################
#plot prevalence with bubble scaled to count----
############################################
prev_highdensp <- ggplot() +
  geom_jitter(data = high, aes(x = time_point, y = prev_tl, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  #geom_text(data = high, aes(x = survey, y = max_tl + 0.02, label = prev_letter)) +
  geom_boxplot(data = high, aes(x = survey, y = prev_tl), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = high, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("TL Prevalence", limits = c(0,0.23), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_size_continuous("Total Colonies", range = c(2,7), limits = c(0,185)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

prev_highdensp

prev_meddensp <- ggplot() +
  geom_jitter(data = med, aes(x = time_point, y = prev_tl, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  #geom_text(data = med, aes(x = survey, y = max_tl + 0.04, label = prev_letter)) +
  geom_boxplot(data = med, aes(x = survey, y = prev_tl), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = med, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("TL Prevalence", limits = c(0,0.55), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_size_continuous("Total Colonies", range = c(2,7), limits = c(0,185)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

prev_meddensp

prev_lowdensp <- ggplot() +
  geom_jitter(data = low, aes(x = time_point, y = prev_tl, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  #geom_text(data = low, aes(x = survey, y = max_tl + 0.01, label = prev_letter)) +
  geom_boxplot(data = low, aes(x = survey, y = prev_tl), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = low, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("TL Prevalence", limits = c(0,1), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_size_continuous("Total Colonies", range = c(2,7), limits = c(0,185)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

prev_lowdensp


prev_vlowdensp <- ggplot() +
  geom_jitter(data = vlow, aes(x = time_point, y = prev_tl, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  #geom_text(data = vlow, aes(x = survey, y = max_tl + 0.01, label = prev_letter)) +
  geom_boxplot(data = vlow, aes(x = survey, y = prev_tl), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = vlow, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("TL Prevalence", limits = c(0,1), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5.5))) + 
  scale_size_continuous("Total Colonies", range = c(2,7), limits = c(0,185)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

prev_vlowdensp

prev_legendp <- ggplot() +
  geom_jitter(data = vlow, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  scale_size_continuous("Total Colonies", range = c(2,7), limits = c(0,185)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

library(gridExtra)
library(cowplot)

prev_legend <- cowplot::get_legend(prev_legendp)

prev_plot <- grid.arrange(prev_highdensp, prev_meddensp, prev_lowdensp, ncol=1, nrow =3)

prevplot2 <- cowplot::plot_grid(prev_plot, prev_legend, rel_widths = c(3/4, 1/4), axis = 't', align = "v")

prevplot2

prevplot3 <- cowplot::plot_grid(prevplot2, prev_vlowdensp, rel_heights = c(7/10,3/10),
                               ncol = 1, axis = 't', align = "v")

prevplot3

#tiff("Figures/prev_tl_plot.tif",width = 9, height = 8, units = "in", res = 400)
#prevplot3
#dev.off()

png("Figures/current/Fig S3.png", width = 9, height = 8, units = "in", res = 400)
prevplot3
dev.off()


