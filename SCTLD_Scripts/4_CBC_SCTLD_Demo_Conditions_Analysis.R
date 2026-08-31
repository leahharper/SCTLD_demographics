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
  TimeLevels <- c("Oct2019","Nov2019","Dec2019",
                  "Jan2020", "Feb2020", "Mar2020" ,"Apr2020","May2020", "Jun2020","Jul2020",
                  "Aug2020" ,"Sep2020", "Oct2020","Nov2020","Dec2020",
                  "Jan2021","Feb2021","Mar2021" ,"Apr2021","May2021", "Jun2021",
                  "Jul2021","Aug2021" ,"Sep2021", "Oct2021","Nov2021","Dec2021",
                  "Jan2022","Feb2022","Mar2022" ,"Apr2022","May2022", "Jun2022","Jul2022",
                  "Aug2022" ,"Sep2022", "Oct2022","Nov2022","Dec2022")
  
  
  sitecolors = c('CBC Central'='#882255','CBC Lagoon'='#117733','CBC30N'='#DDCC77',
                 'Curlew Patch' = '#332288', 'House Reef' ='#44AA99', 
                 'South Reef Central' = '#999933', 'SR30N' = '#CC6677')         # yellow
  
  
  isSig <- function(p) {
    ifelse(p > 0.01 & p < 0.05, "*",
           ifelse(p > 0.001 & p <= 0.01, "**",
                  ifelse(p <= 0.001, "***", "")))
    
  }
  ############################
  #DEMO
  ############################
  demo <- read.csv("CBC_demo_curated.csv") %>% select(-c(X.1, ID))
  cond <- read.csv("CBC_conditions_curated.csv") %>% select(-X)
  
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
  
  srad_check <- demolong %>% subset(time_point == "Oct2019" & species == "Siderastrea radians")
  
  demolong <- demolong %>% 
    mutate(survey = case_when(
      survey == "Oct2019" | survey == "Jan2020" ~ "Nov2019",
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
    subset(survey == "Nov2019") %>%
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
    select(species, Nov2019, May2022, Dec2022) %>%
    rename("Oct 2019/Jan 2020" = Nov2019, "May 2022" = May2022, "Dec 2022" = Dec2022)
    
  # exclude_table <- exclude %>%
  #   kbl(caption = "<span style='color: black;'> <b>Table S2.</b> Total counts by survey timepoint of 
  #       species dropped from models (fewer than 8 total observations for all
  #       except H. cuculatta and M. decactis)<span>",
  #       digits = 3,
  #       format = "html", booktabs = TRUE, longtable = TRUE) %>%
  #   kable_styling(latex_options = c("repeat_header")) %>%
  #   column_spec(column = 1, italic = TRUE)
  # save_kable(exclude_table, file = "Tables/Table S2.pdf")
  
  write.csv(exclude, "Tables/exclude_table.csv")
  
  demo_df <- demo_tl %>% subset(total_all > 7) %>% mutate_at(.vars = vars("time_point"),
                                                   .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) 
  
  
  levels(as.factor(demo_df$species))
  
  #determine thresholds for inclusion
  test <- demo_df %>% 
    group_by(species, survey) %>% summarize(n0 = sum(total == 0)) 
   
  totals <- demo_df %>% group_by(species, total_all) %>% summarize(n = n())

  demo_df <- demo_df %>% subset(species != "Helioseris cucullata" &
                                  species != "Madracis decactis") %>% droplevels()
  # 
  # 
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
  
  # densmod <- glmer(total ~ survey*species + (1|location_name), family = "poisson", 
  #                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
  #                  data = demo_df)
  # 
  # hist(resid(densmod)) 
  # Anova(densmod)
  # summary(densmod)
  # fixed_effects_coefficients <- summary(densmod)$coefficients
  # fixed_effects_coefficients
  # 
  # demo_df <- demo_df %>%
  #   mutate(survey_numeric = as.numeric(factor(survey, levels = c("November19", "May22", "December22")))) %>%
  #   mutate(survey_scaled = as.numeric(scale(as.numeric(survey_numeric))))
  # 
  # nbinommod <- mixed_model(fixed = total ~ survey + species, random = ~1 | location_name, 
  #                          data = demo_df, family = negative.binomial(), 
  #                          max_coef_value = 30,
  #                          control = list(iter_EM = 0))
  #hist(resid(nbinommod))
  
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
  
  
  # demo_table <- demowald %>%
  #   kbl(caption = "<span style='color: black;'> <b>Table S3.</b> Results of
  #       zero inflated negative binomial model testing effects of survey timepoint,
  #       species, and their interaction on total coral count in transect surveys, with
  #       <i>Madracis decactis</i> and <i>Helioseris cucullata</i> included.<span>",
  #       digits = 3,
  #       format = "html", booktabs = TRUE, longtable = TRUE) %>%
  #   kable_styling(latex_options = c("repeat_header"))
  # save_kable(demo_table, file = "Tables/Table S3.pdf")
  # 
  # write.csv(demowald, "Tables/TableS3.csv")
  
  demo_table <- demowald %>%
    kbl(caption = "<span style='color: black;'> <b>Table 5.</b> Results of
        zero inflated negative binomial model testing effects of survey timepoint,
        species, and their interaction on total coral count in transect surveys.<span>",
        digits = 3,
        format = "html", booktabs = TRUE, longtable = TRUE) %>%
    kable_styling(latex_options = c("repeat_header"))
  #save_kable(demo_table, file = "Tables/Table 5.pdf")
  #write.csv(demowald, "Tables/Table 5.csv")
  
  
  # Fit a reduced model for comparison
  # zibinom_reduced <- mixed_model(
  #   fixed = total ~ survey + species,  # No interaction
  #   random = ~1 | location_name, 
  #   data = demo_df, 
  #   family = zi.negative.binomial(), 
  #   zi_fixed = ~ species, 
  #   zi_random = NULL,
  #   max_coef_value = 30,
  #   control = list(iter_EM = 0)
  # )
  
  
  # Compare models
  
  AIC(zibinom)
  AIC(nbinommod)
  #AIC(zibinom_reduced)
  
  #anova(zibinom_reduced, zibinom)
  
  
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
    
    letters$Group <- recode(letters$Group, "Nov219" = "Nov2019",
                            "May222" = "May2022",
                            "Dec222" = "Dec2022")
    
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
    kable_styling(latex_options = c("repeat_header")) %>%
    column_spec(column = 2, italic = TRUE)
  #save_kable(demo_emmefftable, file = "Tables/Table S7.pdf")
  #write.csv(est_eff, "Tables/Table S7.csv")
  
  #####################################
  #Species-specific condition models----
  #####################################
  levels(as.factor(dem3$species))
  
  zero_tl <- dem3 %>% group_by(species) %>% summarize(total_tl = sum(n_tl)) %>% 
     subset(total_tl < 1)
  
  tldf_sub <- dem3 %>% subset(!(species %in% zero_tl$species)) 
  

  
  ##create binomial dataframe where each colony is 0 (healthy) or 1 (diseased)----
  # sum_demo_long <- tldf_sub %>% dplyr::select(event, survey, location_name, species, 
  #                                             prev_tl, n_healthy) %>%
  #   uncount(n_healthy) %>% mutate("tl_val"  = 0)
  # 
  # sum_tl_long <- tldf_sub %>% dplyr::select(event, survey, location_name,
  #                                           species, prev_tl, n_tl) %>%
  #   uncount(n_tl) %>% mutate("tl_val" = 1)
  # 
  # tl_df_long <- rbind(sum_demo_long, sum_tl_long) %>% arrange(prev_tl) %>% 
  #   droplevels()
  # 
  # levels(as.factor(tl_df_long$species))
  
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
  
 check <-  tldf_sub %>%
     group_by(survey, species) %>%
    summarise(total_tl = sum(n_tl), total_healthy = sum(n_healthy), n = n())
  
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
  qqnorm(resid(betabinom)) 
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
    kbl(caption = "<span style='color: black;'> <b>Table 2.</b> Results of
        beta binomial model testing effects of survey timepoint and
        species on disease prevalence in transect surveys.<span>",
        digits = 3,
        format = "html", booktabs = TRUE, longtable = TRUE) %>%
    kable_styling(latex_options = c("repeat_header"))
  #save_kable(cond_table, file = "Tables/Table 2.pdf")
  #write.csv(condwald, "Tables/Table 2.csv")
  
  library(emmeans)
  library(rcompanion)
  
  
  
  #pairwise surveys
  emm_surv <- emmeans(betabinom, ~ survey)
  simple_surv <- pairs(emm_surv, simple = "survey")
  pairwise_surv <- as.data.frame(pairs(emm_surv, simple = "survey")) %>%
    mutate(`p.value` =round(`p.value`, digits = 3)) %>%
    mutate(`p.value` = ifelse(`p.value` < 0.001 ,
                              "<0.001", `p.value` ))
  
  pair_surv_merge <- pairwise_surv %>%
    select(-df) %>% rename("SE_est" = SE)
  
  eff_surv <- as.data.frame(eff_size(emm_surv, sigma = 26, edf = Inf)) 
  
  est_eff_surv <- pair_surv_merge %>% left_join(eff_surv, by = "contrast") %>%
    select(-SE)  %>%
    mutate(sig = isSig(`p.value`))
  
  
  
  
  surv_cond_emmtable <- est_eff_surv %>%
    kbl(caption = "<span style='color: black;'> <b>Table S4.</b> Pairwise contrasts - tissue loss prevalence over time.
        Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
        digits = 3,
    format = "html", booktabs = TRUE, longtable = TRUE) %>%
    kable_styling(latex_options = c("repeat_header"))
  #save_kable(surv_cond_emmtable, file = "Tables/Table S4.pdf")
  #write.csv(est_eff_surv, "Tables/Table S4.csv")
  
  
  #pairwise species
  emm_spp <- emmeans(betabinom, ~ species)
  simple_spp <- pairs(emm_spp, simple = "species")
  pairwise_spp <- as.data.frame(pairs(emm_spp, simple = "species")) %>%
    mutate(`p.value` =round(`p.value`, digits = 3)) %>%
    mutate(`p.value` = ifelse(`p.value` < 0.001 ,
                              "<0.001", `p.value` ))
  
  pair_spp_merge <- pairwise_spp %>%
    select(-df) %>% rename("SE_est" = SE)
  
  eff_spp <- as.data.frame(eff_size(emm_spp, sigma = 26, edf = Inf)) 
  
  est_eff_spp <- pair_spp_merge %>% left_join(eff_spp, by = "contrast") %>%
    select(-SE)  %>%
    mutate(sig = isSig(`p.value`))
  
  
  spp_cond_emmtable <- est_eff_spp %>%
    kbl(caption = "<span style='color: black;'> <b>Table S5.</b> Pairwise contrasts - difference in tissue loss 
    prevalence between pairs of species across all timepoints.",
        digits = 3,
        format = "html", booktabs = TRUE, longtable = TRUE) %>%
    kable_styling(latex_options = c("repeat_header")) %>%
    column_spec(column = 1, italic = TRUE)
  #save_kable(spp_cond_emmtable, file = "Tables/Table S5.pdf")
  #write.csv(est_eff_spp, "Tables/Table S5.csv")
  
  
  desired_order <- c("Nov2019", "May2022", "Dec2022")
  
  lesionsum <- demo_tl %>% group_by(survey, species) %>%
    summarize(n_tl = sum(n_tl), n_healthy = sum(n_healthy)) %>%
    mutate(survey = factor(survey, levels = desired_order)) %>%
    arrange(survey) %>%
    pivot_wider(names_from = survey,
                values_from = c(n_tl, n_healthy)) %>%
    select(species, n_tl_Nov2019, n_healthy_Nov2019,
           n_tl_May2022, n_healthy_May2022, n_tl_Dec2022, n_healthy_Dec2022) %>%
    rename("Oct 2019/Jan 2020 TL" = n_tl_Nov2019, "Oct 2019/Jan 2020 Healthy" = n_healthy_Nov2019,
           "May 2022 TL" = n_tl_May2022, "May 2022 Healthy" = n_healthy_May2022,
           "Dec 2022 TL" = n_tl_Dec2022, "Dec 2022 Healthy" = n_healthy_Dec2022)
  
  lesion_summary_table <- lesionsum %>%
    kbl(caption = "<span style='color: black;'> <b>Table 1.</b> Total counts of colonies with tissue loss and
    healthy colonies (without tissue loss) by survey timepoint and species.
        <span>",
        digits = 3,
        format = "html", booktabs = TRUE, longtable = TRUE) %>%
    kable_styling(latex_options = c("repeat_header")) %>%
    column_spec(column = 1, italic = TRUE)
  #save_kable(lesion_summary_table, file = "Tables/Table 1.pdf")
  #write.csv(lesionsum, "Tables/Table 1.csv")
  
  
  lesionsum2 <- demo_tl %>% group_by(survey, location_name, species) %>%
    summarize(n_tl = sum(n_tl), n_healthy = sum(n_healthy)) %>%
    mutate(survey = factor(survey, levels = desired_order)) %>%
    arrange(survey) 
  
  lesionsum2$survey <- recode(lesionsum2$survey,
                              "Nov2019" = "Oct 2019/Jan 2020",
                              "May2022" = "May 2022",
                              "Dec2022" = "Dec 2022")
  
  #write.csv(lesionsum2, "lesion_summary.csv")

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
  geom_jitter(data = high, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.7, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "Jul2021", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("# Colonies", limits = c(0,205), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('Oct2019','Jan2020', 'Jul2021',
                                                'May2022','Dec2022'),
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
        strip.text = element_text(face = "italic"),
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
  geom_jitter(data = med, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.7, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "Jul2021", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("# Colonies", limits = c(0,135), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('Oct2019','Jan2020', 'Jul2021',
                                                'May2022','Dec2022'),
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
        strip.text = element_text(face = "italic"),
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
  geom_jitter(data = low, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.7, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "Jul2021", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("# Colonies", limits = c(0,14), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('Oct2019','Jan2020', 'Jul2021',
                                                'May2022','Dec2022'),
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
        strip.text = element_text(face = "italic"),
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
  geom_jitter(data = vlow, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.7, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "Jul2021", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("# Colonies", limits = c(0,6), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('Oct2019','Jan2020', 'Jul2021',
                                                'May2022','Dec2022'),
                   expand = expansion(add = c(5, 5.5))) + 
  scale_size_continuous("TL Prevalence", range = c(2,7)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        strip.text = element_text(face = "italic"),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

vlowdensp

legendp <- ggplot() +
  geom_jitter(data = vlow, aes(x = time_point, y = total, fill = location_name, size = prev_tl), alpha = 0.7, pch = 21,  
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
  geom_jitter(data = high, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.7, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "Jul2021", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("TL Prevalence", limits = c(0,0.23), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('Oct2019','Jan2020', 'Jul2021',
                                                'May2022','Dec2022'),
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
        strip.text = element_text(face = "italic"),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

prev_highdensp

prev_meddensp <- ggplot() +
  geom_jitter(data = med, aes(x = time_point, y = prev_tl, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  #geom_text(data = med, aes(x = survey, y = max_tl + 0.04, label = prev_letter)) +
  geom_boxplot(data = med, aes(x = survey, y = prev_tl), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = med, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.7, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "Jul2021", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("TL Prevalence", limits = c(0,0.55), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('Oct2019','Jan2020', 'Jul2021',
                                                'May2022','Dec2022'),
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
        strip.text = element_text(face = "italic"),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

prev_meddensp

prev_lowdensp <- ggplot() +
  geom_jitter(data = low, aes(x = time_point, y = prev_tl, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  #geom_text(data = low, aes(x = survey, y = max_tl + 0.01, label = prev_letter)) +
  geom_boxplot(data = low, aes(x = survey, y = prev_tl), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = low, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.7, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "Jul2021", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("TL Prevalence", limits = c(0,1), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('Oct2019','Jan2020', 'Jul2021',
                                                'May2022','Dec2022'),
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
        strip.text = element_text(face = "italic"),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

prev_lowdensp

pstr_avg <- low %>% subset(time_point == "May22" & species == "Pseudodiploria strigosa") %>%
  summarize(sum = sum(adult), sumtl = sum(n_tl)) %>%
    mutate(prev = sumtl/sum)


prev_vlowdensp <- ggplot() +
  geom_jitter(data = vlow, aes(x = time_point, y = prev_tl, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  #geom_text(data = vlow, aes(x = survey, y = max_tl + 0.01, label = prev_letter)) +
  geom_boxplot(data = vlow, aes(x = survey, y = prev_tl), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = vlow, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.7, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "Jul2021", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~species, nrow = 1, scales = "free") +
  scale_y_continuous("TL Prevalence", limits = c(0,1), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('Oct2019','Jan2020', 'Jul2021',
                                                'May2022','Dec2022'),
                   expand = expansion(add = c(5, 5.5))) + 
  scale_size_continuous("Total Colonies", range = c(2,7), limits = c(0,185)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        strip.text = element_text(face = "italic"),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

prev_vlowdensp

prev_legendp <- ggplot() +
  geom_jitter(data = vlow, aes(x = time_point, y = prev_tl, fill = location_name, size = total), alpha = 0.7, pch = 21,  
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

png("Figures/current/Fig S4.png", width = 9, height = 8, units = "in", res = 400)
prevplot3
dev.off()


#check SSID isolates over time
iso <- demo %>% subset(scientific_name == "Siderastrea siderea" |
                         scientific_name == "Montastraea cavernosa"|
                         scientific_name == "Orbicella spp."|
                       scientific_name == "Pseudidploria strigosa"|
                       scientific_name == "Porites astreoides"|
                         scientific_name == "Agaricia agaricites"|
                         scientific_name == "Agaricia tenuifolia") %>%
  group_by(event, scientific_name, location_name, time_point, survey) %>%
  summarize(isos = sum(isolate_1_to_4_cm), total = sum(total)) %>%
  mutate(prop_iso = isos/total) %>%
  mutate_at(.vars = vars("time_point"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
  mutate_at(.vars = vars("survey"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

iso_p <- ggplot() +
  geom_jitter(data = iso, aes(x = time_point, y = prop_iso, fill = location_name), pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  #geom_text(data = ssid_iso, aes(x = survey, y = max+0.9, label = Letter)) +
  geom_boxplot(data = iso, aes(x = survey, y = prop_iso), fill = "gray80", width = 7, outlier.shape = NA) +
  geom_jitter(data = iso, aes(x = time_point, y = prop_iso, fill = location_name), alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  scale_y_continuous("Proportion Isolates <4cm", limits = c(0,0.3), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5.5))) + 
  facet_wrap(~scientific_name, nrow = 2, scales = "free") +
  #scale_size_continuous("Total Adult Colonies", range = c(2,7)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

iso_p

ssid <- demo %>% subset(scientific_name == "Siderastrea siderea") %>%
  select(-c(X.1, day, month, year, total, juveniles_1_to_4_cm,
            total_lg, total_small, adult)) %>%
  group_by(event, scientific_name, location_name, time_point, survey) %>%
  summarise(across(where(is.numeric), sum, .names = "total_{.col}")) %>%
pivot_longer(
  cols = starts_with("total_"),
  names_to = "size_class",
  values_to = "total"
) %>%
  mutate_at(.vars = vars("survey"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))



ssid$size_class <- recode(ssid$size_class,
                          "total_isolate_1_to_4_cm" = "4",
                          "total_X5_to_10_cm" = "10",
                          "total_X11_to_20_cm" = "20",
                          "total_X21_to_40_cm" = "40",
                          "total_X41_to_80_cm" = "80",
                          "total_over_80_cm" = "120")

ssid$size_class <- as.numeric(ssid$size_class)

ssid_sum <- ssid %>% group_by(survey, size_class) %>%
  summarize(total = sum(total))

ssid_size <- ssid %>% 
  uncount(total) 

#size curve plots----
library(viridis)

sizespecs <- ggplot() +
  geom_bar(data = ssid_size, aes(x= size_class, fill=survey),
               alpha=0.4, adjust = 2) +
  facet_grid(~survey) +
  #scale_color_viridis("", discrete=TRUE, end = 0.92) +  # end = 0.9 makes yellow darker
  #scale_fill_viridis("Survey", discrete=TRUE, end = 0.92) +
  guides(color = "none") + 
  ylab("Number in Class") +
  scale_x_continuous("Size Class (cm)") +
  #scale_x_continuous("Size Class (cm)",breaks = c(4, 10, 20, 40, 80, 120)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

#png("report_Belize/SizeSpec_Belize.png",width = 8, height = 5, units = "in", res = 300)
sizespecs
#()

check <- demo %>% subset(scientific_name == "Siderastrea siderea") %>%
  select(diver, location_name, survey, X21_to_40_cm, X41_to_80_cm, over_80_cm) %>%
  group_by(location_name, survey) %>%
  summarize(total_under = sum(X21_to_40_cm), total_q = sum(X41_to_80_cm), total_over = sum(over_80_cm)) 


