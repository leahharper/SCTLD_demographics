library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)



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



demo <- read.csv("CBC_demo_curated.csv")
cond <- read.csv("CBC_Coral_Conditions.csv")

cond <- cond %>% subset(year < 2023)


#################
#Conditions
#################
levels(as.factor(cond$code))

aata <- cond %>% subset(code == "dlab")

cond <- cond %>% 
  mutate(scientific_name = case_when(
    code == "aaga" ~ "Agaricia agaricites",
    code == "aata" ~ "Agaricia agaricites",
    code == "acer" ~ "Acropora cervicornis",
    code == "aten" ~ "Agaricia tenuifolia",
    code == "cnat" ~ "Colpophyllia natans",
    code == "dlab" ~ "Diploria labyrinthiformis",
    code == "dsto" ~ "Dichocoenia stokesii",
    code == "dstr" ~ "Pseudodiploria strigosa",
    code == "efas" ~ "Eusmilia fastigiata",
    code == "ffra" ~ "Favia fragum",
    code == "hcuc" ~ "Helioseris cucullata",
    code == "mali" ~ "Mycetophyllia aliciae",
    code == "mcav" ~ "Montastraea cavernosa",
    code == "oann" ~ "Orbicella annularis",
    code == "ofav" ~ "Orbicella faveolata",
    code == "past" ~ "Porites astreoides",
    code == "ppor" ~ "Porites porites",
    code == "pstr" ~ "Pseudodiploria strigosa",
    code == "sdstr" ~ "Pseudodiploria strigosa",
    code == "sint" ~ "Stephanocoenia intersepta",
    code == "srad" ~ "Siderastrea radians",
    code == "ssid" ~ "Siderastrea siderea",
    code == "unks" ~ "Unknown coral",
    TRUE ~ as.factor(scientific_name)))


levels(as.factor(cond$year))

cond$year_a <- cond$year
cond$month_a <- cond$month
cond$day_a <- cond$day
cond$location_name_a <- cond$location_name


cond <- cond %>% unite("date", c("year_a","month_a","day_a"), sep = "-")

cond$date <- as.Date(cond$date)

cond$year_a <- cond$year
cond$month_a <- cond$month

cond <- cond %>% unite("time_point", c("year_a","month_a"), sep = "-")

levels(as.factor(cond$time_point))

cond$time_point <- recode(cond$time_point, "2019-10" = "October19",
                          "2022-12" = "December22", 
                          "2022-5" = "May22",
                          "2020-1" = "January20")

cond <- cond %>% mutate("survey" = time_point)
cond$survey <- recode(cond$survey, "October19" = "November19",
                      "January20" = "November19")

cond <- cond %>% subset(location_name != "Tobacco Reef")


cond$condition_code <- recode(cond$condition_code, "CLP;CLB" = "CLP/CLB")

cond <- cond %>% subset(condition_code != "")

cond$condition_code_a <- cond$condition_code

cond2 <- cond %>% separate_wider_delim(condition_code_a, "/", names = c("cond_1", "cond_2"),
                                       too_few = "align_start")

levels(as.factor(cond$percent_affected))

cond3 <- cond2 %>% mutate(percent_affected = ifelse(percent_affected == 1,
                                                    "A", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 10,
                                   "B", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 15,
                                   "B", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 20,
                                   "B", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 25,
                                   "B", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 30,
                                   "B", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 5,
                                   "A", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 50,
                                   "B", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 60,
                                   "C", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 70,
                                   "C", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 75,
                                   "C", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 80,
                                   "C", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 85,
                                   "C", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 90,
                                   "C", percent_affected)) %>%
  mutate(percent_affected = ifelse(percent_affected == 95,
                                   "C", percent_affected))

cond3$percent_affected <- gsub(";", "/", cond3$percent_affected)

cond3$percent_affected_a <- cond3$percent_affected

cond4 <- cond3 %>% separate_wider_delim(percent_affected_a, "/", names = c("perc_1", "perc_2"),
                                        too_few = "align_start")

levels(as.factor(cond4$distribution))

cond4$distribution <- gsub(";", "/", cond4$distribution)

pr <- cond4 %>% subset(distribution == "P"|distribution == "R"|distribution == "DA")

cond4$distribution <- recode(cond4$distribution,
                             "R" = "F",
                             "P" = "D", 
                             "DA" = "D")

cond4$distribution_a <- cond4$distribution

cond4 <- cond4 %>% separate_wider_delim(distribution_a, "/", names = c("dist_1", "dist_2"),
                                        too_few = "align_start")

levels(as.factor(cond4$rate_tissue_loss))

tl <- cond4 %>% subset(rate_tissue_loss == "C"|rate_tissue_loss == "SA,A")


cond4$rate_tissue_loss <- recode(cond4$rate_tissue_loss,
                                 "C" = "SA")

cond4$rate_tissue_loss_a <- cond4$rate_tissue_loss

cond <- cond4 %>% separate_wider_delim(rate_tissue_loss_a, "/", names = c("rate_1", "rate_2"),
                                       too_few = "align_start")


sum_tl <- cond %>% subset(cond_1 == "TL"| cond_2 == "TL") %>%
  #subset(max_diameter > 4) %>%
  group_by(time_point, survey, location_name, scientific_name) %>%
  summarize(n_tl = n())

#lump Orbicellas
orbi <- sum_tl %>% subset(scientific_name == "Orbicella spp."|
                          scientific_name == "Orbicella annularis"|scientific_name == "Orbicella faveolata"|
                          scientific_name == "Orbicella franksi") %>%
  group_by(time_point, survey, location_name) %>%
  summarize(n_tl = sum(n_tl)) %>%
  mutate(scientific_name = "Orbicella spp.")

sum_tl <- sum_tl %>% subset(scientific_name != "Orbicella spp." &
                          scientific_name != "Orbicella annularis" & scientific_name != "Orbicella faveolata" &
                          scientific_name != "Orbicella franksi") %>%
                          bind_rows(orbi)

levels(as.factor(sum_tl$scientific_name))
#sum_tl <- cond %>% subset(cond_1 == "TL"| cond_2 == "TL") %>%
#group_by(time_point, location_name, scientific_name) %>%
#summarize(n_tl = n()) #this makes prev plot A

sum_tl$survey_a <- sum_tl$survey
sum_tl$location_name_a <- sum_tl$location_name
sum_tl$scientific_name_a <- sum_tl$scientific_name

sum_tl <- sum_tl %>% unite("event", c("survey_a", "location_name_a", "scientific_name_a"), sep = "_")

demo <- demo %>% 
  mutate(survey = case_when(
    survey == "October19" | survey == "January20" ~ "November19",
    TRUE ~ as.factor(survey)))

sum_demo <- demo %>%
  group_by(time_point, survey, location_name, scientific_name) %>%
  summarize(total = sum(total))

sum_demo$survey_a <- sum_demo$survey
sum_demo$location_name_a <- sum_demo$location_name
sum_demo$scientific_name_a <- sum_demo$scientific_name

sum_demo <- sum_demo %>% unite("event", c("survey_a", "location_name_a", "scientific_name_a"), sep = "_")


tldf <- sum_demo %>% left_join(sum_tl, by = "event") %>% 
  dplyr::select(-c(time_point.y:scientific_name.y)) %>%
  replace(is.na(.), 0) %>%
  mutate(prev_tl = n_tl/total) %>%
  mutate(n_healthy = total - n_tl) %>%
  subset(!(is.na(prev_tl)))


colnames(tldf) <- gsub("\\.x","",colnames(tldf))


tldf$survey <- recode(tldf$survey, 
                      "January20" = "November19")


tldf <- tldf %>% rename("species" = scientific_name)


write.csv(tldf, "CBC_conditions_curated.csv")


