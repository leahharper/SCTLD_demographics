library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)

# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/Coral Disease")


# Import the survey data
cover <- read.csv("annotations_CBC_19-22.csv")
one <- read.csv("CBCdemo19.csv")
two <- read.csv("CBCdemoMay22.csv")
three <- read.csv("CBCdemoDec22.csv")

split <- strsplit(cover$Date, "/")
split2 <- matrix(unlist(split), ncol=3, byrow=TRUE)
split2 <- as.data.frame(split2)
cover <- cbind(cover, split2)

cover <- cover %>% rename("Year" = "V3", "Month" = "V1") %>%
  select(Name, Date, Habitat, SiteName, Year, Month, Label, Row, Column)

levels(as.factor(cover$SiteName))

cover$SiteName <- recode(cover$SiteName, "CBC 30 N" = "CBC30N",
                         "CBC House Reef" = "House Reef")


cover$Year_a <- cover$Year
cover$Month_a <- cover$Month
cover$SiteName_a <- cover$SiteName

cover <- cover %>% subset(SiteName != "Tobacco Reef") %>% droplevels()

cover <- cover %>% unite("TimePoint", c("Month_a", "Year_a"), sep = "-")

levels(as.factor(cover$TimePoint))

cover$TimePoint <- recode(cover$TimePoint, "1-20" = "10-19")

cover$TimePoint <- recode(cover$TimePoint, "10-19" = "Oct19",
                          "12-22" = "Dec22", 
                          "5-22" = "May22")

cover$TimePoint_a <- cover$TimePoint

cover <- cover %>% unite("Event", c("TimePoint_a", "SiteName_a"), sep = "_")


cover1 <- cover %>% group_by(Event) %>% summarize(npics = (length(unique(Name)))) %>%
  mutate(npoints = npics * 40)

coverx <- left_join(cover, cover1, by = "Event")

cover2 <- coverx %>% 
  group_by(Habitat, SiteName, TimePoint, Label, npoints) %>%
  summarize(n = length(Label)) %>%
  mutate(cover = (n/npoints)*100) %>%
  ungroup() %>%
  complete(Label, nesting(Habitat, SiteName, npoints, TimePoint), 
           fill = list(n = 0, cover = 0)) %>%
  arrange(Habitat, SiteName, npoints, TimePoint, Label) %>%
  select(Habitat, SiteName, npoints, TimePoint, Label, n, cover)

cover2$Label <- as.factor(cover2$Label)
levels(cover2$Label)

cover2$Label <- recode(cover2$Label, "CNAt" = "CNAT",
                       "SMIC" = "SINT",
                       "PSTRI" = "PSTR")


cover <- cover2 %>% mutate('Label_General' = 'x', 'Label_Suscep' = 'x')
cover <- within(cover, Label_General[Label == "THAL"] <- "Thalassia")
cover <- within(cover, Label_General[Label == "AAGA"|Label == "APAL"|Label == "ATEN"|
                                       Label == "Coral"|Label == "MALC"|Label == "MCAV"|
                                       Label == "MCOM"|Label == "OANN"|Label == "OFAV"|
                                       Label == "PAST"|Label == "PPOR"|Label == "PSTR"|
                                       Label == "SINT"|Label == "SSID"|Label == "PCLI"|
                                       Label == "CNAT"|Label == "ACER"|Label == "AAGA"|
                                       Label == "ATEN"|Label == "DLAB"|Label == "MCOM"|
                                       Label == "SRAD"] <- "All Stony Coral")
cover <- within(cover, Label_Suscep[Label == "MCAV"|
                                      Label == "OANN"|Label == "OFAV"|
                                      Label == "PSTR"|Label == "CNAT"|
                                      Label == "SINT"|Label == "SSID"|Label == "PCLI"|
                                      Label == "DLAB"] <- "SCTLD Suscep. Spp.")
cover <- within(cover, Label_General[Label == "ANTI"|Label == "GORGO"|Label == "GVEN"|
                                       Label == "MURO"|Label == "PTER"|Label == "EUNI"|
                                       Label == "BRIA"|Label == "PSEUGORG"|Label == "MURI"]<- "Octocoral")
cover <- within(cover, Label_General[Label == "Hal_spp"|Label == "Macro"|Label == "Turbin"|
                                       Label == "Lob_spp"|Label == "Dicsp"|Label == "BRAN-CALC"] <- "Macroalgae")
cover <- within(cover, Label_General[Label == "Sand"|Label == "LSUB_RUB"] <- "Unconsolidated (Sand/Rubble)")
cover <- within(cover, Label_General[Label == "Tk Tf"|Label == "Tn Tf"] <- "Turf Algae")
cover <- within(cover, Label_General[Label == "Sand"|Label == "LSUB_RUB"] <- "Unconsolidated (Sand/Rubble)")
cover <- within(cover, Label_General[Label == "SpgOth"|Label == "SPONGR"|Label == "SPTU"|
                                       Label == "ENSP"|Label == "CLIONI"|Label == "SPONGB"|
                                       Label == "SPONGV"] <- "Sponge")
cover <- within(cover, Label_General[Label == "PALY"|Label == "Other Inv"] <- "Other Encrusting Invert")

cover <- within(cover, Label_General[Label == "CCA 1"] <- "Hard Substrate")
cover <- within(cover, Label_General[Label == "TAPE"|Label == "Unk"|Label == "SHAD"] <- "Unidentified")

cover <- within(cover, Label_General[Label == "CYAN"|Label == "Cyan red"] <- "Cyanobacteria")


pd <- position_dodge(width = 0.93)
se<-function(x)sqrt(var(x)/length(x))

TimeLevels = c('Oct19','May22','Dec22')
LabLevels = c('Turf Algae','All Stony Coral','Macroalgae','Octocoral')


cover <- cover %>% mutate_at(.vars = vars("TimePoint"), 
                             .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
  mutate_at(.vars = vars("Label_General"),
            .funs = funs(factor(.,levels = LabLevels, ordered = TRUE)))



covgen <- cover %>% group_by(SiteName, Habitat, TimePoint, npoints, Label_General) %>% 
  summarize(sum_gen = sum(n)) %>% 
  mutate(cov_gen = (sum_gen/npoints)*100) %>% 
  mutate(cov_gensq = sqrt(cov_gen)) %>%
  mutate(cov_prop = cov_gen/100) %>%
  subset(!(is.na(Label_General))) 

covgen$TimePoint_a <- covgen$TimePoint
covgen$Label_General_a <- covgen$Label_General

covgen <- covgen %>% unite("GroupEvent", c(Label_General_a, TimePoint_a))

covgroup <- covgen %>% group_by(TimePoint, Label_General) %>%
  summarize(MeanCov = mean(cov_gen), seCov = se(cov_gen), maxCov = max(cov_gen)) %>%
  unite("GroupEvent", c(Label_General,TimePoint))

covgen <- left_join(covgen, covgroup, by = "GroupEvent")


suscep <- cover %>% group_by(SiteName, Habitat, TimePoint, npoints, Label_Suscep) %>% 
  summarize(sum_gen = sum(n)) %>%
  mutate(cov_gen = (sum_gen/npoints)*100) %>% subset(Label_Suscep == "SCTLD Suscep. Spp.")

covsto <- covgen %>% subset(Label_General == "All Stony Coral") %>%
  mutate(cov_gensq = sqrt(cov_gen)) %>%
  mutate(cov_prop = cov_gen/100)


suscep <- suscep %>% 
  mutate(cov_gensq = sqrt(cov_gen)) %>%
  mutate(cov_prop = cov_gen/100) %>%
  unite("GroupEvent", c(Label_Suscep,TimePoint), remove = FALSE)


susgroup <- suscep %>% group_by(TimePoint, Label_Suscep) %>%
  summarize(MeanCov = mean(cov_gen), seCov = se(cov_gen), maxCov = max(cov_gen)) %>%
  unite("GroupEvent", c(Label_Suscep,TimePoint))

suscep <- suscep %>% left_join(susgroup, by = "GroupEvent")

library(betareg)
library(lme4)
library(emmeans)
library(rcompanion)
library(glmmTMB)

mod <- glm(cov_prop ~ TimePoint + SiteName, data = covsto, family = "quasibinomial")
hist(resid(mod))
Anova(mod)
rsquared(mod)

mod2 <- betareg(cov_prop ~ TimePoint + SiteName, data = covsto, link = "logit")
hist(resid(mod2))
Anova(mod2)

mod3 <- lmer(cov_prop ~ TimePoint + (1|SiteName), data = covsto)
hist(resid(mod3))
Anova(mod3)

mod4 <- glmer(cov_prop ~ TimePoint +
                (1|SiteName),
              weights=npoints,
              data=covsto,family="binomial")
hist(resid(mod4))
Anova(mod4)

#mod5 <- glmmTMB(cov_prop ~ TimePoint +
                 # (1|SiteName),
               # data=covsto,beta_family())
#hist(resid(mod4))
#Anova(mod4)

isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}

Labels <- levels(as.factor(covgen2$Label_General))

df <- data.frame()
mod_df <- data.frame()

for(current_Label in Labels) {
  
  Lab_df <- covgen %>% subset(Label_General == current_Label) 
  mod <- glmer(cov_prop ~ TimePoint +
                 (1|SiteName),
               weights=npoints,
               data=Lab_df,family="binomial")
  
  modresult <- as.data.frame(Anova(mod)) %>%
    mutate("Label_General" = current_Label)
  modresult$sig <- isSig(modresult$"Pr(>Chisq)")
  
  mod_df <- mod_df %>%
    bind_rows(modresult)
  
  pairwise <- as.data.frame(pairs(emmeans(mod, ~ TimePoint)))
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = pairwise,
                                   threshold = 0.05)) %>%
    mutate("Label_General" = current_Label)
  
  df <- df %>%
    bind_rows(letters)
}

MainLetters <- df %>% unite("GroupEvent", c(Label_General,Group))


covgen <- covgen %>% 
  left_join(MainLetters, by = "GroupEvent")

covgen1 <- covgen %>% subset(Label_General == "Turf Algae"|
                               Label_General == "All Stony Coral")

covgen2 <- covgen %>% subset(Label_General == "Macroalgae"|
                               Label_General == "Octocoral")

modsus <- glm(cov_prop ~ TimePoint + SiteName, data = suscep, family = "quasibinomial")
modsus <- glmer(cov_prop ~ TimePoint +
                  (1|SiteName),
                weights=npoints,
                data=suscep,family="binomial")
hist(resid(modsus))
Anova(modsus)
rsquared(modsus)

susresult <- as.data.frame(Anova(modsus))
susresult$sig <- isSig(susresult$"Pr(>Chisq)")

suspairwise <- as.data.frame(pairs(emmeans(modsus, ~ TimePoint)))

suspairwise$sig <- isSig(suspairwise$p.value)

susletters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = suspairwise,
                                 threshold = 0.05))

suscep <- suscep %>% left_join(susletters, by = c("TimePoint" = "Group"))
 
suscep <- suscep %>% mutate(TimePoint = fct_relevel(TimePoint,
                                  "Oct19", "May22", "Dec22")) 

library(viridis)
library(cowplot)
library(gridExtra)

levels(as.factor(covgen$SiteName))

sitecolors = c('CBC Central'='red3','CBC Lagoon'='darkorchid','CBC30N'='gold1',
               'Curlew Patch' = 'blue', 'House Reef' ='seagreen4', 
               'South Reef Central' = 'orange', 'SR30N' = 'pink')


labyearp1 <- ggplot() +
  geom_boxplot(data = covgen1, aes(x = TimePoint, y = cov_gen), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = covgen1, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 2, pch = 21,  
              width = 0.25, 
              height = 0) +
  geom_text(data = covgen1, 
            aes(x = TimePoint, y = maxCov, label = Letter), nudge_y = 5) + 
  facet_wrap(~Label_General, scales = "free") +
  scale_y_continuous("Percent Cover") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
labyearp1

labyearp2 <- ggplot() +
  geom_boxplot(data = covgen2, aes(x = TimePoint, y = cov_gen), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = covgen2, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 2, pch = 21,  
              width = 0.25, 
              height = 0) +
  geom_text(data = covgen2, 
            aes(x = TimePoint, y = maxCov, label = Letter), nudge_y = 5) +
  facet_wrap(~Label_General, scales = "free") +
  scale_y_continuous("Percent Cover") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
labyearp2

suscepyearp <- ggplot() +
  geom_boxplot(data = suscep, aes(x = TimePoint, y = cov_gen), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = suscep, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 2, pch = 21,  
              width = 0.25, 
              height = 0) +
  geom_text(data = suscep, 
            aes(x = TimePoint, y = maxCov, label = Letter), nudge_y = 5) +
  facet_wrap(~Label_Suscep, scales = "free") +
  scale_y_continuous("") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
suscepyearp


legendp <- ggplot() +
  geom_boxplot(data = suscep, aes(x = TimePoint, y = cov_gen), alpha = 0) +
  #geom_point(data = covgen, aes(x = TimePoint, y = cov_gen, fill = SiteName)) +
  geom_jitter(data = suscep, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 2, pch = 21,  
              width = 0.25, 
              height = 0) +
  facet_wrap(~Label_Suscep, scales = "free") +
  scale_y_continuous("") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10))
legendp

legend <- cowplot::get_legend(legendp)

plot1 <- cowplot::plot_grid(labyearp1, suscepyearp, rel_widths = c(2/3, 1/3), axis = 'tblr', align = "v")
plot2 <- cowplot::plot_grid(labyearp2, legend, rel_widths = c(2/3, 1/3), axis = 't', align = "v")


tiff("CBCCoverBoxplot.tif",width = 7, height = 6, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol=1, nrow =2)
dev.off()

png("CBCCoverBoxplot.png", width = 7, height = 6, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol=1, nrow =2)
dev.off()




sto <- cover %>% subset(Label_General == "All Stony Coral") %>% droplevels()
levels(as.factor(sto$Label))

sto <- sto %>% mutate(Category = ifelse(Label == "ACER" | Label == "ATEN" |
                                          Label == "PPOR" | Label == "PAST" |
                                          Label == "AAGA"| Label == "SRAD",
                                        "SCTLD Resistant", "x")) %>%
  mutate(Category = ifelse(Label == "DLAB"|Label == "PSTR"|Label == "CNAT",
                           "Highly SCTLD Susceptible", Category)) %>%
  mutate(Category = ifelse(Label == "MCAV"|Label == "SSID"|Label == "OANN"|
                             Label == "OFAV"|Label == "SINT",
                           "Int. SCTLD Susceptible", Category)) %>%
  subset(Category != "x") %>%
  subset(Label != "ACER" & Label != "CNAT" & Label != "DLAB" & Label != "SRAD" & Label != "SINT")

category <- sto %>% group_by(SiteName, Habitat, TimePoint, npoints, Category) %>% 
  summarize(sum_cat = sum(n)) %>%
  mutate(cov_cat = (sum_cat/npoints)*100) %>%
  mutate(cov_prop = cov_cat/100)



Cats <- levels(as.factor(category$Category))

df <- data.frame()
mod_df <- data.frame()

for(current_Cat in Cats) {
  
  Cat_df <- category %>% subset(Category == current_Cat) 
  mod <- glmer(cov_prop ~ TimePoint +
                 (1|SiteName),
               weights=npoints,
               data=Cat_df,family="binomial")  
  
  modresult <- as.data.frame(Anova(mod)) %>%
    mutate("Category" = current_Cat)
  modresult$sig <- isSig(modresult$"Pr(>Chisq)")
  
  mod_df <- mod_df %>%
    bind_rows(modresult)
  
  pairwise <- as.data.frame(pairs(emmeans(mod, ~ TimePoint)))
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = pairwise,
                                   threshold = 0.05)) %>%
    mutate("Category" = current_Cat)
  
  df <- df %>%
    bind_rows(letters)
}

sto2 <- sto %>% rename("Species" = Label) %>%
  mutate(cov_prop = cover/100) %>%
  subset(Species != "ACER" & Species != "CNAT" &
           Species != "DLAB" & Species != "MALC" &
           Species != "MCOM" & Species != "SINT" & 
           Species != "SRAD") %>% droplevels()

SpecList <- levels(as.factor(sto2$Species))



df <- data.frame()
mod_df <- data.frame()

for(current_Spec in SpecList) {
  
  Sp_df <- sto2 %>% subset(Species == current_Spec) 
  mod <- glmer(cov_prop ~ TimePoint +
             (1|SiteName),
              weights=npoints,
              data=Sp_df,family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
              #found this online as a suggestion to fix convergence issues,
              #not totally sure if appropriate but it works
  
  modresult <- as.data.frame(Anova(mod)) %>%
    mutate("Species" = current_Spec)
  modresult$sig <- isSig(modresult$"Pr(>Chisq)")
  
  mod_df <- mod_df %>%
    bind_rows(modresult)
  
  pairwise <- as.data.frame(pairs(emmeans(mod, ~ TimePoint)))
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = pairwise,
                                   threshold = 0.05)) %>%
    mutate("Species" = current_Spec)
  
  df <- df %>%
    bind_rows(letters)
}


speclet <- df %>% unite("Event", c("Species", "Group"))

stomeans <- sto2 %>% group_by(Species, TimePoint) %>%
  summarize(MeanCov = mean(cover), seCov = se(cover), maxCov = max(cover)) %>%
  unite("Event", c(Species,TimePoint))

sto3 <- sto2 %>% unite("Event", c("Species", "TimePoint"), remove = FALSE) %>%
  left_join(speclet, by = "Event") %>%
  left_join(stomeans, by = "Event")


library(ggh4x)

targetp1 <- ggplot() +
  geom_boxplot(data = sto3, aes(x = TimePoint, y = cover), fill = "gray80") +
  geom_jitter(data = sto3, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21,  
              width = 0.15, 
              height = 0) +
  geom_text(data = sto3, 
            aes(x = TimePoint, y = maxCov, label = Letter), nudge_y = 0.7) +
  facet_nested_wrap(~Category + Species, scales = "free",
                    ncol = 5, labeller = labeller(Category = label_wrap_gen(width = 16))) +
  scale_y_continuous("Percent Cover") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10))

tiff("CBCTargetsBoxplot.tif",width = 10, height = 8, units = "in", res = 300)
targetp1
dev.off()

png("CBCTargetsBoxplot.png", width = 10, height = 8, units = "in", res = 300)
targetp1
dev.off()

############################
#DEMO
############################
demo <- rbind.fill(one, two, three) %>%
  select(-c(transect:coral_demographics_notes))


demo <- demo %>%
  replace(is.na(.), 0) %>%
  mutate(total_count = rowSums(.[8:14])) %>%
  mutate(total_small = rowSums(.[10:12])) %>%
  mutate(total_lg = rowSums(.[13:14])) %>%
  mutate(non_juv = rowSums(.[10:14]))

levels(as.factor(demo$location_name))

demo <- demo %>% subset(location_name != "")

sitecolors = c('CBC Central'='red3','CBC Lagoon'='darkorchid','CBC30N'='gold1',
               'Curlew Patch' = 'blue', 'House Reef' ='seagreen4', 
               'South Reef Central' = 'orange', 'SR30N' = 'pink')

demo$location_name <- recode(demo$location_name, 
                             "CBC 30 C" = "CBC Central",
                             "South Reef 30 C" = "South Reef Central",
                             "CBC 30 North" = "CBC30N",
                             "South Reef 30 North" = "SR30N",
                             "CBC House Reef" = "House Reef",
                             "CBC Reef Centra" = "CBC Central",
                             "Curlew Patch Reef" = "Curlew Patch",
                             "CBC Lagoon Reef" = "CBC Lagoon", 
                             "CBC Reef Central" = "CBC Central")

demo <- demo %>% subset(location_name != "Tobacco Reef")

demo$sample_collection_year <- as.factor(demo$sample_collection_year)
levels(demo$sample_collection_year)

demo <- demo %>% subset(sample_collection_year != 0) %>% droplevels()

demo <- demo %>% mutate(time_point = "x") %>%
  mutate(time_point = case_when(
    sample_collection_year == "2019"|sample_collection_year == "2020" ~ "Oct19",
    TRUE ~ as.character(time_point))) %>% 
  mutate(time_point = case_when(
    sample_collection_year == "2022" & sample_collection_month == "5" ~ "May22",
    TRUE ~ as.character(time_point))) %>% 
  mutate(time_point = case_when(
    sample_collection_year == "2022" & sample_collection_month == "12" ~ "Dec22",
    TRUE ~ as.character(time_point)))

demo$scientific_name <- as.factor(demo$scientific_name)
levels(demo$scientific_name)

c <- demo %>% subset(scientific_name == "MYCE")

#remove unknown (LA?), lumped ORBIs and AGARs (uncommon, mainly a Leah thing), and MCOM

demo <- demo %>% subset(scientific_name != "LA(?)" &
                          scientific_name != "ORBI" &
                          scientific_name != "AGAR" &
                          scientific_name != "MCOM" &
                          scientific_name != "MYCE") %>% droplevels()
demo$scientific_name <- recode(demo$scientific_name, "PSTRI" = "PSTR",
                               "DSTR" = "PSTR",
                               "ISRI" = "IRIG")


########################
#NMDS
########################
library(vegan)

demo$location_name_a <- demo$location_name
demo$time_point_a <- demo$time_point


demo <- demo %>%
  unite("event", c("location_name_a", "time_point_a"), sep = "_")

cast <- demo %>% group_by(location_name, time_point, event, scientific_name) %>%
  summarize(total = sum(total_count)) %>%
  pivot_wider(id_cols = c(location_name, time_point, event), names_from = scientific_name, 
              values_from = total)

mat = cast[,3:ncol(cast)]
mat <- mat %>% remove_rownames %>% column_to_rownames(var="event")
mat[is.na(mat)] <- 0
mat[mat > 0] <- 1
mat <- as.matrix(mat)
#mat <- sqrt(mat)
set.seed(123456)


NMDS1 <-
  metaMDS(mat,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

goodness(NMDS1)
stressplot(NMDS1)
plot(NMDS1, type = "t")

#set up grouping variables
data.scores1 = as.data.frame(scores(NMDS1)$sites)

data.scores1$location_name = cast$location_name
data.scores1$time_point = cast$time_point
data.scores1$event = cast$event

data.scores1 <- data.scores1 %>% 
  mutate(time_point = fct_relevel(time_point,
                                  "Oct19", "May22", "Dec22")) 



species.scores1 <- as.data.frame(scores(NMDS1, "species")) 
species.scores1$species <- rownames(species.scores1)  # create a column of species, from the rownames of species.scores

uno <- data.scores1[data.scores1$time_point == "Oct19", ][chull(data.scores1[data.scores1$time_point == 
                                                                               "Oct19", c("NMDS1", "NMDS2")]), ]
dos <- data.scores1[data.scores1$time_point == "May22", ][chull(data.scores1[data.scores1$time_point == 
                                                                               "May22", c("NMDS1", "NMDS2")]), ]
tres <- data.scores1[data.scores1$time_point == "Dec22", ][chull(data.scores1[data.scores1$time_point == 
                                                                                "Dec22", c("NMDS1", "NMDS2")]), ]

hull.data1 <- rbind(uno, dos, tres) %>% 
  mutate(time_point = fct_relevel(time_point,
                                  "Oct19", "May22", "Dec22")) 


hull.data1

library(ggrepel)
library(viridis)

tiff("allNMDS.tif",width = 7, height = 7, units = "in", res = 400)
ggplot() + 
  geom_polygon(data=hull.data1,aes(x=NMDS1,y=NMDS2,fill=time_point,group=time_point),alpha=0.30, color = "black") + 
  geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2,colour = time_point, fill = time_point, shape = location_name),size=4, alpha = 0.8) +
  geom_point(data=species.scores1,aes(x=NMDS1,y=NMDS2), pch = 20, color = "black", size=2) +
  geom_text_repel(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),size = 2, min.segment.length = 0.5, max.overlaps = 18) + 
  #geom_text_repel(data=data.scores1,aes(x=NMDS1,y=NMDS2,label=location_name),size = 4, colour = "black", min.segment.length = 0.5, max.overlaps = 18) + 
  scale_fill_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Site",values=c("CBC Central" = 21, "CBC Lagoon" = 22, "House Reef" = 23, "CBC30N" = 24,
                                     "South Reef Central" = 25, "Curlew Patch" = 3, "SR30N" = 9)) +
  #xlim(-1,2) +  
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()


#large corals

lgcast <- demo %>% group_by(location_name, time_point, event, scientific_name) %>%
  summarize(total = sum(total_lg)) %>%
  pivot_wider(id_cols = c(location_name, time_point, event), names_from = scientific_name, 
              values_from = total)

lgmat = lgcast[,3:ncol(lgcast)]
lgmat <- lgmat %>% remove_rownames %>% column_to_rownames(var="event")
lgmat[is.na(lgmat)] <- 0
lgmat[lgmat > 0] <- 1
lgmat <- as.matrix(lgmat)
#lgmat <- sqrt(lgmat)
set.seed(123456)



lgNMDS1 <-
  metaMDS(lgmat,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

goodness(lgNMDS1)
stressplot(lgNMDS1)
plot(lgNMDS1, type = "t")

#set up grouping variables
data.scores1 = as.data.frame(scores(lgNMDS1)$sites)

data.scores1$location_name = lgcast$location_name
data.scores1$time_point = lgcast$time_point
data.scores1$event = lgcast$event

data.scores1 <- data.scores1 %>% 
  mutate(time_point = fct_relevel(time_point,
                                  "Oct19", "May22", "Dec22")) 



species.scores1 <- as.data.frame(scores(lgNMDS1, "species")) 
species.scores1$species <- rownames(species.scores1)  # create a column of species, from the rownames of species.scores

uno <- data.scores1[data.scores1$time_point == "Oct19", ][chull(data.scores1[data.scores1$time_point == 
                                                                               "Oct19", c("NMDS1", "NMDS2")]), ]
dos <- data.scores1[data.scores1$time_point == "May22", ][chull(data.scores1[data.scores1$time_point == 
                                                                               "May22", c("NMDS1", "NMDS2")]), ]
tres <- data.scores1[data.scores1$time_point == "Dec22", ][chull(data.scores1[data.scores1$time_point == 
                                                                                "Dec22", c("NMDS1", "NMDS2")]), ]

hull.data1 <- rbind(uno, dos, tres) %>% 
  mutate(time_point = fct_relevel(time_point,
                                  "Oct19", "May22", "Dec22")) 


hull.data1


lgNMDS <- ggplot() + 
  geom_polygon(data=hull.data1,aes(x=NMDS1,y=NMDS2,fill=time_point,group=time_point),alpha=0.30, color = "black") + 
  geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2,colour = time_point, fill = time_point, shape = location_name),size=4, alpha = 0.8) +
  geom_point(data=species.scores1,aes(x=NMDS1,y=NMDS2), pch = 20, color = "black", size=2) +
  geom_text_repel(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),size = 2, min.segment.length = 0.5, max.overlaps = 18) + 
  #geom_text_repel(data=data.scores1,aes(x=NMDS1,y=NMDS2,label=location_name),size = 4, colour = "black", min.segment.length = 0.5, max.overlaps = 18) + 
  ggtitle("Large Corals (>40cm)") +
  scale_fill_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Site",values=c("CBC Central" = 21, "CBC Lagoon" = 22, "House Reef" = 23, "CBC30N" = 24,
                                     "South Reef Central" = 25, "Curlew Patch" = 3, "SR30N" = 9)) +
  # xlim(-1,2) +  
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5))

tiff("lgNMDS.tif",width = 7, height = 7, units = "in", res = 400)
lgNMDS
dev.off()

#small corals


smcast <- demo %>% group_by(location_name, time_point, event, scientific_name) %>%
  summarize(total = sum(total_small)) %>%
  pivot_wider(id_cols = c(location_name, time_point, event), names_from = scientific_name, 
              values_from = total)

smmat = smcast[,3:ncol(smcast)]
smmat <- smmat %>% remove_rownames %>% column_to_rownames(var="event")
smmat[is.na(smmat)] <- 0
smmat[smmat > 0] <- 1
smmat <- as.matrix(smmat)
#smmat <- sqrt(smmat)
set.seed(123456)


smNMDS1 <-
  metaMDS(smmat,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

goodness(smNMDS1)
stressplot(smNMDS1)
plot(smNMDS1, type = "t")

#set up grouping variables
data.scores1 = as.data.frame(scores(smNMDS1)$sites)

data.scores1$location_name = smcast$location_name
data.scores1$time_point = smcast$time_point
data.scores1$event = smcast$event

data.scores1 <- data.scores1 %>% 
  mutate(time_point = fct_relevel(time_point,
                                  "Oct19", "May22", "Dec22")) 



species.scores1 <- as.data.frame(scores(smNMDS1, "species")) 
species.scores1$species <- rownames(species.scores1)  # create a column of species, from the rownames of species.scores

uno <- data.scores1[data.scores1$time_point == "Oct19", ][chull(data.scores1[data.scores1$time_point == 
                                                                               "Oct19", c("NMDS1", "NMDS2")]), ]
dos <- data.scores1[data.scores1$time_point == "May22", ][chull(data.scores1[data.scores1$time_point == 
                                                                               "May22", c("NMDS1", "NMDS2")]), ]
tres <- data.scores1[data.scores1$time_point == "Dec22", ][chull(data.scores1[data.scores1$time_point == 
                                                                                "Dec22", c("NMDS1", "NMDS2")]), ]

hull.data1 <- rbind(uno, dos, tres) %>% 
  mutate(time_point = fct_relevel(time_point,
                                  "Oct19", "May22", "Dec22")) 


hull.data1


smNMDS <- ggplot() + 
  geom_polygon(data=hull.data1,aes(x=NMDS1,y=NMDS2,fill=time_point,group=time_point),alpha=0.30, color = "black") + 
  geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2,colour = time_point, fill = time_point, shape = location_name),size=4, alpha = 0.8) +
  geom_point(data=species.scores1,aes(x=NMDS1,y=NMDS2), pch = 20, color = "black", size=2) +
  geom_text_repel(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),size = 2, min.segment.length = 0.5, max.overlaps = 18) + 
  #geom_text_repel(data=data.scores1,aes(x=NMDS1,y=NMDS2,label=location_name),size = 4, colour = "black", min.segment.length = 0.5, max.overlaps = 18) + 
  ggtitle("Small Corals (5-40cm)") +
  scale_fill_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Site",values=c("CBC Central" = 21, "CBC Lagoon" = 22, "House Reef" = 23, "CBC30N" = 24,
                                     "South Reef Central" = 25, "Curlew Patch" = 3, "SR30N" = 9)) +
  # xlim(-1,2) +  
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5))

tiff("smNMDS.tif",width = 7, height = 7, units = "in", res = 400)
smNMDS
dev.off()



njcast <- demo %>% group_by(location_name, time_point, event, scientific_name) %>%
  summarize(total_nj = sum(non_juv)) %>%
  pivot_wider(id_cols = c(location_name, time_point, event), names_from = scientific_name, 
              values_from = total_nj) %>%
  mutate_all(~replace(., is.na(.), 0))

#PERMANOVA using adonis function
adonis2(formula = lgmat ~ location_name + time_point, data = lgcast, permutations = 10000)

adonis2(formula = smmat ~ location_name + time_point, data = smcast, permutations = 10000)


njlong <- melt(njcast, id.vars =c(
  "location_name","time_point","event")) %>%
  rename("scientific_name" = variable, "count_all" = value) %>%
  mutate(time_point = fct_relevel(time_point,
                                  "Oct19", "May22", "Dec22")) 


TimeLevels = c('Oct19','May22','Dec22')

njlong <- njlong %>% mutate_at(.vars = vars("time_point"), 
                               .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))



suscep1 <- njlong %>% subset(scientific_name == "MMEA" |
                               scientific_name == "EFAS" |
                               scientific_name == "DSTO" |
                               scientific_name == "DCYL") %>%
  group_by(location_name, time_point, scientific_name) %>%
  summarize(density = (sum(count_all))/30)

susmeans1 <- suscep1 %>% group_by(time_point, scientific_name) %>%
  summarize(MeanDens = mean(density), seDens = se(density), maxDens = max(density)) %>%
  unite("GroupEvent", c(scientific_name,time_point))

suscep1 <- suscep1 %>% unite("GroupEvent", c(scientific_name, time_point), remove = FALSE) %>%
  left_join(susmeans1, by = "GroupEvent")

library(MASS)

mcav <- njlong %>% subset(scientific_name == "MCAV")

mod <- glm(count_all ~ time_point + location_name, data = mcav, family = "poisson")
hist(resid(mod))
Anova(mod)
rsquared(mod)

mod2 <- glm.nb(count_all ~ time_point + location_name, 
               data = mcav, link = "log", maxit = 500)
hist(resid(mod2))
Anova(mod2)
rsquared(mod2)

mod3 <- glmer(count_all ~ time_point +
                (1|location_name),
              data=mcav,family="poisson")
hist(resid(mod3))
Anova(mod3)
rsquared(mod3)

mod4 <- lmer(count_all ~ time_point + (1|location_name), data = mcav)
hist(resid(mod4))
Anova(mod4)

AIC(mod3)
AIC(mod2)
AIC(mod)
AIC(mod4)

dens <- njlong %>% rename("Species" = scientific_name) %>%
  subset(Species == "DSTO" | Species == "EFAS" |
           Species == "MMEA" | Species == "DCYL" |
           Species == "PSTR" | Species == "MCAV" |
           Species == "DLAB" | Species == "SSID") %>% droplevels()

SpecList <- levels(as.factor(dens$Species))

df <- data.frame()
mod_df <- data.frame()

for(current_Spec in SpecList) {
  
  Sp_df <- dens %>% subset(Species == current_Spec) 
  # <- lmer(count_all ~ time_point + (1|location_name), data = Sp_df)  
  mod <- glm(count_all ~ time_point + location_name, data = Sp_df, family = "poisson")
  
  modresult <- as.data.frame(Anova(mod)) %>%
    mutate("Species" = current_Spec)
  modresult$sig <- isSig(modresult$"Pr(>Chisq)")
  
  mod_df <- mod_df %>%
    bind_rows(modresult)
  
  pairwise <- as.data.frame(pairs(emmeans(mod, ~ time_point)))
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = pairwise,
                                   threshold = 0.05)) %>%
    mutate("Species" = current_Spec)
  
  df <- df %>%
    bind_rows(letters)
}

susletters <- df %>%
  unite("GroupEvent", c(Species, Group))

suscep1 <- suscep1 %>% left_join(susletters, by = "GroupEvent")


suscep1p <- ggplot() +
  #geom_bar(stat = "identity", color = "black", position = pd) +
  geom_violin(data = suscep1, aes(x = time_point, y = density), fill = "gray80") +
  geom_jitter(data = suscep1, aes(x = time_point, y = density, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0) +
  geom_text(data = suscep1, 
            aes(x = time_point, y = maxDens, label = Letter), nudge_y = 0.03) + 
  facet_wrap(~scientific_name, scales = "free", ncol = 4) +
  scale_y_continuous("Mean Density/m2") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10))

png("dens_highsuscep.png", width = 8, height = 4, units = "in", res = 400)
suscep1p
dev.off()

suscep2 <- njlong %>% subset(scientific_name == "MCAV" |
                               scientific_name == "PSTR" |
                               scientific_name == "DLAB" |
                               scientific_name == "SSID") %>%
  group_by(location_name, time_point, scientific_name) %>%
  summarize(density = (sum(count_all))/30) 

susmeans2 <- suscep2 %>% group_by(time_point, scientific_name) %>%
  summarize(MeanDens = mean(density), seDens = se(density), maxDens = max(density)) %>%
  unite("GroupEvent", c(scientific_name,time_point))

suscep2 <- suscep2 %>% unite("GroupEvent", c(scientific_name, time_point), remove = FALSE) %>%
  left_join(susmeans2, by = "GroupEvent") %>%
  left_join(susletters, by = "GroupEvent")


suscep2p <- ggplot() +
  geom_violin(data = suscep2, aes(x = time_point, y = density), fill = "gray80") +
  geom_jitter(data = suscep2, aes(x = time_point, y = density, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0) +
  geom_text(data = suscep2, 
            aes(x = time_point, y = maxDens, label = Letter), nudge_y = 0.07) + 
  facet_wrap(~scientific_name, scales = "free", ncol = 4) +
  scale_y_continuous("Mean Density/m2") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10))
        #legend.position = "none")



tiff("dens_modsuscep.tif",width = 8, height = 4, units = "in", res = 400)
suscep2p
dev.off()

png("dens_modsuscep.png", width = 8, height = 4, units = "in", res = 400)
plot(suscep2p)
dev.off()



#################
#Conditions
#################
cond1 <- read.csv("CBCcond19.csv") %>%
  rename("scientific_name" = "taxon_id", "rate_tissue_loss" = "rate_tl",
         "coral_condition_notes" = "sample_collection_notes")
cond2 <- read.csv("CBCcondMay22.csv")
cond3 <- read.csv("CBCcondDec22.csv")

cond <- rbind.fill(cond1, cond2, cond3)

levels(as.factor(cond$location_name))

cond$location_name <- recode(cond$location_name,
                             "CBC 30 C" = "CBC Central",
                             "CBC 30 N" = "CBC30N",
                             "CBC 30 North" = "CBC30N",
                             "CBC House Reef" = "House Reef",
                             "CBC Lagoon Reef" = "CBC Lagoon",
                             "CBC Reef Central" = "CBC Central",
                             "South Reef 30 C" = "South Reef Central",
                             "South Reef 30 Central" = "South Reef Central",
                             "South Reef 30 North" = "SR30N")

cond <- cond %>% mutate(time_point = "x") %>%
  mutate(time_point = case_when(
    sample_collection_year == "2019"|sample_collection_year == "2020" ~ "Oct19",
    TRUE ~ as.character(time_point))) %>% 
  mutate(time_point = case_when(
    sample_collection_year == "2022" & sample_collection_month == "5" ~ "May22",
    TRUE ~ as.character(time_point))) %>% 
  mutate(time_point = case_when(
    sample_collection_year == "2022" & sample_collection_month == "12" ~ "Dec22",
    TRUE ~ as.character(time_point)))

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

sum_tl <- cond %>% subset(cond_1 == "TL" &  perc_1 != "A"| cond_2 == "TL" & perc_2 != "A") %>%
  group_by(time_point, location_name, scientific_name) %>%
  summarize(n_tl = n())

A <- cond %>% subset(perc_1 == "A")

levels(as.factor(cond$perc_1))

#sum_tl <- cond %>% subset(cond_1 == "TL"| cond_2 == "TL") %>%
  #group_by(time_point, location_name, scientific_name) %>%
  #summarize(n_tl = n()) #this makes prev plot A

sum_tl$time_point_a <- sum_tl$time_point
sum_tl$location_name_a <- sum_tl$location_name
sum_tl$scientific_name_a <- sum_tl$scientific_name

sum_tl <- sum_tl %>% unite("Event", c("time_point_a", "location_name_a", "scientific_name_a"), sep = "_")


sum_demo <- demo %>%
  group_by(time_point, location_name, scientific_name) %>%
  summarize(total_nj = sum(non_juv))

sum_demo$time_point_a <- sum_demo$time_point
sum_demo$location_name_a <- sum_demo$location_name
sum_demo$scientific_name_a <- sum_demo$scientific_name

sum_demo <- sum_demo %>% unite("Event", c("time_point_a", "location_name_a", "scientific_name_a"), sep = "_")

tldf <- sum_demo %>% left_join(sum_tl, by = "Event") %>% 
  dplyr::select(-c(time_point.y:scientific_name.y)) %>%
  replace(is.na(.), 0) %>%
  mutate(prev_tl = n_tl/total_nj)

tldf <- na.omit(tldf)

colnames(tldf) <- gsub("\\.x","",colnames(tldf))

tldf <- tldf %>% mutate(scaled_prevtl = prev_tl)

mod <- glm(prev_tl ~ time_point*location_name, data = tldf, family = "quasibinomial")
hist(resid(mod))
Anova(mod)
rsquared(mod)




mod3 <- lmer(prev_tl ~ time_point + (1|location_name), data = tldf)
hist(resid(mod3))
Anova(mod3)


tldf$scaled_prevtl <- scale(tldf$scaled_prevtl)

mod4 <- glmer(prev_tl ~ time_point +
                (1|location_name),
              weights=total_nj,
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=tldf,family="binomial")
hist(resid(mod4))
Anova(mod4)

###########
#most species for supplement (those that appear every year and/or aren't high suscep)
TimeLevels = c('Oct19','May22','Dec22')

tldf <- tldf %>% mutate_at(.vars = vars("time_point"), 
                                         .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


sup <- tldf %>% subset(scientific_name != "ALAM" & scientific_name != "CNAT" &
                         scientific_name != "IRIG" & scientific_name != "ISIN" &
                         scientific_name != "MALI" & scientific_name != "OFRA" &
                         scientific_name != "PCLI" & scientific_name != "ACER")

prevp_spp_all <- ggplot() +
  geom_violin(data = sup, aes(x = time_point, y = prev_tl), fill = "gray80") +
  geom_jitter(data = sup, aes(x = time_point, y = prev_tl, fill = location_name), size = 2.5, pch = 21,  
              width = 0.1, 
              height = 0) +
  facet_wrap(~scientific_name, ncol = 6) +
  scale_y_continuous("Tissue Loss Prevalence") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))


tiff("prevplotall.tif", width = 10, height = 8, units = "in", res = 400)
prevp_spp_all
dev.off()

png("prevplotall.png", width = 10, height = 8, units = "in", res = 400)
prevp_spp_all
dev.off()

##########
#target species

target_prev <- tldf %>% subset(
  scientific_name == "SSID"|
    scientific_name == "PSTR"|scientific_name == "MCAV"|
    scientific_name == "OANN"|scientific_name == "OFAV") %>%
  droplevels() %>% unite("GroupEvent", c(scientific_name, time_point), remove = FALSE)

prev_group <- target_prev %>% group_by(time_point, scientific_name) %>%
  summarize(MeanPrev = mean(prev_tl), sePrev = se(prev_tl), maxPrev = max(prev_tl)) %>%
  unite("GroupEvent", c(scientific_name,time_point))

target_prev <- left_join(target_prev, prev_group, by = "GroupEvent")


TimeLevels = c('Oct19','May22','Dec22')

target_prev <- target_prev %>% mutate_at(.vars = vars("time_point"), 
                                         .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


SpecList <- levels(as.factor(target_prev$scientific_name))

df <- data.frame()
mod_df <- data.frame()

for(current_Specie in SpecList) {
  
  Spec_df <- target_prev %>% subset(scientific_name == current_Specie) 
  
  mod <- lmer(prev_tl ~ time_point + (1|location_name), data = Spec_df)
  
  modresult <- as.data.frame(Anova(mod)) %>%
    mutate("Species" = current_Specie)
  modresult$sig <- isSig(modresult$"Pr(>Chisq)")
  
  mod_df <- mod_df %>%
    bind_rows(modresult)
  
  pairwise <- as.data.frame(pairs(emmeans(mod, ~ time_point)))
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = pairwise,
                                   threshold = 0.05)) %>%
    mutate("Species" = current_Specie)
  
  df <- df %>%
    bind_rows(letters)
}

PrevLetters <- df %>%
  unite("GroupEvent", c(Species,Group))

target_prev <- left_join(target_prev, PrevLetters, by = "GroupEvent")



prevp_spp <- ggplot() +
  geom_violin(data = target_prev, aes(x = time_point, y = prev_tl), fill = "gray80") +
  geom_jitter(data = target_prev, aes(x = time_point, y = prev_tl, fill = location_name), size = 2.5, pch = 21,  
              width = 0.1, 
              height = 0) +
  geom_text(data = target_prev, 
            aes(x = time_point, y = maxPrev, label = Letter), nudge_y = 0.05) + 
  geom_vline(xintercept=c(3.6), linetype="dotted") +
  facet_wrap(~scientific_name, ncol = 5) +
  scale_y_continuous("Tissue Loss Prevalence") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))


#tiff("prevplot.tif", width = 8, height = 4, units = "in", res = 400)
prevp_spp
#dev.off()



#all stony coral

sto_tl <- sum_tl %>% group_by(time_point, location_name) %>%
  summarize(n_tl = sum(n_tl))

sto_tl$time_point_a <- sto_tl$time_point
sto_tl$location_name_a <- sto_tl$location_name

sto_tl <- sto_tl %>% unite("Event", "time_point_a", "location_name_a", sep = "_")

sto_demo <- demo %>%
  group_by(time_point, location_name) %>%
  summarize(total_nj = sum(non_juv))

sto_demo$time_point_a <- sto_demo$time_point
sto_demo$location_name_a <- sto_demo$location_name

sto_demo <- sto_demo %>% unite("Event", "time_point_a", "location_name_a", sep = "_")


prev_sto <- sto_demo %>% left_join(sto_tl, by = "Event") %>% 
  dplyr::select(-c(time_point.y:location_name.y)) %>%
  replace(is.na(.), 0) %>%
  mutate(prev_tl = n_tl/total_nj) %>%
  mutate(label = "All Stony Coral")

colnames(prev_sto) <- gsub("\\.x","",colnames(prev_sto))

sto_groups <- prev_sto %>% group_by(time_point) %>%
  summarize(MeanPrev = mean(prev_tl), sePrev = se(prev_tl), maxPrev = max(prev_tl))

prev_sto <- left_join(prev_sto, sto_groups, by = "time_point")


TimeLevels = c('Oct19','May22','Dec22')

prev_sto <- prev_sto %>% mutate_at(.vars = vars("time_point"), 
                                   .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))



stomod <- lmer(prev_tl ~ time_point + (1|location_name), data = prev_sto)
hist(resid(stomod))
Anova(stomod)
rsquared(stomod)

#this model gave weird results, can't be appropriate
#stomod2 <- glm(n_tl ~ time_point,
                # weights=total_nj,
                # data=prev_sto,family="poisson") 
#hist(resid(stomod2))
#Anova(stomod2)


pairwise <- as.data.frame(pairs(emmeans(stomod, ~ time_point)))

pairwise$sig <- isSig(pairwise$p.value)

letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))

prev_sto <- left_join(prev_sto, letters, by = c("time_point" = "Group"))

TimeLevels = c('Oct19','May22','Dec22')

prev_sto <- prev_sto %>% mutate_at(.vars = vars("time_point"), 
                                   .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


prevp_sto <- ggplot() +
  geom_violin(data = prev_sto, aes(x = time_point, y = prev_tl), fill = "gray80") +
  geom_jitter(data = prev_sto, aes(x = time_point, y = prev_tl, fill = location_name), size = 2.5, pch = 21,  
              width = 0.1, 
              height = 0) +
  geom_text(data = prev_sto, 
            aes(x = time_point, y = maxPrev, label = Letter), nudge_y = 0.05) + 
  facet_grid(~label, scales = 'free') +
  scale_y_continuous("", limits = c(0,1.05)) +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

prevp_sto


#SCTLD suscep stony coral

sus_tl <- sum_tl %>% subset(scientific_name == "MCAV"|
                              scientific_name == "OANN"|scientific_name == "OFAV"|
                              scientific_name == "PSTR"|scientific_name == "CNAT"|
                              scientific_name == "SINT"|scientific_name == "SSID"|scientific_name == "PCLI"|
                              scientific_name == "DLAB"|scientific_name == "DSTO"|
                              scientific_name == "EFAS"|scientific_name == "MMEA"|scientific_name == "DCYL") %>% droplevels() %>%
  group_by(time_point, location_name) %>%
  summarize(n_tl = sum(n_tl))

sus_tl$time_point_a <- sus_tl$time_point
sus_tl$location_name_a <- sus_tl$location_name

sus_tl <- sus_tl %>% unite("Event", "time_point_a", "location_name_a", sep = "_")

sus_demo <- demo %>% subset(scientific_name == "MCAV"|
                              scientific_name == "OANN"|scientific_name == "OFAV"|
                              scientific_name == "PSTR"|scientific_name == "CNAT"|
                              scientific_name == "SINT"|scientific_name == "SSID"|scientific_name == "PCLI"|
                              scientific_name == "DLAB"|scientific_name == "DSTO"|
                              scientific_name == "EFAS"|scientific_name == "MMEA"|scientific_name == "DCYL") %>%
  group_by(time_point, location_name) %>%
  summarize(total_nj = sum(non_juv))

sus_demo$time_point_a <- sus_demo$time_point
sus_demo$location_name_a <- sus_demo$location_name

sus_demo <- sus_demo %>% unite("Event", "time_point_a", "location_name_a", sep = "_")


prev_sus <- sus_demo %>% left_join(sus_tl, by = "Event") %>% 
  dplyr::select(-c(time_point.y:location_name.y)) %>%
  replace(is.na(.), 0) %>%
  mutate(prev_tl = n_tl/total_nj) %>%
  mutate(label = "SCTLD Suscep. Spp.")

colnames(prev_sus) <- gsub("\\.x","",colnames(prev_sus))

sus_groups <- prev_sus %>% group_by(time_point) %>%
  summarize(MeanPrev = mean(prev_tl), sePrev = se(prev_tl), maxPrev = max(prev_tl))

prev_sus <- left_join(prev_sus, sus_groups, by = "time_point")

TimeLevels = c('Oct19','May22','Dec22')

prev_sus <- prev_sus %>% mutate_at(.vars = vars("time_point"), 
                                   .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


susmod <- lmer(prev_tl ~ time_point + (1|location_name), data = prev_sus)
hist(resid(susmod))
Anova(susmod)
rsquared(susmod)

#susmod2 <- glmer(prev_tl ~ time_point +
                  # (1|location_name),
                # weights=total_nj,
                # data=prev_sus,family="binomial")
#hist(resid(susmod2))
#Anova(susmod2)

pairwise <- as.data.frame(pairs(emmeans(susmod, ~ time_point)))

pairwise$sig <- isSig(pairwise$p.value)

letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))

prev_sus <- left_join(prev_sus, letters, by = c("time_point" = "Group"))


TimeLevels = c('Oct19','May22','Dec22')

prev_sus <- prev_sus %>% mutate_at(.vars = vars("time_point"), 
                                   .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


prevp_sus <- ggplot() +
  geom_violin(data = prev_sus, aes(x = time_point, y = prev_tl), fill = "gray80") +
  geom_jitter(data = prev_sus, aes(x = time_point, y = prev_tl, fill = location_name), size = 2.5, pch = 21,  
              width = 0.1, 
              height = 0) +
  geom_text(data = prev_sus,
            aes(x = time_point, y = maxPrev, label = Letter), nudge_y = 0.05) + 
  facet_grid(~label, scales = 'free') +
  scale_y_continuous("", limits = c(0,1.05)) +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
prevp_sus


plot3 <- grid.arrange(prevp_sus, prevp_sto, legend, ncol=3, nrow =1)

plot4 <- cowplot::plot_grid(prevp_spp, plot3, rel_widths = c(4/7, 3/7), axis = 't', align = "v")

tiff("prevplot_tl_complete.tif", width = 12, height = 5, units = "in", res = 400)
plot4
dev.off()

png("prevplot_tl_complete.png", width = 12, height = 5, units = "in", res = 400)
plot4
dev.off()

#discoloration
sum_dc <- cond %>% subset(cond_1 == "DC" &  perc_1 != "A"| cond_2 == "DC" & perc_2 != "A") %>%
  group_by(time_point, location_name, scientific_name) %>%
  summarize(n_dc = n())

sum_dc$time_point_a <- sum_dc$time_point
sum_dc$location_name_a <- sum_dc$location_name
sum_dc$scientific_name_a <- sum_dc$scientific_name

sum_dc <- sum_dc %>% unite("Event", c("time_point_a", "location_name_a", "scientific_name_a"), sep = "_")


prev_dc <- sum_demo %>% left_join(sum_dc, by = "Event") %>% 
  dplyr::select(-c(time_point.y:scientific_name.y)) %>%
  replace(is.na(.), 0) %>%
  mutate(prev_dc = n_dc/total_nj)

prev_dc <- na.omit(prev_dc)

colnames(prev_dc) <- gsub("\\.x","",colnames(prev_dc))

target_prev_dc <- prev_dc %>% subset(
  scientific_name == "SSID")

TimeLevels = c('Oct19','May22','Dec22')

target_prev_dc <- target_prev_dc %>% mutate_at(.vars = vars("time_point"), 
                                               .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


prevp_dc <- ggplot() +
  geom_violin(data = target_prev_dc, aes(x = time_point, y = prev_dc), fill = "gray80") +
  geom_jitter(data = target_prev_dc, aes(x = time_point, y = prev_dc, fill = location_name), size = 2.5, pch = 21,  
              width = 0.1, 
              height = 0) +
  facet_grid(~scientific_name, scales = 'free') +
  scale_y_continuous("Discoloration Prevalence") +
  scale_x_discrete("") +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10))

tiff("prevplot_ssid_dc.tif", width = 8, height = 4, units = "in", res = 400)
prevp_dc
dev.off()

#################
#Fate Tracking
#################
fate <- read_csv("CBC_ColonyData.csv")

mmea <- fate %>% subset(Species == "MMEA")

fate <- fate %>% separate_wider_delim(Date_InitialTag, "/",
                  names = c("Month_Initial", "Day_Initial", "Year_Initial"))

fate2 <- fate %>% subset(Year_Initial == "19")

fate2 <- fate2 %>% mutate("Date_Disease" = Date_DocumentedDisease, 
                          "Date_Dead" = Date_DocumentedMortality)

fate2 <- fate2 %>% mutate(Date_Disease = ifelse(is.na(Date_Disease),
                                 "Unknown", Date_Disease))

words <- c("Healthy","Unknown")

for (i in 1:length(fate2$Date_Disease)) {
  string <- fate2$Date_Disease[i]
  findWords <- intersect(unlist(strsplit(string," ")),words)
  if (!is.null(findWords)) {
    for (j in findWords) {
      fate2$Date_Disease[i] <- gsub(j,paste0(j,"/ / "),string)
    }
  }
}

split <- strsplit(fate2$Date_Disease, "/")
split2 <- matrix(unlist(split), ncol=3, byrow=TRUE)
split2 <- as.data.frame(split2)

split2 <- split2 %>% dplyr::select(-c(V2)) %>%
  rename("Month_Disease" = V1, "Year_Disease" = V3)

fate3 <- cbind(fate2, split2)


fate3 <- fate3 %>% mutate(Date_Dead = ifelse(is.na(Date_Dead),
                                                "Unknown", Date_Dead))

words <- c("Diseased","Unknown", "Healthy")

for (i in 1:length(fate3$Date_Dead)) {
  string <- fate3$Date_Dead[i]
  findWords <- intersect(unlist(strsplit(string," ")),words)
  if (!is.null(findWords)) {
    for (j in findWords) {
      fate3$Date_Dead[i] <- gsub(j,paste0(j,"/ / "),string)
    }
  }
}

split3 <- strsplit(fate3$Date_Dead, "/")
split4 <- matrix(unlist(split3), ncol=3, byrow=TRUE)
split4 <- as.data.frame(split4)

split4 <- split4 %>% dplyr::select(-c(V2)) %>%
  rename("Month_Dead" = V1, "Year_Dead" = V3)

fate4 <- cbind(fate3, split4)

fate5 <- fate4 %>% mutate(Cond_May22 = Month_Disease) %>%
                         mutate(Cond_May22 = ifelse(Month_Disease == "5" &
                                                  Year_Disease == "22",
                                                   "Diseased", Cond_May22)) %>%
                        mutate(Cond_May22 = ifelse(Month_Dead == 5 &
                                                  Year_Dead == 22, 
                                                  "Dead", Cond_May22)) %>%
                        mutate(Cond_May22 = ifelse(Month_Disease == 12 &
                                                   Year_Disease == 22, 
                                                  "Healthy", Cond_May22))

mmea <- fate5 %>% subset(Species == "MMEA")

fate6 <- fate5 %>% mutate(Cond_Dec22 = Month_Dead) %>%
  mutate(Cond_Dec22 = ifelse(Month_Disease == "12" &
                               Year_Disease == "22",
                             "Diseased", Cond_Dec22)) %>%
  mutate(Cond_Dec22 = ifelse(Month_Dead == "12" &
                               Year_Dead == "22",
                             "Dead", Cond_Dec22)) %>%
  mutate(Cond_Dec22 = ifelse(Cond_May22 == "Dead", 
                             "Dead", Cond_Dec22)) %>%
  mutate(Cond_Dec22 = ifelse(Month_Disease == "Healthy" &
                               Year_Dead == "23", 
                             "Healthy", Cond_Dec22)) %>%
  mutate(Cond_Dec22 = ifelse(Month_Disease == "Healthy" &
                               Year_Dead == "24", 
                             "Healthy", Cond_Dec22)) %>%
mutate(Cond_Dec22 = ifelse(Year_Disease == "22" &
                             Year_Dead == "23", 
                           "Diseased", Cond_Dec22)) %>%
  mutate(Cond_Dec22 = ifelse(Year_Disease == "22" &
                               Year_Dead == "24", 
                             "Diseased", Cond_Dec22)) %>%
  mutate(Cond_19 = "Healthy")

  
Long <-melt(fate6, id.vars =c("Species","OldTagNuim","NewTagNum",
                                  "MaxDiameter","Height","Transect"), measure.vars = 
              c("Cond_May22", "Cond_Dec22", "Cond_19")) %>%
  rename("Time_Point" = variable, "Condition" = value) 

Long$Time_Point <- recode(Long$Time_Point, "Cond_May22" = "May22",
                          "Cond_Dec22" = "Dec22",
                          "Cond_19" = "Oct19")


SpecLevels = c('SSID','MCAV','PSTR','PAST','MMEA','CNAT')

Long <- Long %>% mutate_at(.vars = vars("Time_Point"), 
                           .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
  mutate_at(.vars = vars("Species"),
            .funs = funs(factor(.,levels = SpecLevels, ordered = TRUE)))

check <- Long %>% subset(Species == "SSID" & Time_Point == "Oct19"|Species == "SSID" & Time_Point =="May22") %>%
  subset(Transect == "CBC30N") %>% arrange(OldTagNuim)

check2 <- Long %>% subset(Species == "SSID" & Time_Point == "Dec22" & Transect == "CBC30N")

sum <- check %>% group_by(Transect, Time_Point) %>% summarize(n = n())

sumsimple <- Long %>% group_by(Time_Point,Species) %>% count(Condition) %>% ungroup()

condcolors2 = c('Dead'='coral3','Diseased'='gold1','Healthy'='springgreen4',
                'Increased Old Mortality'='gray60',
                'Not Found' = 'black')


bysp <- ggplot(sumsimple, aes(x = Time_Point, y = n, fill = Condition)) +
  facet_wrap(~Species) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0,30)) +
  scale_x_discrete("Species") +
  scale_fill_manual("Condition", values = c(condcolors2)) +
  labs(y="Count") +
  labs(fill = "Condition") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank())

tiff("ConditionbySp.tif",width = 6, height = 6, units = "in", res = 300)
bysp
dev.off()

percents_all <- fate %>% group_by(Time_Point) %>% mutate(Total = n()) %>%
  ungroup() %>%
  group_by(Time_Point, Total, Condition) %>% summarize(Count = n()) %>%
  mutate(Prop = Count/Total) %>% arrange(Time_Point, Condition)

percents_sus <- fate %>% subset(Species != "PAST") %>%
  group_by(Time_Point) %>% mutate(Total = n()) %>%
  ungroup() %>%
  group_by(Time_Point, Total, Condition) %>% summarize(Count = n()) %>%
  mutate(Prop = Count/Total) %>% arrange(Time_Point, Condition)

percents_sp <- fate %>% group_by(Time_Point, Species) %>% mutate(Total = n()) %>%
  ungroup() %>%
  group_by(Time_Point, Species, Total, Condition) %>% summarize(Count = n()) %>%
  mutate(Prop = Count/Total) %>% arrange(Species, Time_Point, Condition)







########################
#Retired Code from old .csvs
########################

cbcMay <- read_csv("CBC Fate Tracking.csv")
cbcDec <- read_csv("Dec22fate.csv")

cbcMay$Transect <- recode(cbcMay$Transect,"Lagoon" = "CBC Lagoon")
cbcDec$TransectName <- recode(cbcDec$TransectName, "CBCLAGOON" = "CBC Lagoon",
                              "CURLEW" = "Curlew")

cbcMay$Date_Initial <- as.character(cbcMay$Date_Initial)
split <- strsplit(cbcMay$Date_Initial, "-")
split2 <- matrix(unlist(split), ncol=3, byrow=TRUE)
split2 <- as.data.frame(split2)
cbcMay2 <- cbind(cbcMay, split2) %>%
  subset(V1 != "2022") %>% droplevels()

cbc19 <- cbcMay2 %>% mutate(Time_Point = "Oct19", Condition = "Healthy") %>% 
  dplyr::select(Time_Point,Date_Initial,Transect,Species,Tag_Orig,New_Tag,Max_Diameter,Height,Condition)

cbc19$Date_Initial <- as.character(cbc19$Date_Initial)

cbcMay <- cbcMay %>% mutate(Time_Point = "May22") %>% 
  rename("Condition1May" = "Condition_May_2022") %>%
  dplyr::select(Time_Point,Date_Initial,Transect,Species,Tag_Orig,New_Tag,Max_Diameter,Height,Condition1May) %>%
  mutate(New_Tag = if_else(is.na(New_Tag),
                                Tag_Orig, New_Tag))
  #correct tag number for MCAV24 (double check this)
cbcMay <- cbcMay %>%
  mutate(New_Tag = case_when(
    Tag_Orig == "380" & Transect == "CBC30N" & Species == "MCAV" ~ "24",
    TRUE ~ as.character(New_Tag)))

cbcMay$New_Tag_a <- cbcMay$New_Tag
cbcMay$Transect_a <- cbcMay$Transect
cbcMay <- cbcMay %>% unite("Coral_ID", c("Transect_a", "New_Tag_a"), sep = "_") %>%
  rename("Date_Initial1" = "Date_Initial") 


cbcDec <- cbcDec %>% mutate(Time_Point = "Dec22") %>%
  rename("Date_Initial" = "Date_InitialTag", "Transect" = "TransectName",
         "New_Tag" = "NewTag", "Max_Diameter" = "MaxDiameter",
         "Condition2Dec" = "Condition") %>%
  dplyr::select(Time_Point,Date_Initial,Transect,Species,New_Tag,Max_Diameter,Height,Condition2Dec) %>%
  distinct() 


#make corrections (based on photos and notes from fate script V2 - double check)

cbcDec$New_Tag_a <- cbcDec$New_Tag
cbcDec$Transect_a <- cbcDec$Transect

cbcDec <- cbcDec %>% unite("Coral_ID", c("Transect_a", "New_Tag_a"), sep = "_") %>%
  arrange(Coral_ID) %>%
  subset(!(is.na(Condition2Dec)))

#unite May and Dec
cbcMayDec <- cbcDec %>%
    right_join(cbcMay, by = "Coral_ID")

cbcMayDec <- cbcMayDec[,order(colnames(cbcMayDec))]
  
cbcMayDec$Date_Initial1 <- as.character(cbcMayDec$Date_Initial1)
split <- strsplit(cbcMayDec$Date_Initial1, "-")
split2 <- matrix(unlist(split), ncol=3, byrow=TRUE)
split2 <- as.data.frame(split2)
cbcMayDec <- cbind(cbcMayDec, split2) %>%
  subset(V1 != "2022") %>% droplevels()

#correct errors (December conditions that were missing or wrong)
cbcMayDec <- cbcMayDec %>%
  mutate(Condition2Dec = case_when(
    Coral_ID == "SR30N_51" ~ "Not Found",
    TRUE ~ as.character(Condition2Dec))) %>%
  mutate(Condition2Dec = case_when(
    Coral_ID == "CBC30N_11" ~ "Diseased",
    TRUE ~ as.character(Condition2Dec))) %>%
  mutate(Condition2Dec = case_when(
    Coral_ID == "CBC Lagoon_13" ~ "Increased Old Mortality",
    TRUE ~ as.character(Condition2Dec))) %>%
  mutate(Condition2Dec = case_when(
    Coral_ID == "CBC30N_19" ~ "Dead",
    TRUE ~ as.character(Condition2Dec))) %>%
  mutate(Condition2Dec = case_when(
    Coral_ID == "SR30N_59" ~ "Increased Old Mortality",
    TRUE ~ as.character(Condition2Dec))) %>%
  mutate(Condition2Dec = case_when(
    Coral_ID == "CBC30N_5" ~ "Diseased",
    TRUE ~ as.character(Condition2Dec)))

cbcMayDec$Condition1May <- recode(cbcMayDec$Condition1May, "Damage, mostly dead" = "Dead",
                        "Dead - recent dead" = "Dead",
                        "Diseased - mostly dead" = "Diseased",
                        "Diseased/Arrested TL" = "Diseased", #need to check this one but I had it as diseased
                        "Gone - storm?" = "Physical damage - no remains",
                        "Healthy - increased partial mortality" = "Increased Old Mortality",
                        "Healthy?" = "Healthy",
                        "Healthyish" = "Diseased",
                        "minor paling" = "Color Loss",
                        "Mostly dead - old mortality" = "Increased Old Mortality",
                        "Multifocal Color Loss" = "Color Loss",
                        "Physical damage?" = "Physical damage",
                        "Not found" = "Not Found",
                        "Increased Old Mort" = "Increased Old Mortality")

#simplify conditions
cbcMayDec$Condition1May <- recode(cbcMayDec$Condition1May, "Physical damage - no remains" = "Not Found",
                         "Increased Old Mort" = "Increased Old Mortality",
                         "Physical damage" = "Healthy",                         
                         "Color Loss" = "Healthy",
                         "Presumed dead" = "Not Found")

cbcMayDec$Condition2Dec <- recode(cbcMayDec$Condition2Dec, "DISEASED" = "Diseased",
                                  "DEAD" = "Dead",
                                  "HEALTHY" = "Healthy")

#Replace NA December with May condition 

nas <- cbcMayDec %>% subset(is.na(Condition2Dec)) #check that all remaining are dead or not found

cbcMayDec <- cbcMayDec %>%
  mutate(Condition2Dec = case_when(
    is.na(Condition2Dec) ~ Condition1May,
    TRUE ~ as.character(Condition2Dec)))

Long <-melt(cbcMayDec, id.vars =c("Species.y","Date_Initial1","Max_Diameter.y",
                                  "Height.y","New_Tag.y","Tag_Orig","Transect.y"), measure.vars = 
              c("Condition1May", "Condition2Dec")) %>%
  rename("Time_Point" = variable, "Condition" = value,
         "Date_Initial" = Date_Initial1) %>%
  mutate(Time_Point = if_else(Time_Point == "Condition1May",
                                "May22", Time_Point)) %>%
  mutate(Time_Point = if_else(Time_Point == "Condition2Dec",
                              "Dec22", Time_Point))
  
colnames(Long) <- gsub("\\.y","",colnames(Long))

fate <- rbind(Long, cbc19)

fate$Species <- recode(fate$Species, "DSTR" = "PSTR")
SpecLevels = c('SSID','MCAV','PSTR','PAST','MMEA','CNAT')

fate <- fate %>% mutate_at(.vars = vars("Time_Point"), 
                             .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
  mutate_at(.vars = vars("Species"),
            .funs = funs(factor(.,levels = SpecLevels, ordered = TRUE)))

check <- fate %>% subset(Species == "SSID" & Time_Point == "Oct19"|Species == "SSID" & Time_Point =="May22") %>%
  subset(Transect == "CBC30N") %>% arrange(Tag_Orig)

check2 <- fate %>% subset(Species == "SSID" & Time_Point == "Dec22" & Transect == "CBC30N")

sum <- check %>% group_by(Transect, Time_Point) %>% summarize(n = n())

sumsimple <- fate %>% group_by(Time_Point,Species) %>% count(Condition) %>% ungroup()

condcolors2 = c('Dead'='coral3','Diseased'='gold1','Healthy'='springgreen4',
                'Increased Old Mortality'='gray60',
                'Not Found' = 'black')


bysp <- ggplot(sumsimple, aes(x = Time_Point, y = n, fill = Condition)) +
  facet_wrap(~Species) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0,30)) +
  scale_x_discrete("Species") +
  scale_fill_manual("Condition", values = c(condcolors2)) +
  labs(y="Count") +
  labs(fill = "Condition") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank())

tiff("ConditionbySp.tif",width = 6, height = 6, units = "in", res = 300)
bysp
dev.off()

percents_all <- fate %>% group_by(Time_Point) %>% mutate(Total = n()) %>%
  ungroup() %>%
  group_by(Time_Point, Total, Condition) %>% summarize(Count = n()) %>%
  mutate(Prop = Count/Total) %>% arrange(Time_Point, Condition)

percents_sus <- fate %>% subset(Species != "PAST") %>%
  group_by(Time_Point) %>% mutate(Total = n()) %>%
  ungroup() %>%
  group_by(Time_Point, Total, Condition) %>% summarize(Count = n()) %>%
  mutate(Prop = Count/Total) %>% arrange(Time_Point, Condition)

percents_sp <- fate %>% group_by(Time_Point, Species) %>% mutate(Total = n()) %>%
  ungroup() %>%
  group_by(Time_Point, Species, Total, Condition) %>% summarize(Count = n()) %>%
  mutate(Prop = Count/Total) %>% arrange(Species, Time_Point, Condition)
  
