library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 

#To do

#GET RID OF LOCATION FOR DEMO DATA
#arcsin transformation for prevalence data
#try GLMTBB for zero inflated density data
#try GAM for disease prevalence
#model peak disease prevalence against starting colony density


# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/GitHub/SCTLD_demographics/")


# Import the survey data
cover <- read.csv("annotations_CBC_19-22.csv")
#one <- read.csv("CBCdemo19.csv")
#two <- read.csv("CBCdemoMay22.csv")
#three <- read.csv("CBCdemoDec22.csv")

split <- strsplit(cover$Date, "/")
split2 <- matrix(unlist(split), ncol=3, byrow=TRUE)
split2 <- as.data.frame(split2)
cover <- cbind(cover, split2)

cover <- cover %>% rename("Year" = "V3", "Day" = "V2", "Month" = "V1") %>%
  select(Name, Date, Habitat, SiteName, Year, Month, Day, Label, Row, Column)


levels(as.factor(cover$SiteName))

cover$SiteName <- recode(cover$SiteName, "CBC 30 N" = "CBC30N",
                         "CBC House Reef" = "House Reef")

cover$Year <-  paste("20", cover$Year, sep="")

for (i in 1:nrow(cover)) {
  
  if (nchar(cover$Month[i]) == 1) {
    
    cover$Month[i] <- paste0(0, cover$Month[i])
  } 
}

cover$Year_a <- cover$Year
cover$Month_a <- cover$Month
cover$Day_a <- cover$Day
cover$SiteName_a <- cover$SiteName


cover <- cover %>% subset(SiteName != "Tobacco Reef") %>% droplevels()

cover <- cover %>% unite("Date", c("Year_a","Month_a","Day_a"), sep = "-")

cover$Date <- as.Date(cover$Date)

cover$Year_a <- cover$Year
cover$Month_a <- cover$Month

cover <- cover %>% unite("TimePoint", c("Year_a","Month_a"), sep = "-")

levels(as.factor(cover$TimePoint))



check <- cover %>% subset(TimePoint == "2019-11")

#cover$TimePoint <- recode(cover$TimePoint, "20-01" = "19-10")

cover$TimePoint <- recode(cover$TimePoint, "2019-10" = "October19",
                          "2022-12" = "December22", 
                          "2022-05" = "May22",
                          "2020-01" = "January20")

#cover <- cover %>% mutate("TimeDate" = TimePoint)

#cover$TimeDate <- recode(cover$TimeDate,
#  "Oct19" = "2019-11-20",
# "May22" = "2022-05-23",
# "Dec22" = "2022-12-03")

#cover$TimeDate <- as.Date(cover$TimeDate)

cover$TimePoint_a <- cover$TimePoint

cover <- cover %>% unite("Event", c("TimePoint_a", "SiteName_a"), sep = "_")


cover1 <- cover %>% group_by(Event) %>% summarize(npics = (length(unique(Name)))) %>%
  mutate(npoints = npics * 40)

coverx <- left_join(cover, cover1, by = "Event")

cover2 <- coverx %>% 
  group_by(Habitat, SiteName, Date, TimePoint, Label, npoints) %>%
  summarize(n = length(Label)) %>%
  mutate(cover = (n/npoints)*100) %>%
  ungroup() %>%
  complete(Label, nesting(Habitat, SiteName, Date, npoints, TimePoint), 
           fill = list(n = 0, cover = 0)) %>%
  arrange(Habitat, SiteName, npoints, TimePoint, Label) %>%
  select(Habitat, SiteName, Date, npoints, TimePoint, Label, n, cover)

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
##################
#fill out time points
##################

TimeLevels <- c("October19","November19","December19",
                "January20", "February20", "March20" ,"April20","May20", "June20","July20",
                "August20" ,"September20", "October20","November20","December20",
                "January21","February21","March21" ,"April21","May21", "June21",
                "July21","August21" ,"September21", "October21","November21","December21",
                "January22","February22","March22" ,"April22","May22", "June22","July22",
                "August22" ,"September22", "October22","November22","December22")


LabLevels = c('Turf Algae','All Stony Coral','Macroalgae','Octocoral')


cover <- cover %>% 
  mutate_at(.vars = vars("Label_General"),
            .funs = funs(factor(.,levels = LabLevels, ordered = TRUE))) %>% 
  mutate_at(.vars = vars("TimePoint"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

cover <- cover %>% mutate(Survey = TimePoint)
cover$Survey <- recode(cover$Survey, "October19" = "November19",
                       "January20" = "November19")

covgen <- cover %>% group_by(SiteName, Habitat, Date, Survey, npoints, Label_General) %>% 
  summarize(sum_gen = sum(n)) %>% 
  mutate(cov_gen = (sum_gen/npoints)*100) %>% 
  mutate(cov_gensq = sqrt(cov_gen)) %>%
  mutate(cov_prop = cov_gen/100) %>%
  subset(!(is.na(Label_General))) 

covgen$Survey_a <- covgen$Survey
covgen$Label_General_a <- covgen$Label_General

covgen <- covgen %>% unite("GroupEvent", c(Label_General_a, Survey_a))

covgroup <- covgen %>% group_by(Survey, Label_General) %>%
  summarize(MeanCov = mean(cov_gen), seCov = se(cov_gen), maxCov = max(cov_gen)) %>%
  unite("GroupEvent", c(Label_General,Survey))

covgen <- left_join(covgen, covgroup, by = "GroupEvent")


suscep <- cover %>% group_by(SiteName, Habitat, Date, Survey, npoints, Label_Suscep) %>% 
  summarize(sum_gen = sum(n)) %>%
  mutate(cov_gen = (sum_gen/npoints)*100) %>% subset(Label_Suscep == "SCTLD Suscep. Spp.")

covsto <- covgen %>% subset(Label_General == "All Stony Coral") %>%
  mutate(cov_gensq = sqrt(cov_gen)) %>%
  mutate(cov_prop = cov_gen/100)


suscep <- suscep %>% 
  mutate(cov_gensq = sqrt(cov_gen)) %>%
  mutate(cov_prop = cov_gen/100) %>%
  unite("GroupEvent", c(Label_Suscep,Survey), remove = FALSE)


susgroup <- suscep %>% group_by(Survey, Label_Suscep) %>%
  summarize(MeanCov = mean(cov_gen), seCov = se(cov_gen), maxCov = max(cov_gen)) %>%
  unite("GroupEvent", c(Label_Suscep,Survey))

suscep <- suscep %>% left_join(susgroup, by = "GroupEvent")

library(betareg)
library(lme4)
library(emmeans)
library(rcompanion)
library(glmmTMB)



#mod5 <- glmmTMB(cov_prop ~ TimePoint +
# (1|SiteName),
# data=covsto,beta_family())
#hist(resid(mod4))
#Anova(mod4)

#########
#Full Model
#########
mod1 <- glmer(cov_prop ~ Survey*Label_General +
                (1|SiteName),
              weights=npoints,
              data=covgen,family="binomial")
hist(resid(mod1))
Anova(mod1)

emm <- emmeans(mod1, ~ Survey*Label_General)
simple <- pairs(emm, simple = "Survey")
pairwise <- as.data.frame(pairs(emm, simple = "Survey"))


isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}



df <- data.frame()

Labels <- levels(as.factor(pairwise$Label_General))

for(current_Label in Labels) {
  
  Lab_df <- pairwise %>% subset(Label_General == current_Label) 
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = Lab_df,
                                   threshold = 0.05)) %>%
    mutate("Label_General" = current_Label)
  
  df <- df %>%
    bind_rows(letters)
}




#devtools::install_github("ianmoran11/mmtable2")
library(mmtable2)
library(gt)
library(kableExtra)




MainLetters <- df %>% unite("GroupEvent", c(Label_General,Group))


covgenx <- covgen %>% 
  left_join(MainLetters, by = "GroupEvent") %>%
  mutate(TimePoint = case_when(
    SiteName == "SR30N" & Survey == "November19" ~ "January20",
    SiteName == "CBC30N" & Survey == "November19" ~ "January20",
    SiteName == "CBC Central" & Survey == "November19" ~ "October19",
    SiteName == "South Reef Central" & Survey == "November19" ~ "October19",
    SiteName == "House Reef" & Survey == "November19" ~ "October19",
    SiteName == "Curlew Patch" & Survey == "November19" ~ "October19",
    SiteName == "CBC Lagoon" & Survey == "November19" ~ "October19",
    TRUE ~ as.factor(Survey))) %>%
  mutate_at(.vars = vars("TimePoint"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))



covgen1 <- covgenx %>% subset(Label_General == "Turf Algae"|
                                Label_General == "All Stony Coral")

covgen2 <- covgenx %>% subset(Label_General == "Macroalgae"|
                                Label_General == "Octocoral")

modsus <- glm(cov_prop ~ Survey + SiteName, data = suscep, family = "quasibinomial")
modsus <- glmer(cov_prop ~ Survey +
                  (1|SiteName),
                weights=npoints,
                data=suscep,family="binomial")
hist(resid(modsus))
Anova(modsus)
rsquared(modsus)

susresult <- as.data.frame(Anova(modsus)) %>%
  mutate(Label_General = "SCTLD Susceptible Coral")
susresult$sig <- isSig(susresult$"Pr(>Chisq)")
mod_df <- data.frame()
mod_df <- mod_df %>% rbind(susresult)

covtable <- mod_df %>%
  remove_rownames %>% column_to_rownames(var="Label_General") %>%
  select(-c(sig)) %>%
  rename('p-value' = 'Pr(>Chisq)') %>%
  mutate('p-value' = '<0.001') %>%
  kbl() %>%
  kable_styling()

png("CovTable.png", width = 8, height = 6, units = "in", res = 300)
covtable
dev.off()

suspairwise <- as.data.frame(pairs(emmeans(modsus, ~ Survey)))

suspairwise$sig <- isSig(suspairwise$p.value)

susletters <- as.data.frame(cldList(p.value ~ contrast,
                                    data = suspairwise,
                                    threshold = 0.05))

suscep <- suscep %>% left_join(susletters, by = c("Survey" = "Group")) %>%
  mutate(TimePoint = case_when(
    SiteName == "SR30N" & Survey == "November19" ~ "January20",
    SiteName == "CBC30N" & Survey == "November19" ~ "January20",
    SiteName == "CBC Central" & Survey == "November19" ~ "October19",
    SiteName == "South Reef Central" & Survey == "November19" ~ "October19",
    SiteName == "House Reef" & Survey == "November19" ~ "October19",
    SiteName == "Curlew Patch" & Survey == "November19" ~ "October19",
    SiteName == "CBC Lagoon" & Survey == "November19" ~ "October19",
    TRUE ~ as.factor(Survey))) %>%
  mutate_at(.vars = vars("TimePoint"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

#suscep <- suscep %>% mutate(Survey = fct_relevel(Survey,
# "Novemer19", "May22", "Deccember22")) 


library(viridis)
library(cowplot)
library(gridExtra)

levels(as.factor(covgen$SiteName))

sitecolors = c('CBC Central'='red3','CBC Lagoon'='darkorchid','CBC30N'='gold1',
               'Curlew Patch' = 'blue', 'House Reef' ='seagreen4', 
               'South Reef Central' = 'orange', 'SR30N' = 'pink')

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}


labyearp1 <- ggplot() +
  geom_jitter(data = covgen1, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 2, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  geom_text(data = covgen1, 
            aes(x = Survey, y = maxCov, label = Letter), nudge_y = 3) +
  geom_violin(data = covgen1, aes(x = Survey, y = cov_gen), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = covgen1, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 4, alpha = 0.5, pch = 21,  
              width = 0.25, 
              height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~Label_General, scales = "free") +
  scale_y_continuous("Percent Cover") +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(2, 2))) +   
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
  geom_jitter(data = covgen2, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 4, alpha = 0, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  geom_text(data = covgen2, 
            aes(x = Survey, y = maxCov, label = Letter), nudge_y = 3) +
  geom_violin(data = covgen2, aes(x = Survey, y = cov_gen), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = covgen2, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 4, alpha = 0.5, pch = 21,  
              width = 0.25, 
              height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~Label_General, scales = "free") +
  scale_y_continuous("Percent Cover") +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20','July21',
                                                'May22','December22'),
                   expand = expansion(add = c(2, 2))) +   
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
  geom_jitter(data = suscep, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 4, alpha = 0, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  geom_text(data = suscep, 
            aes(x = Survey, y = maxCov, label = Letter), nudge_y = 1) +
  geom_violin(data = suscep, aes(x = Survey, y = cov_gen), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = suscep, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 4, alpha = 0.5, pch = 21,  
              width = 0.25, 
              height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~Label_Suscep, scales = "free") +
  scale_y_continuous("") +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20','July21',
                                                'May22','December22'),
                   expand = expansion(add = c(2, 2))) +   
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

plot1 <- grid.arrange(labyearp1, labyearp2, ncol=1, nrow =2)

plot2 <- cowplot::plot_grid(plot1, legend, rel_widths = c(4/5, 1/5), axis = 't', align = "v")


tiff("Figures/CBCCoverViolin.tif",width = 7, height = 7, units = "in", res = 300)
plot2
dev.off()

png("Figures/CBCCoverViolin.png", width = 7, height = 7, units = "in", res = 300)
plot2
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

category <- sto %>% group_by(SiteName, Habitat, TimePoint, npoints, Survey, Category) %>% 
  summarize(sum_cat = sum(n)) %>%
  mutate(cov_cat = (sum_cat/npoints)*100) %>%
  mutate(cov_prop = cov_cat/100)



Cats <- levels(as.factor(category$Category))


#Single Model
mod2 <- glmer(cov_prop ~ Survey*Category +
                (1|SiteName),
              weights=npoints,
              data=category,family="binomial")

hist(resid(mod2))
Anova(mod2)


emm <- emmeans(mod2, ~ Survey*Category)
simple <- pairs(emm, simple = "Survey")
pairwise <- as.data.frame(pairs(emm, simple = "Survey"))

df <- data.frame()

Cats <- levels(as.factor(pairwise$Category))

for(current_Cat in Cats) {
  
  Cat_df <- pairwise %>% subset(Category == current_Cat) 
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = Cat_df,
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

orbi <- sto2 %>% subset(Species == "ORBI"|Species == "OFAV"|
                          Species == "OANN"|Species == "OANN") %>%
  group_by(Habitat, SiteName, Date, npoints, TimePoint, 
           Label_General, Label_Suscep, Survey, Category) %>%
  summarize(n = sum(n)) %>%
  mutate("Species" = "ORBI", cov_prop = n/npoints, cover = cov_prop*100)

sto2 <- sto2 %>% rbind(orbi) %>% subset(Species != "OFAV" & Species != "OANN") %>% droplevels()

SpecList <- levels(as.factor(sto2$Species))

#models
mod3 <- glmer(cov_prop ~ Survey*Species +
                (1|SiteName),
              weights=npoints,
              data=sto2,family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
hist(resid(mod3))
Anova(mod3)

df <- data.frame()


emm <- emmeans(mod3, ~ Survey*Species)
simple <- pairs(emm, simple = "Survey")
pairwise <- as.data.frame(pairs(emm, simple = "Survey"))

df <- data.frame()

Spec <- levels(as.factor(pairwise$Species))

for(current_Spec in Spec) {
  
  Spec_df <- pairwise %>% subset(Species == current_Spec) 
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = Spec_df,
                                   threshold = 0.05)) %>%
    mutate("Species" = current_Spec)
  
  df <- df %>%
    bind_rows(letters)
}




speclet <- df %>% unite("Event", c("Species", "Group"))

stomeans <- sto2 %>% group_by(Species, Survey) %>%
  summarize(MeanCov = mean(cover), seCov = se(cover), maxCov = max(cover)) %>%
  unite("Event", c(Species,Survey))

sto3 <- sto2 %>% 
  unite("Event", c("Species", "Survey"), remove = FALSE) %>%
  left_join(speclet, by = "Event") %>%
  left_join(stomeans, by = "Event") %>%
  mutate(TimePoint = case_when(
    SiteName == "SR30N" & Survey == "November19" ~ "January20",
    SiteName == "CBC30N" & Survey == "November19" ~ "January20",
    SiteName == "CBC Central" & Survey == "November19" ~ "October19",
    SiteName == "South Reef Central" & Survey == "November19" ~ "October19",
    SiteName == "House Reef" & Survey == "November19" ~ "October19",
    SiteName == "Curlew Patch" & Survey == "November19" ~ "October19",
    SiteName == "CBC Lagoon" & Survey == "November19" ~ "October19",
    TRUE ~ as.factor(Survey))) %>%
  mutate_at(.vars = vars("TimePoint"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
  mutate(Letter = ifelse(Species == "AAGA" |
                           Species == "PPOR", 
                         "", Letter))

#sto3 <- sto2 %>% unite("Event", c("Species", "TimePoint"), remove = FALSE) %>%
#left_join(speclet, by = "Event") %>%
#left_join(stomeans, by = "Event")


library(ggh4x)

targetp1 <- ggplot() +
  geom_jitter(data = sto3, aes(x = TimePoint, y = cover, fill = SiteName), size = 4, pch = 21, alpha=0,  
              width = 0.25, 
              height = 0) +
  geom_text(data = sto3, 
            aes(x = Survey, y = maxCov, label = Letter), nudge_y = 1.5) +
  geom_violin(data = sto3, aes(x = Survey, y = cover), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = sto3, aes(x = TimePoint, y = cover, fill = SiteName), size = 4, pch = 21, alpha=0.5,  
              width = 0.25, 
              height = 0)+
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_nested_wrap(~Category + Species, scales = "free",
                    ncol = 4, labeller = labeller(Category = label_wrap_gen(width = 16))) +
  scale_y_continuous("Percent Cover") +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20','July21',
                                                'May22','December22'),
                   expand = expansion(add = c(4, 4))) +   
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10))

tiff("Figures/CBCTargetsBoxplot.tif",width = 10, height = 8, units = "in", res = 300)
targetp1
dev.off()

png("Figures/CBCTargetsBoxplot.png", width = 10, height = 8, units = "in", res = 300)
targetp1
dev.off()

############################
#DEMO
############################
demo <- read.csv("CBC_Coral_Demographics.csv")
cond <- read.csv("CBC_Coral_Conditions.csv")

demo <- demo %>% subset(year < 2023)

dead <- demo %>% subset(X0_recent_total_mortality > 0) #all from May22 - LH put recent mort in size bins, need to remove

dead2 <- cond %>% subset(percent_mortality == 100)
#3 EFAS and 1 PSTR at SR Central were totally dead!
#PSTR = 6cm
#EFAS at 6, 10, and 13cm

#this code was to remove the dead colonies from demo but
#they are not in the GitHub file
#demo <- demo %>%
# mutate(X5_to_10_cm = ifelse(data_collector == "Leah Harper" &
#                              sample_collection_month == 5 &
#                             sample_collection_year == 2022 & 
#                            location_name == "South Reef Central" &
#                           scientific_name == "EFAS", 
#                        X5_to_10_cm - 2, X5_to_10_cm)) %>%
#mutate(X11_to_20_cm = ifelse(data_collector == "Leah Harper" &
#                              sample_collection_month == 5 &
#                             sample_collection_year == 2022 & 
#                            location_name == "South Reef Central" &
#                           scientific_name == "EFAS", 
#                        X11_to_20_cm - 1, X11_to_20_cm)) %>%
#mutate(X5_to_10_cm = ifelse(data_collector == "Leah Harper" &
#                             sample_collection_month == 5 &
#                            sample_collection_year == 2022 & 
#                           location_name == "South Reef Central" &
#                          scientific_name == "PSTR", 
#                       X5_to_10_cm - 1, X5_to_10_cm))

#demo <- demo %>% dplyr::select(-c(transect:coral_demographics_notes))

demo <- demo %>%
  replace(is.na(.), 0) %>%
  mutate(total_small = rowSums(.[17:19])) %>% #5-40cm
  mutate(total_lg = rowSums(.[20:21])) %>% #>40cm
  mutate(adult = rowSums(.[17:21])) #>5cm

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

demo$year <- as.factor(demo$year)
levels(demo$year)

demo$Year_a <- demo$year
demo$Month_a <- demo$month
demo$Day_a <- demo$day
demo$location_name_a <- demo$location_name


demo <- demo %>% subset(location_name != "Tobacco Reef") %>% droplevels()

demo <- demo %>% unite("date", c("Year_a","Month_a","Day_a"), sep = "-")

demo$date <- as.Date(demo$date)

demo$Year_a <- demo$year
demo$Month_a <- demo$month

demo <- demo %>% unite("time_point", c("Year_a","Month_a"), sep = "-")

levels(as.factor(demo$time_point))

demo$time_point <- recode(demo$time_point, "2019-10" = "October19",
                          "2022-12" = "December22", 
                          "2022-5" = "May22",
                          "2020-1" = "January20")

#cover <- cover %>% mutate("TimeDate" = TimePoint)


demo$scientific_name <- as.factor(demo$scientific_name)
levels(demo$scientific_name)

check <- demo %>% subset(scientific_name == "") 

demo$scientific_name <- as.character(demo$scientific_name)


demo <- demo %>% mutate(scientific_name = ifelse(code == "agar",
                                                 "Agaricia spp.", scientific_name))
demo$scientific_name <- as.factor(demo$scientific_name)

levels(demo$scientific_name)

########################
#NMDS
########################
library(vegan)

demo$time_point_a <- demo$time_point


demo <- demo %>%
  unite("event", c("location_name_a", "time_point_a"), sep = "_") %>%
  mutate("survey" = time_point)

demo$survey <- recode(demo$survey, "October19" = "January20")

cast <- demo %>% group_by(location_name, time_point, survey, event, scientific_name) %>%
  summarize(total = sum(adult)) %>%
  pivot_wider(id_cols = c(location_name, survey, event), names_from = scientific_name, 
              values_from = total)

#Jaccard (P/A)

mat = cast[,3:ncol(cast)]
mat <- mat %>% remove_rownames %>% column_to_rownames(var="event")
mat[is.na(mat)] <- 0
mat[mat > 0] <- 1
mat <- as.matrix(mat)
#mat <- sqrt(mat)
set.seed(123456)


NMDS1 <-
  metaMDS(mat,
          distance = "jaccard",
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
data.scores1$survey = cast$survey
data.scores1$event = cast$event

data.scores1 <- data.scores1 %>% 
  mutate(survey = fct_relevel(survey,
                              "January20", "May22", "December22")) 



species.scores1 <- as.data.frame(scores(NMDS1, "species")) 
species.scores1$species <- rownames(species.scores1)  # create a column of species, from the rownames of species.scores

uno <- data.scores1[data.scores1$survey == "January20", ][chull(data.scores1[data.scores1$survey == 
                                                                               "January20", c("NMDS1", "NMDS2")]), ]
dos <- data.scores1[data.scores1$survey == "May22", ][chull(data.scores1[data.scores1$survey == 
                                                                           "May22", c("NMDS1", "NMDS2")]), ]
tres <- data.scores1[data.scores1$survey == "December22", ][chull(data.scores1[data.scores1$survey == 
                                                                                 "December22", c("NMDS1", "NMDS2")]), ]

hull.data1 <- rbind(uno, dos, tres) %>% 
  mutate(survey = fct_relevel(survey,
                              "January20", "May22", "December22")) 


hull.data1

library(ggrepel)
library(viridis)

tiff("Figures/adultNMDS_jaccard.tif",width = 7, height = 7, units = "in", res = 400)
ggplot() + 
  #geom_polygon(data=hull.data1,aes(x=NMDS1,y=NMDS2,fill=survey,group=survey),alpha=0.30, color = "black") + 
  geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2,colour = survey, fill = survey, shape = location_name),size=4, alpha = 0.8) +
  geom_point(data=species.scores1,aes(x=NMDS1,y=NMDS2), pch = 20, color = "black", size=2) +
  geom_text_repel(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),size = 2, min.segment.length = 0.5, max.overlaps = 18) + 
  #geom_text_repel(data=data.scores1,aes(x=NMDS1,y=NMDS2,label=location_name),size = 4, colour = "black", min.segment.length = 0.5, max.overlaps = 18) + 
  scale_fill_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Site",values=c("CBC Central" = 21, "CBC Lagoon" = 22, "House Reef" = 23, "CBC30N" = 24,
                                     "South Reef Central" = 25, "Curlew Patch" = 3, "SR30N" = 9)) +
  #xlim(-1,2) +  
  ggtitle("Jaccard Dissimilarity") +
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


jacdf <- adonis2(formula = mat ~ location_name + survey, data = cast, permutations = 10000)
jactable <- jacdf %>%
  kbl() %>%
  kable_styling()

jactable

##Bray-Curtis

mat = cast[,3:ncol(cast)]
mat <- mat %>% remove_rownames %>% column_to_rownames(var="event")
mat[is.na(mat)] <- 0
#mat[mat > 0] <- 1
mat <- as.matrix(mat)
mat <- sqrt(mat)
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
data.scores1$survey = cast$survey
data.scores1$event = cast$event

data.scores1 <- data.scores1 %>% 
  mutate(survey = fct_relevel(survey,
                              "January20", "May22", "December22")) 



species.scores1 <- as.data.frame(scores(NMDS1, "species")) 
species.scores1$species <- rownames(species.scores1)  # create a column of species, from the rownames of species.scores

uno <- data.scores1[data.scores1$survey == "January20", ][chull(data.scores1[data.scores1$survey == 
                                                                               "January20", c("NMDS1", "NMDS2")]), ]
dos <- data.scores1[data.scores1$survey == "May22", ][chull(data.scores1[data.scores1$survey == 
                                                                           "May22", c("NMDS1", "NMDS2")]), ]
tres <- data.scores1[data.scores1$survey == "December22", ][chull(data.scores1[data.scores1$survey == 
                                                                                 "December22", c("NMDS1", "NMDS2")]), ]

hull.data1 <- rbind(uno, dos, tres) %>% 
  mutate(survey = fct_relevel(survey,
                              "January20", "May22", "December22")) 


hull.data1

library(ggrepel)
library(viridis)

tiff("Figures/adultNMDS_bray.tif",width = 7, height = 7, units = "in", res = 400)
ggplot() + 
  #geom_polygon(data=hull.data1,aes(x=NMDS1,y=NMDS2,fill=survey,group=survey),alpha=0.30, color = "black") + 
  geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2,colour = survey, fill = survey, shape = location_name),size=4, alpha = 0.8) +
  geom_point(data=species.scores1,aes(x=NMDS1,y=NMDS2), pch = 20, color = "black", size=2) +
  geom_text_repel(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),size = 2, min.segment.length = 0.5, max.overlaps = 18) + 
  #geom_text_repel(data=data.scores1,aes(x=NMDS1,y=NMDS2,label=location_name),size = 4, colour = "black", min.segment.length = 0.5, max.overlaps = 18) + 
  scale_fill_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Time Point", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Site",values=c("CBC Central" = 21, "CBC Lagoon" = 22, "House Reef" = 23, "CBC30N" = 24,
                                     "South Reef Central" = 25, "Curlew Patch" = 3, "SR30N" = 9)) +
  #xlim(-1,2) +  
  ggtitle("Bray-Curtis Dissimilarity") +
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


braydf <- adonis2(formula = mat ~ location_name + survey, data = cast, permutations = 10000)
braytable <- braydf %>%
  kbl() %>%
  kable_styling()

braytable



njlong <- melt(cast, id.vars =c(
  "location_name", "survey", "event")) %>%
  rename("scientific_name" = variable, "count_adult" = value) %>%
  mutate("time_point" = survey) %>%
  mutate(time_point = case_when(
    location_name == "SR30N" & survey == "January20" ~ "January20",
    location_name == "CBC30N" & survey == "January20" ~ "January20",
    location_name == "CBC Central" & survey == "January20" ~ "October19",
    location_name == "South Reef Central" & survey == "January20" ~ "October19",
    location_name == "House Reef" & survey == "January20" ~ "October19",
    location_name == "Curlew Patch" & survey == "January20" ~ "October19",
    location_name == "CBC Lagoon" & survey == "January20" ~ "October19",
    TRUE ~ as.factor(survey))) %>%
  mutate_at(.vars = vars("time_point"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


njlong <- njlong %>% 
  mutate(survey = case_when(
    survey == "October19" | survey == "January20" ~ "November19",
    TRUE ~ as.factor(survey)))



suscep1 <- njlong %>% subset(scientific_name == "Meandrina meandrites" |
                               scientific_name == "Eusmilia fastigiata" |
                               scientific_name == "Dichocoenia stokesii" |
                               scientific_name == "Dendrogyra cylindrus") %>%
  group_by(location_name, time_point, survey, scientific_name) %>%
  summarize(N_site = sum(count_adult), density = (sum(count_adult))/30)

susmeans1 <- suscep1 %>% group_by(survey, scientific_name) %>%
  summarize(MeanDens = mean(density), seDens = se(density), maxDens = max(density), N = sum(N_site)) %>%
  unite("GroupEvent", c(scientific_name,survey))

suscep1 <- suscep1 %>% unite("GroupEvent", c(scientific_name, survey), remove = FALSE) %>%
  left_join(susmeans1, by = "GroupEvent") %>%
  mutate_at(.vars = vars("time_point"),
            .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

library(MASS)

library(detectseparation)

#family level df
levels(as.factor(njlong$scientific_name))

efas <- demo %>% subset(scientific_name == "Eusmilia fastigiata")

njlong$scientific_name <- recode(njlong$scientific_name, 
                                 "Porites porities" = "Porites porites")

njfam <- njlong %>% mutate("Family" = case_when(
  scientific_name == "Agaricia spp."|scientific_name == "Agaricia agaricites"|
    scientific_name == "Agaricia tenuifolia"|scientific_name == "Agaricia lamarckiana"|
    scientific_name == "Helioseris cucullata" ~ "Agariciidae",
  scientific_name == "Pseudodiploria strigosa"|scientific_name == "Orbicella faveolata"|
    scientific_name == "Montastraea cavernosa"|scientific_name == "Diploria labyrinthiformis"|
    scientific_name == "Orbicella spp."|scientific_name == "Orbicella annularis"|
    scientific_name == "Orbicella franksi"|scientific_name == "Favia fragum"|
    scientific_name == "Colpophyllia natans"|scientific_name == "Pseudodiploria clivosa" ~ "Faviidae",
  scientific_name == "Porites porites"|scientific_name == "Porites astreoides" ~ "Poritidae",
  scientific_name == "Stephanocoenia intersepta" ~ "Astrocoeniidae",
  scientific_name == "Eusmilia fastigiata"|scientific_name == "Dendrogyra cylindrus"|
    scientific_name == "Meandrina meandrites"|scientific_name == "Dichocoenia stokesii" ~ "Meandrinidae",
  scientific_name == "Siderastrea radians"|scientific_name == "Siderastrea siderea" ~ "Siderastreidae",
  scientific_name == "Mycetophyllia aliciae"|scientific_name == "Isophyllia rigida"|
    scientific_name == "Mycetophyllia spp."|scientific_name == "Scolymia cubensis"|
    scientific_name == "Isophyllia sinuosa" ~ "Mussidae",
  scientific_name == "Madracis decactis" ~ "Pocilloporidae",
  scientific_name == "Acropora cervicornis" ~ "Acroporidae"))

njfam$count_adult[is.na(njfam$count_adult)] <- 0 

mea <- njfam %>% subset(Family == "Meandrinidae")
famdens <- njfam %>% 
  group_by(location_name, survey, event, time_point, Family) %>%
  summarize(count_adult = sum(count_adult)) %>%
  mutate(count = count_adult+1, dens = count_adult/30)

hist(famdens$count)

fammod <- glmer(count ~ survey*Family + (1|location_name), family = "poisson", 
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                data = famdens)

hist(resid(fammod)) #nice
Anova(fammod)

fammod2 <- glmer.nb(count ~ survey*Family + (1|location_name), 
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    data = famdens)
hist(resid(fammod2))
Anova(fammod2)

AIC(fammod)
AIC(fammod2)

emm <- emmeans(fammod2, ~ survey*Family)
simple <- pairs(emm, simple = "survey")
pairwise <- as.data.frame(pairs(emm, simple = "survey"))

df <- data.frame()

Fam <- levels(as.factor(pairwise$Family))

for(current_Fam in Fam) {
  
  Fam_df <- pairwise %>% subset(Family == current_Fam) 
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = Fam_df,
                                   threshold = 0.05)) %>%
    mutate("Family" = current_Fam)
  
  df <- df %>%
    bind_rows(letters)
}



############################
#species specific densities
dens <- njlong %>% rename("Species" = scientific_name)
dens$count_adult[is.na(dens$count_adult)] <- 0 


denssum <- dens %>% group_by(Species) %>% summarize(total = sum(count_adult)) %>%
  arrange(total)

dens <- dens %>% right_join(denssum, by = "Species")

dens <- dens %>% subset(total > 22)


hist(dens$count_adult)

densmod <- glmer(count_adult ~ survey*Species + (1|location_name), family = "poisson", 
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 data = dens)

hist(resid(densmod)) 
Anova(densmod)



densmod2 <- glmer.nb(count_adult ~ survey*Species + (1|location_name), 
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                     data = dens)
hist(resid(densmod2))
Anova(densmod2)

AIC(densmod)
AIC(densmod2)

emm <- emmeans(densmod2, ~ survey*Species)
simple <- pairs(emm, simple = "survey")
pairwise <- as.data.frame(pairs(emm, simple = "survey"))

df <- data.frame()

Spec <- levels(as.factor(pairwise$Species))

for(current_Spec in Spec) {
  
  Spec_df <- pairwise %>% subset(Species == current_Spec) 
  
  pairwise$sig <- isSig(pairwise$p.value)
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = Spec_df,
                                   threshold = 0.05)) %>%
    mutate("Species" = current_Spec)
  
  df <- df %>%
    bind_rows(letters)
}

denslet <- df %>% unite("Event", c("Species", "Group"))

dens <- dens %>% unite("Event", c("Species", "survey"), remove = FALSE)

dens2 <- dens %>% left_join(denslet, by = "Event")

denssum <- dens2 %>% group_by(Event) %>%
  summarize(max = max(count_adult), N = sum(count_adult))

dens3 <- dens2 %>% left_join(denssum, by = "Event")

suscep1p <- ggplot() +
  geom_jitter(data = dens3, aes(x = time_point, y = count_adult, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  #geom_text(data = suscep1, 
  #aes(x = survey, y = maxDens + 0.05, label = Letter)) +
  geom_text(data = dens3, 
            aes(x = survey, y = max*1.4, label = Letter)) +
  geom_text(data = dens3, 
            aes(x = survey, y = max*1.7, label = paste0("N=", N)), size = 2, fontface = "italic") +
  geom_violin(data = dens3, aes(x = survey, y = count_adult), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = dens3, aes(x = time_point, y = count_adult, fill = location_name), size = 2, pch = 21,  
              width = 0.7, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  scale_y_continuous("Count in 30m2", expand = expansion(mult = c(0, 0.2))) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22')) +
  coord_cartesian(xlim = c(-5, 45)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  facet_wrap(~Species, scales = "free") +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

tiff("dens_fullmod_negbinom.tif",width = 9, height = 8, units = "in", res = 400)
suscep1p
dev.off()



suscep2 <- njlong %>% subset(scientific_name == "MCAV" |
                               scientific_name == "PSTR" |
                               scientific_name == "DLAB" |
                               scientific_name == "SSID") %>%
  group_by(location_name, survey, time_point, scientific_name) %>%
  summarize(N_site = sum(count_all), density = (sum(count_all))/30) 

susmeans2 <- suscep2 %>% group_by(survey, scientific_name) %>%
  summarize(MeanDens = mean(density), seDens = se(density), maxDens = max(density), N = sum(N_site)) %>%
  unite("GroupEvent", c(scientific_name,survey))

suscep2 <- suscep2 %>% unite("GroupEvent", c(scientific_name, survey), remove = FALSE) %>%
  left_join(susmeans2, by = "GroupEvent") %>%
  left_join(Letters, by = "GroupEvent") %>%
  mutate(Letter = ifelse(scientific_name == "DLAB",  
                         "", Letter))


suscep2p <- ggplot() +
  geom_jitter(data = suscep2, aes(x = time_point, y = density, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  geom_text(data = suscep2, 
            aes(x = survey, y = maxDens*1.4, label = Letter)) +
  geom_text(data = suscep2, 
            aes(x = survey, y = maxDens*1.6, label = paste0("N=", N)), size = 2, fontface = "italic") +
  geom_violin(data = suscep2, aes(x = survey, y = density), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = suscep2, aes(x = time_point, y = density, fill = location_name), size = 2, pch = 21,  
              width = 0.7, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~scientific_name, nrow = 1, scales = "free") +
  scale_y_continuous("Mean Density/m2", expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22')) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
suscep2p


plotsus1 <- grid.arrange(suscep1p, suscep2p, ncol=1, nrow =2)

plotsus2 <- cowplot::plot_grid(plotsus1, legend, rel_widths = c(6/7, 1/7), axis = 't', align = "v")

tiff("dens_suscep.tif",width = 9, height = 7, units = "in", res = 400)
plotsus2
dev.off()

png("dens_suscep.png", width = 9, height = 7, units = "in", res = 400)
plotsus2
dev.off()

dens <- dens %>% mutate(Category = ifelse(Species == "ACER" | Species == "ATEN" |
                                            Species == "PPOR" | Species == "PAST" |
                                            Species == "AAGA"| Species == "SRAD",
                                          "SCTLD Resistant", "x")) %>%
  mutate(Category = ifelse(Species == "DLAB"|Species == "PSTR"|Species == "CNAT"|
                             Species == "EFAS"|Species == "DSTO"|Species == "MMEA",
                           "Highly SCTLD Susceptible", Category)) %>%
  mutate(Category = ifelse(Species == "MCAV"|Species == "SSID"|Species == "OANN"|
                             Species == "OFAV"|Species == "SINT",
                           "Int. SCTLD Susceptible", Category)) %>%
  subset(Category != "x") %>%
  mutate(density = count/30)



###########
#most species for supplement (those that appear every year and/or aren't high suscep)


densall <- njlong %>% subset(scientific_name != "ALAM" & scientific_name != "CNAT" &
                               scientific_name != "IRIG" & scientific_name != "ISIN" &
                               scientific_name != "MALI" & scientific_name != "OFRA" &
                               scientific_name != "PCLI" & scientific_name != "ACER" &
                               scientific_name != "FFRA" & scientific_name != "SCUB" &
                               scientific_name != "AGAR" & scientific_name != "ORBI") %>%
  group_by(location_name, survey ,time_point, scientific_name) %>%
  summarize(N = sum(count_all), density = (sum(count_all))/30, maxDens = max(density)) 

densall <- densall %>% mutate_at(.vars = vars("time_point"), 
                                 .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


dens_spp_all <- ggplot() +
  geom_jitter(data = densall, aes(x = time_point, y = density, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  #geom_text(data = densall, 
  # aes(x = survey, y = maxDens*1.6, label = paste0("N=", N)), size = 2) +
  geom_violin(data = densall, aes(x = survey, y = density), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = densall, aes(x = time_point, y = density, fill = location_name), size = 2, pch = 21,  
              width = 0.7, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~scientific_name, scales = "free") +
  scale_y_continuous("Mean Density/m2", expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22')) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))


tiff("densplotall.tif", width = 12, height = 9, units = "in", res = 400)
dens_spp_all
dev.off()

png("densplotall.png", width = 12, height = 9, units = "in", res = 400)
dens_spp_all
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

cond$sample_collection_year <- as.factor(cond$sample_collection_year)
levels(cond$sample_collection_year)

cond <- cond %>% subset(sample_collection_year != 0) %>% droplevels()

cond$Year_a <- cond$sample_collection_year
cond$Month_a <- cond$sample_collection_month
cond$Day_a <- cond$sample_collection_day
cond$location_name_a <- cond$location_name


cond <- cond %>% subset(location_name != "Tobacco Reef") %>% droplevels()

cond <- cond %>% unite("date", c("Year_a","Month_a","Day_a"), sep = "-")

cond$date <- as.Date(cond$date)

cond$Year_a <- cond$sample_collection_year
cond$Month_a <- cond$sample_collection_month

cond <- cond %>% unite("time_point", c("Year_a","Month_a"), sep = "-")

levels(as.factor(cond$time_point))

cond$time_point <- recode(cond$time_point, "2019-10" = "October19",
                          "2022-12" = "December22", 
                          "2022-5" = "May22",
                          "2020-1" = "January20")

cond <- cond %>% mutate("survey" = time_point)
cond$survey <- recode(cond$survey, "October19" = "November19",
                      "January20" = "November19")

cond$condition_code <- recode(cond$condition_code, "CLP;CLB" = "CLP/CLB")

cond <- cond %>% subset(condition_code != "")

cond$condition_code_a <- cond$condition_code

cond2 <- cond %>% separate_wider_delim(condition_code_a, "/", names = c("cond_1", "cond_2"),
                                       too_few = "align_start")
condprint <- cond2 %>% dplyr::select(-c(cond_1, cond_2))
write.csv(condprint, "cond19to22.csv")

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

#sum_tl <- cond %>% subset(cond_1 == "TL" &  perc_1 != "A"| cond_2 == "TL" & perc_2 != "A") %>%
# group_by(time_point, location_name, scientific_name) %>%
#summarize(n_tl = n())

sum_tl <- cond %>% subset(cond_1 == "TL"| cond_2 == "TL") %>%
  subset(max_diameter > 4) %>%
  group_by(time_point, survey, location_name, scientific_name) %>%
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
  group_by(time_point, survey, location_name, scientific_name) %>%
  summarize(total_nj = sum(non_juv))

sum_demo$time_point_a <- sum_demo$time_point
sum_demo$location_name_a <- sum_demo$location_name
sum_demo$scientific_name_a <- sum_demo$scientific_name

sum_demo <- sum_demo %>% unite("Event", c("time_point_a", "location_name_a", "scientific_name_a"), sep = "_")

#need to add back the corals that were total recent dead-include in TL

tldf <- sum_demo %>% left_join(sum_tl, by = "Event") %>% 
  dplyr::select(-c(time_point.y:scientific_name.y)) %>%
  mutate(total_nj = ifelse(Event == "May22_South Reef Central_EFAS", 
                           total_nj + 2, total_nj)) %>%
  mutate(total_nj = ifelse(Event == "May22_South Reef Central_PSTR", 
                           total_nj + 1, total_nj)) %>%
  replace(is.na(.), 0) %>%
  mutate(prev_tl = n_tl/total_nj) %>%
  mutate(n_healthy = total_nj - n_tl) %>%
  subset(!(is.na(prev_tl)))


#tldf <- na.omit(tldf)

colnames(tldf) <- k gsub("\\.x","",colnames(tldf))


tldf$survey <- recode(tldf$survey, 
                      "January20" = "November19")

#lump Orbicellas
orbi <- tldf %>% subset(scientific_name == "ORBI"|
                          scientific_name == "OANN"|scientific_name == "OFAV"|
                          scientific_name == "OFRA") %>%
  group_by(time_point, survey, location_name, Event) %>%
  summarize(total_nj = sum(total_nj), n_tl = sum(n_tl), n_healthy = sum(n_healthy)) %>%
  mutate(scientific_name = "ORBI") %>%
  mutate(prev_tl = n_tl/total_nj)

tldf <- tldf %>% subset(scientific_name != "ORBI" &
                          scientific_name != "OANN" & scientific_name != "OFAV" &
                          scientific_name != "OFRA") %>%
  bind_rows(orbi)

#convert to binomial (this didn't help)
#sum_demo_long <- tldf %>% dplyr::select(Event, time_point, location_name, scientific_name, 
# prev_tl, n_healthy) %>%
#uncount(n_healthy) %>% mutate("tl_val"  = 0)

#sum_tl_long <- tldf %>% dplyr::select(Event, time_point, location_name,
#scientific_name, prev_tl, n_tl) %>%
#uncount(n_tl) %>% mutate("tl_val" = 1)

#tl_df_long <- rbind(sum_demo_long, sum_tl_long)

tldf <- tldf %>% mutate(scaled_prevtl = prev_tl)

check <- tldf %>% arrange(prev_tl)

mod <- glm(prev_tl ~ survey*location_name, data = tldf, family = "quasibinomial")
hist(resid(mod))
Anova(mod)
rsquared(mod)


mod3 <- lmer(prev_tl ~ survey + (1|location_name), data = tldf)
hist(resid(mod3))
Anova(mod3)


tldf$scaled_prevtl <- scale(tldf$scaled_prevtl)


mod4 <- glmer(scaled_prevtl ~ survey +
                (1|location_name),
              weights=total_nj,
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=tldf,family="binomial")
hist(resid(mod4))
shapiro.test(resid(mod4))
Anova(mod4)



###########
#most species for supplement (those that appear every year and/or aren't high suscep)

tldf <- tldf %>% mutate_at(.vars = vars("time_point"), 
                           .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


sup <- tldf %>% subset(scientific_name != "ALAM" & scientific_name != "CNAT" &
                         scientific_name != "IRIG" & scientific_name != "ISIN" &
                         scientific_name != "MALI" & scientific_name != "OFRA" &
                         scientific_name != "PCLI" & scientific_name != "ACER" &
                         scientific_name != "AGAR" & scientific_name != "MYCE")

prevp_spp_all <- ggplot() +
  geom_jitter(data = sup, aes(x = time_point, y = prev_tl, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  #geom_text(data = densall, 
  # aes(x = survey, y = maxDens*1.6, label = paste0("N=", N)), size = 2) +
  geom_violin(data = sup, aes(x = survey, y = prev_tl), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = sup, aes(x = time_point, y = prev_tl, fill = location_name), size = 2, pch = 21,  
              width = 0.7, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~scientific_name, ncol = 6, scales = "free") +
  scale_y_continuous("Tissue Loss Prevalence", limits = function(x){c(0, max(0.0001, x))}) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22')) +  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 8, angle = 45),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))


tiff("prevplotall.tif", width = 12, height = 6, units = "in", res = 400)
prevp_spp_all
dev.off()

png("prevplotall.png", width = 12, height = 9, units = "in", res = 400)
prevp_spp_all
dev.off()

##########
#target species

target_prev <- tldf %>% subset(
  scientific_name == "SSID"|
    scientific_name == "PSTR"|scientific_name == "MCAV"|
    scientific_name == "ORBI") %>%
  droplevels() %>% unite("GroupEvent", c(scientific_name, survey), remove = FALSE)

prev_group <- target_prev %>% group_by(survey, scientific_name) %>%
  summarize(MeanPrev = mean(prev_tl), sePrev = se(prev_tl), maxPrev = max(prev_tl)) %>%
  unite("GroupEvent", c(scientific_name,survey))

target_prev <- left_join(target_prev, prev_group, by = "GroupEvent")

#library('brms')

#write.csv(tl_df_long, "binomial_survey_tl.csv")

#binom_prev <- tl_df_long %>% subset(
#scientific_name == "SSID"|
#scientific_name == "PSTR"|scientific_name == "MCAV"|
# scientific_name == "OANN"|scientific_name == "OFAV") %>%
#droplevels() %>% unite("GroupEvent", c(scientific_name, time_point), remove = FALSE)

#TimeLevels = c('Oct19','May22','Dec22')

target_prev <- target_prev %>% mutate_at(.vars = vars("time_point"), 
                                         .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


SpecList <- levels(as.factor(target_prev$scientific_name))

#spec <- binom_prev %>% subset(scientific_name == "MCAV")

#mod <- glmer(tl_val ~ time_point + scientific_name +
#(1|location_name),
#control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
# data=binom_prev, family = "binomial") 

#hist(resid(mod))
#Anova(mod)
#rsquared(mod)


#formula <- bf(tl_val ~ time_point + (time_point|scientific_name))

#brm_optimized <- brm(formula = formula, data = binom_prev, family = bernoulli(),
#chains = 4, cores = min(12, parallel::detectCores()),
#threads = threading(threads=min(12, parallel::detectCores())), iter = 8000,
#warmup = 2000, thin = 1, control = list(adapt_delta = 0.99,
#max_treedepth = 12), save_pars = save_pars(all = TRUE))

#pairwise <- as.data.frame(pairs(emmeans(mod, ~ time_point)))

#pairwise$sig <- isSig(pairwise$p.value)

#letters <- as.data.frame(cldList(p.value ~ contrast,
#data = pairwise,
#threshold = 0.05)) %>%
#mutate("Species" = current_Spec


#####################################
#Individual Species Prevalence Models

#MCAV
mcav <- target_prev %>% subset(scientific_name == "MCAV")

hist(mcav$prev_tl)

tlmod1 <- glm(prev_tl ~ survey*location_name, data = mcav, family = "quasibinomial")
hist(resid(mod1)) #bad
Anova(tlmod1) #does not produce p values

options(contrasts = c("contr.sum","contr.poly"))
tlmod2 <- lmer(prev_tl ~ survey + (1|location_name), data = mcav)
hist(resid(tlmod2))
shapiro.test(resid(tlmod2)) #not normal
Anova(tlmod2) #not significant

tlmod3 <- glmer(prev_tl ~ survey +
                  (1|location_name),
                weights=total_nj,
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                data=mcav,family="binomial") #issingular
hist(resid(tlmod3))
Anova(tlmod3)

tlmod3.1 <- glmer.nb(n_tl ~ survey +
                       (1|location_name),
                     weights=total_nj,
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                     data=mcav) #issingular
hist(resid(tlmod3))
Anova(tlmod3)

tlmod4 <- glm(prev_tl ~ survey, data = mcav, weights = total_nj, family = "quasibinomial")
hist(resid(tlmod4))
Anova(tlmod4)


tlmod4.1 <- glm(prev_tl ~ survey, data = mcav, weights = total_nj, family = "binomial")
hist(resid(tlmod4.1))
Anova(tlmod4.1)

library(PROreg)
#mod4.2 <- BBmm(fixed.formula = prev_tl~survey,random.formula = ~location_name,m=10,data=mcav)
#model

AIC(mod4)
AIC(mod4.1)

mcav_pairwise <- as.data.frame(pairs(emmeans(tlmod4.1, ~ survey)))

mcav_pairwise$sig <- isSig(mcav_pairwise$p.value)

mcav_letters <- as.data.frame(cldList(p.value ~ contrast,
                                      data = mcav_pairwise,
                                      threshold = 0.05)) %>%
  mutate("Species" = "MCAV")

#SSID
ssid <- target_prev %>% subset(scientific_name == "SSID")

hist(ssid$prev_tl)

tlmod5 <- glm(prev_tl ~ survey*location_name, data = ssid, family = "quasibinomial")
hist(resid(tlmod5)) #one big block
Anova(tlmod5) #does not produce p values

tlmod6 <- lmer(prev_tl ~ survey + (1|location_name), data = ssid)
hist(resid(tlmod6)) #not terrible
Anova(tlmod6) #significant

tlmod7 <- glmer(prev_tl ~ survey +
                  (1|location_name),
                weights=total_nj,
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                data=ssid,family="binomial") #issingular
hist(resid(tlmod7))
Anova(tlmod7)

tlmod8 <- glm(prev_tl ~ survey, data = ssid, weights = total_nj, family = "quasibinomial")
hist(resid(tlmod8))
Anova(tlmod8)

tlmod8.1 <- glm(prev_tl ~ survey, data = ssid, weights = total_nj, family = "binomial")
hist(resid(tlmod8.1))
Anova(tlmod8.1)

AIC(tlmod8)
AIC(tlmod8.1)


ssid_pairwise <- as.data.frame(pairs(emmeans(tlmod6, ~ survey)))

ssid_pairwise$sig <- isSig(ssid_pairwise$p.value)

ssid_letters <- as.data.frame(cldList(p.value ~ contrast,
                                      data = ssid_pairwise,
                                      threshold = 0.05)) %>%
  mutate("Species" = "ssid")

#Tissue Loss Loop
df <- data.frame()
mod_df <- data.frame()

for(current_Specie in SpecList) {
  
  Spec_df <- target_prev %>% subset(scientific_name == current_Specie) 
  
  mod <- lmer(prev_tl ~ survey + (1|location_name), data = Spec_df)
  
  modresult <- as.data.frame(Anova(mod)) %>%
    mutate("Species" = current_Specie)
  modresult$sig <- isSig(modresult$"Pr(>Chisq)")
  
  mod_df <- mod_df %>%
    bind_rows(modresult)
  
  pairwise <- as.data.frame(pairs(emmeans(mod, ~ survey)))
  
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
  geom_jitter(data = target_prev, aes(x = time_point, y = prev_tl, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  #geom_text(data = densall, 
  # aes(x = survey, y = maxDens*1.6, label = paste0("N=", N)), size = 2) +
  geom_violin(data = target_prev, aes(x = survey, y = prev_tl), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = target_prev, aes(x = time_point, y = prev_tl, fill = location_name), size = 2, pch = 21,  
              width = 0.7, height = 0) +
  #geom_text(data = target_prev, 
  #aes(x = survey, y = maxPrev, label = Letter), nudge_y = 0.05) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~scientific_name, ncol = 6, scales = "free") +
  scale_y_continuous("Tissue Loss Prevalence", limits = function(x){c(0, max(0.0001, x))}) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22')) +  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 8, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))


#tiff("prevplot.tif", width = 8, height = 4, units = "in", res = 400)
prevp_spp
#dev.off()



#all stony coral

#sto_tl <- sum_tl %>% group_by(time_point, survey, location_name) %>%
#summarize(n_tl = sum(n_tl))

sto_tl <- tldf %>% group_by(time_point, survey, location_name) %>%
  summarize(n_tl = sum(n_tl), total_nj = sum(total_nj), n_healthy = sum(n_healthy)) %>%
  mutate(prev_tl = n_tl/total_nj)

sto_tl$survey_a <- sto_tl$survey
sto_tl$location_name_a <- sto_tl$location_name

sto_tl <- sto_tl %>% unite("Event", "survey_a", "location_name_a", sep = "_")

#sto_demo <- demo %>%
# group_by(time_point, survey ,location_name) %>%
#summarize(total_nj = sum(non_juv))

#sto_demo$time_point_a <- sto_demo$time_point
#sto_demo$location_name_a <- sto_demo$location_name

#sto_demo <- sto_demo %>% unite("Event", "time_point_a", "location_name_a", sep = "_")


#prev_sto <- sto_demo %>% left_join(sto_tl, by = "Event")# %>% 
#dplyr::select(-c(time_point.y:location_name.y)) %>%
#replace(is.na(.), 0) %>%
#mutate(prev_tl = n_tl/total_nj) %>%
prev_sto <- sto_tl %>% mutate(label = "All Stony Coral")

#colnames(prev_sto) <- gsub("\\.x","",colnames(prev_sto))

sto_groups <- prev_sto %>% group_by(survey) %>%
  summarize(MeanPrev = mean(prev_tl), sePrev = se(prev_tl), maxPrev = max(prev_tl))

prev_sto <- left_join(prev_sto, sto_groups, by = "survey")


#TimeLevels = c('Oct19','May22','Dec22')

#prev_sto <- prev_sto %>% mutate_at(.vars = vars("time_point"), 
#.funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

hist(prev_sto$prev_tl)

stomod4 <- lmer(prev_tl ~ survey + (1|location_name), data = prev_sto)
hist(resid(stomod4))
Anova(stomod4)
rsquared(stomod4)

#this model gave weird results, can't be appropriate
stomod2 <- glm(n_tl ~ survey,
               weights=total_nj,
               data=prev_sto,family="poisson") 
hist(resid(stomod2))
shapiro.test(resid(stomod2))
Anova(stomod2)




stomod <- glmer(prev_tl ~ survey +
                  (1|location_name),
                weights=total_nj,
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                data=prev_sto,family="binomial")
hist(resid(stomod))
shapiro.test(resid(stomod)) #good
Anova(stomod)


pairwise <- as.data.frame(pairs(emmeans(stomod, ~ survey)))

pairwise$sig <- isSig(pairwise$p.value)

letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))

prev_sto <- left_join(prev_sto, letters, by = c("survey" = "Group"))

#TimeLevels = c('Oct19','May22','Dec22')

prev_sto <- prev_sto %>% mutate_at(.vars = vars("time_point"), 
                                   .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))


prevp_sto <- ggplot() +
  geom_jitter(data = prev_sto, aes(x = time_point, y = prev_tl, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  geom_violin(data = prev_sto, aes(x = survey, y = prev_tl), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = prev_sto, aes(x = time_point, y = prev_tl, fill = location_name), size = 2, pch = 21,  
              width = 0.7, height = 0) +
  geom_text(data = prev_sto, 
            aes(x = survey, y = maxPrev, label = Letter), nudge_y = 0.05) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_grid(~label, scales = 'free') +
  scale_y_continuous("Tissue Loss Prevalence", limits = function(x){c(0, max(0.0001, x))}) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22')) +  scale_fill_manual("Site",values=c(sitecolors)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

prevp_sto


sums <- prev_sto %>% group_by(survey) %>%
  summarize(n_tl = sum(n_tl), total_nj = sum(total_nj))

#SCTLD suscep stony coral

sus_tl <- tldf %>% subset(scientific_name == "MCAV"|
                            scientific_name == "OANN"|scientific_name == "OFAV"|
                            scientific_name == "PSTR"|scientific_name == "CNAT"|
                            scientific_name == "SINT"|scientific_name == "SSID"|scientific_name == "PCLI"|
                            scientific_name == "DLAB"|scientific_name == "DSTO"|
                            scientific_name == "EFAS"|scientific_name == "MMEA"|scientific_name == "DCYL") %>% droplevels() %>%
  group_by(time_point, survey, location_name) %>%
  summarize(n_tl = sum(n_tl), total_nj = sum(total_nj), n_healthy = sum(n_healthy)) %>%
  mutate(prev_tl = n_tl/total_nj)

sus_tl$survey_a <- sus_tl$survey
sus_tl$location_name_a <- sus_tl$location_name

sus_tl <- sus_tl %>% unite("Event", "survey_a", "location_name_a", sep = "_")

#sus_demo <- demo %>% subset(scientific_name == "MCAV"|
#scientific_name == "OANN"|scientific_name == "OFAV"|
#scientific_name == "PSTR"|scientific_name == "CNAT"|
#scientific_name == "SINT"|scientific_name == "SSID"|scientific_name == "PCLI"|
#scientific_name == "DLAB"|scientific_name == "DSTO"|
#scientific_name == "EFAS"|scientific_name == "MMEA"|scientific_name == "DCYL") %>%
#group_by(time_point, location_name) %>%
#summarize(total_nj = sum(non_juv))

#sus_demo$time_point_a <- sus_demo$time_point
#sus_demo$location_name_a <- sus_demo$location_name

#sus_demo <- sus_demo %>% unite("Event", "time_point_a", "location_name_a", sep = "_")


#prev_sus <- sus_demo %>% left_join(sus_tl, by = "Event") %>% 
#dplyr::select(-c(time_point.y:location_name.y)) %>%
#replace(is.na(.), 0) %>%
#mutate(prev_tl = n_tl/total_nj) %>%
prev_sus <- sus_tl %>% mutate(label = "SCTLD Suscep. Spp.")

#colnames(prev_sus) <- gsub("\\.x","",colnames(prev_sus))

sus_groups <- prev_sus %>% group_by(survey) %>%
  summarize(MeanPrev = mean(prev_tl), sePrev = se(prev_tl), maxPrev = max(prev_tl))

prev_sus <- left_join(prev_sus, sus_groups, by = "survey")

#TimeLevels = c('Oct19','May22','Dec22')

prev_sus <- prev_sus %>% mutate_at(.vars = vars("time_point"), 
                                   .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

hist(prev_sus$prev_tl)
shapiro.test(prev_sus$prev_tl)

susmod <- lmer(prev_tl ~ time_point + (1|location_name), data = prev_sus)
hist(resid(susmod))
Anova(susmod)
rsquared(susmod)

susmod2 <- glmer(prev_tl ~ survey +
                   (1|location_name),
                 weights=total_nj,
                 data=prev_sus,family="binomial")
#hist(resid(susmod2))
#Anova(susmod2)

prev_sus$scaled_prev_tl <- scale(prev_sus$prev_tl)
prev_sus$log_prev_tl <- log(prev_sus$prev_tl + 1)
prev_sus$sqrt_prev_tl <- sqrt(prev_sus$prev_tl)


susmod3 <- glmer(log_prev_tl ~ survey +
                   (1|location_name),
                 weights=total_nj,
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 data=prev_sus,family="binomial")

pairwise <- as.data.frame(pairs(emmeans(susmod, ~ time_point)))

pairwise$sig <- isSig(pairwise$p.value)

letters <- as.data.frame(cldList(p.value ~ contrast,
                                 data = pairwise,
                                 threshold = 0.05))

prev_sus <- left_join(prev_sus, letters, by = c("time_point" = "Group"))


#TimeLevels = c('Oct19','May22','Dec22')

prev_sus <- prev_sus %>% mutate_at(.vars = vars("time_point"), 
                                   .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))



prevp_sus <- ggplot() +
  geom_jitter(data = prev_sus, aes(x = time_point, y = prev_tl, fill = location_name), size = 2, pch = 21,  
              width = 0.25, 
              height = 0, alpha = 0) +
  geom_violin(data = prev_sus, aes(x = survey, y = prev_tl), fill = "gray80", outlier.shape = NA) +
  geom_jitter(data = prev_sus, aes(x = time_point, y = prev_tl, fill = location_name), size = 2, pch = 21,  
              width = 0.7, height = 0) +
  #geom_text(data = prev_sus, 
  #aes(x = survey, y = maxPrev, label = Letter), nudge_y = 0.05) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_grid(~label, scales = 'free') +
  scale_y_continuous("Tissue Loss Prevalence", limits = function(x){c(0, max(0.0001, x))}) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22')) +  scale_fill_manual("Site",values=c(sitecolors)) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

prevp_sus


sums <- prev_sus %>% group_by(survey) %>%
  summarize(n_tl = sum(n_tl), total_nj = sum(total_nj))

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

fate2 <-fate %>% rename("Condition_062019" = "062019_Condition",
                        "Condition_052022" = "052022_Condition",
                        "Condition_122022" = "122022_Condition")

fate2 <- fate2 %>% dplyr::select(Date_InitialTag, NewTagNum, Transect, MaxDiameter,
                                 Species, Condition_062019, Condition_052022,
                                 Condition_122022)

fate2 <- fate2 %>% separate_wider_delim(Date_InitialTag, "/",
                                        names = c("Month_Initial", "Day_Initial", "Year_Initial"))

fate3 <- fate2 %>% mutate(Condition_052022 = if_else(Species == "PAST"&
                                                       Transect == "SR30N" & NewTagNum == "47",
                                                     "Not_Visited", Condition_052022))

fate4 <- fate3 %>% subset(Year_Initial == "19") %>%
  mutate(Condition_122022 = if_else(Condition_052022 == "Dead",
                                    "Dead", Condition_122022)) %>%
  subset(Condition_052022 != "Not_Visited") %>%
  subset(Species != "CNAT")

rec <- fate4 %>% subset(Condition_052022 == "Diseased" &
                          Condition_122022 == "Healthy"|Condition_122022 == "Healthy?")

checkSSID <- fate4 %>% subset(Species == "SSID")

Long <-melt(fate4, id.vars =c("Species","NewTagNum",
                              "Transect"), measure.vars = 
              c("Condition_052022", "Condition_122022", "Condition_062019")) %>%
  rename("Time_Point" = variable, "Condition" = value) 

Long$Condition <- recode(Long$Condition, "Healthy?" = "Healthy")

Long$Time_Point <- recode(Long$Time_Point, "Condition_052022" = "May22",
                          "Condition_122022" = "Dec22",
                          "Condition_062019" = "Oct19")


SpecLevels = c('SSID','MCAV','PSTR','PAST','MMEA','CNAT')

Long <- Long %>% mutate_at(.vars = vars("Time_Point"), 
                           .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
  mutate_at(.vars = vars("Species"),
            .funs = funs(factor(.,levels = SpecLevels, ordered = TRUE)))


sumsimple <- Long %>% group_by(Time_Point,Species) %>% count(Condition) %>% ungroup()

condcolors2 = c('Dead'='coral3','Diseased'='gold1','Healthy'='springgreen4',
                'Increased Old Mortality'='gray60',
                'Not Found' = 'black')


bysp <- ggplot(sumsimple, aes(x = Time_Point, y = n, fill = Condition)) +
  facet_wrap(~Species, as.table = FALSE) +
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
        axis.title.x = element_blank(),
        legend.position = c(0.8,0.8))

tiff("ConditionbySp.tif",width = 6, height = 6, units = "in", res = 300)
bysp
dev.off()

png("ConditionbySp.png",width = 6, height = 6, units = "in", res = 300)
bysp
dev.off()

fate4$Condition_122022 <- recode(fate4$Condition_122022, 
                                 "Healthy?" = "Healthy")

td_df <- fate4 %>%
  subset(Condition_052022 == "Diseased"|Condition_052022 == "Dead"|
           Condition_122022 == "Diseased"|Condition_122022 == "Dead") %>%
  subset(!(is.na(MaxDiameter))) %>%
  mutate(timedead = ifelse(Condition_052022 == "Dead",
                           1, 'timedead')) %>%
  mutate(timedead = ifelse((Condition_052022 == "Healthy"|
                              Condition_052022 == "Diseased") &
                             Condition_122022 == "Dead",
                           2, timedead)) %>%
  mutate(timedead = ifelse(Condition_122022 == "Healthy"|
                             Condition_122022 == "Diseased",
                           3, timedead))

td_df <- td_df %>% mutate(Species = fct_relevel(Species,
                                                "MMEA", "SSID", "PAST", "MCAV", "PSTR")) 


td_df$timedead <- as.numeric(td_df$timedead)

model <- glm(timedead ~ MaxDiameter + Species, family=quasibinomial(link='logit'),data=td_df)
hist(resid(model))
Anova(model)
summary(model)


sizep <- ggplot() +
  geom_point(data = td_df, aes(x = MaxDiameter, y = timedead, fill = Species), pch = 21, size = 3, alpha = 0.5) +
  facet_wrap(~Species) +
  scale_y_continuous("Time pds to death", breaks = c(1,2,3)) +
  #scale_x_discrete("") +
  #scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
sizep

#SSID 3 (CBC30N) and SSID 36 (Lagoon) both were diseased in May and healthy in Dec
#photo evidence

remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)



counts <- fate4 %>% group_by(Species) %>% mutate(count = n()) %>%
  select(Species, count) %>% unique()

SpecList <- levels(as.factor(fate4$Species))



sankeydf <- data.frame()

for(current_Specie in SpecList) {
  
  Spec_df <- fate4 %>% subset(Species == current_Specie) 
  
  sankdf <- Spec_df %>% 
    make_long(Condition_062019, Condition_052022, Condition_122022) %>%
    mutate("Species" = current_Specie) 
  
  sankeydf <- sankeydf %>%
    bind_rows(sankdf)
}


sankeydf <- sankeydf %>% left_join(counts, by = "Species")

sankeydf$x <- recode(sankeydf$x, "Condition_062019" = "Oct19",
                     "Condition_052022" = "May22",
                     "Condition_122022" = "Dec22")

sankeydf <- sankeydf %>% mutate(Species = fct_relevel(Species,
                                                      "SSID", "PAST", "MMEA", "MCAV", "PSTR")) 

sankey <- ggplot(sankeydf, aes(x = x, 
                               next_x = next_x, 
                               node = node, 
                               next_node = next_node,
                               fill = factor(node))) +
  facet_wrap(~Species, as.table = FALSE) +
  geom_sankey(flow.alpha = 0.6, node.color = 'black', flow.color = 'black') +
  geom_sankey_label(
    aes(
      x = as.numeric(x) - 0.2,
      label = after_stat(freq)),
    size = 7 / .pt, color = "black", fill = "white") +  
  scale_fill_manual("Condition", values = c(condcolors2)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        strip.text = element_text(size = 13),
        axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 12),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.position = c(0.8,0.8),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
sankey

png("ConditionSankey.png",width = 9, height = 6, units = "in", res = 300)
sankey
dev.off()
###############3
#retired code from long from date_diseased version


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


Long <-melt(fate6, id.vars =c("Species","OldTagNum","NewTagNum",
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
  subset(Transect == "CBC30N") %>% arrange(OldTagNum)

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

#make corrections (based on photos and notes from fate script V2 - double check)
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

write.csv(fate, "CBC_Fate_2019to2022.csv")

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


##############################################
#Can't Use

#Tissue Loss Loop
df <- data.frame()
mod_df <- data.frame()

for(current_Specie in SpecList) {
  
  Spec_df <- target_prev %>% subset(scientific_name == current_Specie) 
  
  mod <- glmer(prev_tl ~ time_point +
                 (1|location_name),
               weights=total_nj,
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
               data=Spec_df,family="binomial") 
  
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
