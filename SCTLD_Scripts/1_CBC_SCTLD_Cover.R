library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)


#To do


# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/GitHub/SCTLD_demographics/")


# Import the survey data
cover <- read.csv("annotations_CBC_19-22.csv")

pd <- position_dodge(width = 0.93)
se<-function(x)sqrt(var(x)/length(x))

isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}

# Establish time levels
TimeLevels <- c("October19","November19","December19",
                "January20", "February20", "March20" ,"April20","May20", "June20","July20",
                "August20" ,"September20", "October20","November20","December20",
                "January21","February21","March21" ,"April21","May21", "June21",
                "July21","August21" ,"September21", "October21","November21","December21",
                "January22","February22","March22" ,"April22","May22", "June22","July22",
                "August22" ,"September22", "October22","November22","December22")


#clean up dates & site names ----
split <- strsplit(cover$Date, "/")
split2 <- matrix(unlist(split), ncol=3, byrow=TRUE)
split2 <- as.data.frame(split2)
cover <- cbind(cover, split2)

cover <- cover %>% rename("Year" = "V3", "Day" = "V2", "Month" = "V1") %>%
  select(Name, Date, Habitat, SiteName, Year, Month, Day, Label, Row, Column)

npics <- cover %>% group_by(Date, SiteName) %>%
  summarize(n = length(unique(Name)))

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

#calculate cover ----
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

#assign general cover label categories ----
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

levels(as.factor(cover$Label_General))




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


#mod5 <- glmmTMB(cov_prop ~ TimePoint +
# (1|SiteName),
# data=covsto,beta_family())
#hist(resid(mod4))
#Anova(mod4)

#########
#model general cover categories ----
#########
mod1 <- glmer(cov_prop ~ Survey*Label_General +
                (1|SiteName),
              weights=npoints,
              data=covgen,family="binomial")
hist(resid(mod1))
Anova(mod1)
AIC(mod1)


gentable <- as.data.frame(Anova(mod1)) %>%
  mutate(`Pr(>Chisq)` = ifelse(`Pr(>Chisq)` < 0.001 ,
                            "<0.001", `Pr(>Chisq)` )) %>%
  mutate(sig = isSig(`Pr(>Chisq)`)) %>%
  kbl(caption = "<span style='color: black;'> <b>Table 6.</b> Model results testing for
      variation in cover of major benthic functional groups over time.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(gentable, file = "Tables/Table 6.pdf")


emm <- emmeans(mod1, ~ Survey*Label_General)
simple <- pairs(emm, simple = "Survey")
pairwise <- as.data.frame(pairs(emm, simple = "Survey"))

eff <- as.data.frame(eff_size(simple, sigma = 26, edf = Inf)) 


df <- data.frame()

Labels <- levels(as.factor(pairwise$Label_General))

for(current_Label in Labels) {
  
  Lab_df <- pairwise %>% subset(Label_General == current_Label) 
  
  letters <- as.data.frame(cldList(p.value ~ contrast,
                                   data = Lab_df,
                                   threshold = 0.05)) %>%
    mutate("Label_General" = current_Label)
  
  df <- df %>%
    bind_rows(letters)
}

#make tables
library(kableExtra)

pairwise <- pairwise %>%
mutate(`p.value` =round(`p.value`, digits = 3)) %>%
  mutate(`p.value` = ifelse(`p.value` < 0.001 ,
                            "<0.001", `p.value` )) %>%
  mutate(sig = isSig(`p.value`))


gen_emmtable <- pairwise %>%
  kbl(caption = "<span style='color: black;'> <b>Table S11.</b> Pairwise contrasts - cover of benthic
  functional groups over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(gen_emmtable, file = "Tables/Table S11.pdf")


gen_efftable <- eff %>%
  kbl(caption = "<span style='color: black;'> <b>Table S12.</b> Effect sizes - cover of benthic functional
  groups over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(gen_efftable, file = "Tables/Table S12.pdf")




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

# ##model susceptible coral cover ----
# modsus <- glm(cov_prop ~ Survey + SiteName, data = suscep, family = "quasibinomial")
# modsus <- glmer(cov_prop ~ Survey +
#                   (1|SiteName),
#                 weights=npoints,
#                 data=suscep,family="binomial")
# hist(resid(modsus))
# Anova(modsus)
# rsquared(modsus)
# 
# susresult <- as.data.frame(Anova(modsus)) %>%
#   mutate(Label_General = "SCTLD Susceptible Coral")
# susresult$sig <- isSig(susresult$"Pr(>Chisq)")
# mod_df <- data.frame()
# mod_df <- mod_df %>% rbind(susresult)
# 
# covtable <- mod_df %>%
#   remove_rownames %>% column_to_rownames(var="Label_General") %>%
#   select(-c(sig)) %>%
#   rename('p-value' = 'Pr(>Chisq)') %>%
#   mutate('p-value' = '<0.001') %>%
#   kbl() %>%
#   kable_styling()
# 
# #png("Tables/CovTable.png", width = 8, height = 6, units = "in", res = 300)
# covtable
# #dev.off()
# 
# suspairwise <- as.data.frame(pairs(emmeans(modsus, ~ Survey)))
# 
# suspairwise$sig <- isSig(suspairwise$p.value)
# 
# susletters <- as.data.frame(cldList(p.value ~ contrast,
#                                     data = suspairwise,
#                                     threshold = 0.05))
# 
# suscep <- suscep %>% left_join(susletters, by = c("Survey" = "Group")) %>%
#   mutate(TimePoint = case_when(
#     SiteName == "SR30N" & Survey == "November19" ~ "January20",
#     SiteName == "CBC30N" & Survey == "November19" ~ "January20",
#     SiteName == "CBC Central" & Survey == "November19" ~ "October19",
#     SiteName == "South Reef Central" & Survey == "November19" ~ "October19",
#     SiteName == "House Reef" & Survey == "November19" ~ "October19",
#     SiteName == "Curlew Patch" & Survey == "November19" ~ "October19",
#     SiteName == "CBC Lagoon" & Survey == "November19" ~ "October19",
#     TRUE ~ as.factor(Survey))) %>%
#   mutate_at(.vars = vars("TimePoint"),
#             .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

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
  geom_boxplot(data = covgen1, aes(x = Survey, y = cov_gen), fill = "gray80",width = 4, outlier.shape = NA) +
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
  geom_boxplot(data = covgen2, aes(x = Survey, y = cov_gen), fill = "gray80", width = 4, outlier.shape = NA) +
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
# 
# suscepyearp <- ggplot() +
#   geom_jitter(data = suscep, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 4, alpha = 0, pch = 21,  
#               width = 0.25, 
#               height = 0, alpha = 0) +
#   geom_text(data = suscep, 
#             aes(x = Survey, y = maxCov, label = Letter), nudge_y = 1) +
#   geom_violin(data = suscep, aes(x = Survey, y = cov_gen), fill = "gray80", outlier.shape = NA) +
#   geom_jitter(data = suscep, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 4, alpha = 0.5, pch = 21,  
#               width = 0.25, 
#               height = 0) +
#   geom_vline(xintercept = "July21", 
#              color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
#   facet_wrap(~Label_Suscep, scales = "free") +
#   scale_y_continuous("") +
#   #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
#   scale_x_discrete("", drop = FALSE, breaks = c('October19','January20','July21',
#                                                 'May22','December22'),
#                    expand = expansion(add = c(2, 2))) +   
#   scale_fill_manual("Site",values=c(sitecolors)) +
#   theme(plot.title = element_text(size = 16,hjust = 0.5),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
#         axis.title = element_text(size = 10),
#         legend.position = 'none',
#         plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
# suscepyearp


legendp <- ggplot() +
  geom_boxplot(data = covgen1, aes(x = TimePoint, y = cov_gen), alpha = 0) +
  #geom_point(data = covgen, aes(x = TimePoint, y = cov_gen, fill = SiteName)) +
  geom_jitter(data = covgen1, aes(x = TimePoint, y = cov_gen, fill = SiteName), size = 2, pch = 21,  
              width = 0.25, 
              height = 0) +
  facet_wrap(~Label_General, scales = "free") +
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


#tiff("Figures/current/CBCCover_General.tif",width = 7, height = 7, units = "in", res = 300)
plot2
#dev.off()

png("Figures/current/Fig 7.png", width = 7, height = 7, units = "in", res = 300)
plot2
dev.off()

#Mean cover values over time ----
covsum <- covgen %>% group_by(Survey, Label_General) %>%
  summarize(mean_cov = mean(cov_gen))

#species specific cover ----
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
  subset(Category != "x") #%>%
#subset(Label != "ACER" & Label != "CNAT" & Label != "DLAB" & Label != "SRAD" & Label != "SINT")




sto2 <- sto %>% rename("Species" = Label) %>%
  mutate(cov_prop = cover/100) %>%
  subset(Species != "ACER" & Species != "CNAT" &
           Species != "DLAB" & Species != "MALC" &
           Species != "MCOM" & 
           Species != "SRAD") %>% droplevels()

orbi <- sto2 %>% subset(Species == "ORBI"|Species == "OFAV"|
                          Species == "OANN"|Species == "OANN") %>%
  group_by(Habitat, SiteName, Date, npoints, TimePoint, 
           Label_General, Label_Suscep, Survey, Category) %>%
  summarize(n = sum(n)) %>%
  mutate("Species" = "ORBI", cov_prop = n/npoints, cover = cov_prop*100)

sto2 <- sto2 %>% rbind(orbi) %>% subset(Species != "OFAV" & Species != "OANN") %>% droplevels()

SpecList <- levels(as.factor(sto2$Species))

#models species specific cover ----
mod3 <- glmer(cov_prop ~ Survey*Species +
                (1|SiteName),
              weights=npoints,
              data=sto2,family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
hist(resid(mod3))
Anova(mod3)

spptable <- as.data.frame(Anova(mod3)) %>%
  mutate(`Pr(>Chisq)` = ifelse(`Pr(>Chisq)` < 0.001 ,
                            "<0.001", `Pr(>Chisq)` )) %>%
  mutate(sig = isSig(`Pr(>Chisq)`)) %>%
  kbl(caption = "<span style='color: black;'> <b>Table 5.</b> Model results testing for
      variation in cover of target scleractinian coral species over time.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(spptable, file = "Tables/Table 5.pdf")


emm <- emmeans(mod3, ~ Survey*Species)
simple <- pairs(emm, simple = "Survey")
pairwise <- as.data.frame(pairs(emm, simple = "Survey"))

eff <- as.data.frame(eff_size(simple, sigma = 26, edf = Inf)) 


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

#make tables
pairwise <- pairwise %>%
  mutate(`p.value` =round(`p.value`, digits = 3)) %>%
  mutate(`p.value` = ifelse(`p.value` < 0.001 ,
                            "<0.001", `p.value` )) %>%
  mutate(sig = isSig(`p.value`))


spp_emmtable <- pairwise %>%
  kbl(caption = "<span style='color: black;'> <b>Table S9.</b> Pairwise contrasts - cover of scleractinian coral
  species over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(spp_emmtable, file = "Tables/Table S9.pdf")


spp_efftable <- eff %>%
  kbl(caption = "<span style='color: black;'> <b>Table S10.</b> Effect sizes - cover of coral species over time.
      Note that the October 2019/January 2020 timepoint is represented as November 2019.<span>",
      digits = 3,
      format = "html", booktabs = TRUE, longtable = TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))
save_kable(spp_efftable, file = "Tables/Table S10.pdf")

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
  mutate(Letter = ifelse(Species == "AAGA" |
                           Species == "PPOR" | 
                           Species == "SINT",
                         "", Letter))

zeros <- read.csv("add_zero_cover.csv")

zeros$TimePoint <- recode(zeros$TimePoint,
                          "Dec-22" = "December22",
                          "May-22" = "May22",
                          "Oct-19" = "October19",
                          "Jan-20" = "January20")

zeros$Survey <- recode(zeros$Survey,
                          "Dec-22" = "December22",
                          "May-22" = "May22",
                          "Nov-19" = "November19")


sto3 <- rbind(sto3, zeros) %>% mutate_at(.vars = vars("TimePoint"),
                                       .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE)))

sto3$Species <- recode(sto3$Species,
                       "AAGA" = "Agaricia agaricites",
                       "ATEN" = "Agaricia tenuifolia",
                       "MCAV" = "Montastraea cavernosa",
                       "PAST" = "Porites astreoides",
                       "PPOR" = "Porites porites",
                       "PSTR" = "Pseudodiploria strigosa",
                       "SINT" = "Stephanocoenia intersepta",
                       "SSID" = "Siderastrea siderea",
                       "ORBI" = "Orbicella spp.",
                       "DSTO" = "Dichocoenia stokesii",
                       "EFAS" = "Eusmilia fastigiata",
                       "MMEA" = "Meandrina meandrites",
                       "DLAB" = "Diploria labyrinthiformis")

library(ggh4x)

highcov <- sto3 %>% subset(Species == "Agaricia tenuifolia" | Species == "Orbicella spp."|
                          Species == "Siderastrea siderea")

highcovp <- ggplot() +
  geom_jitter(data = highcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21, alpha=0,  
              width = 0.25, 
              height = 0) +
  geom_text(data = highcov, 
            aes(x = Survey, y = maxCov+2, label = Letter)) +
  geom_boxplot(data = highcov, aes(x = Survey, y = cover), fill = "gray80", width = 3, outlier.shape = NA) +
  geom_jitter(data = highcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21, alpha=0.5,  
              width = 0.25, 
              height = 0)+
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~ Species, nrow = 1, scales = "free") +
  scale_y_continuous("Percent Cover", limits = c(0,22)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20','July21',
                                                'May22','December22'),
                   expand = expansion(add = c(4, 4))) +   
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

highcovp

medcov <- sto3 %>% subset(Species == "Porites astreoides" | Species == "Porites porites"|
                         Species == "Agaricia agaricites")

medcovp <- ggplot() +
  geom_jitter(data = medcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21, alpha=0,  
              width = 0.25, 
              height = 0) +
  geom_text(data = medcov, 
            aes(x = Survey, y = maxCov+0.4, label = Letter), nudge_y = 0.1) +
  geom_boxplot(data = medcov, aes(x = Survey, y = cover), fill = "gray80", width = 3, outlier.shape = NA) +
  geom_jitter(data = medcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21, alpha=0.5,  
              width = 0.25, 
              height = 0)+
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~ Species, nrow = 1, scales = "free") +
  scale_y_continuous("Percent Cover", limits = c(0,4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20','July21',
                                                'May22','December22'),
                   expand = expansion(add = c(4, 4))) +   
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))
medcovp


lowcov <- sto3 %>% subset(Species == "Montastraea cavernosa" | Species == "Pseudodiploria strigosa"|
                            Species == "Stephanocoenia intersepta")

lowcovp <- ggplot() +
  geom_jitter(data = lowcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21, alpha=0,  
              width = 0.25, 
              height = 0) +
  geom_text(data = lowcov, 
            aes(x = Survey, y = maxCov+0.25, label = Letter), nudge_y = 0.1) +
  geom_boxplot(data = lowcov, aes(x = Survey, y = cover), fill = "gray80", width = 3, outlier.shape = NA) +
  geom_jitter(data = lowcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21, alpha=0.5,  
              width = 0.25, 
              height = 0)+
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~ Species, nrow = 1, scales = "free") +
  scale_y_continuous("Percent Cover", limits = c(0,2.1)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20','July21',
                                                'May22','December22'),
                   expand = expansion(add = c(4, 4))) +   
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

lowcovp


vlowcov <- sto3 %>% subset(Species == "Diploria labyrinthiformis" | Species == "Meandrina meandrites"|
                            Species == "Eusmilia fastigiata" | Species == "Dichocoenia stokesii")

vlowcovp <- ggplot() +
  geom_jitter(data = vlowcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21, alpha=0,  
              width = 0.25, 
              height = 0) +
  geom_text(data = vlowcov, 
            aes(x = Survey, y = maxCov*1.1, label = Letter), nudge_y = 0.1) +
  geom_boxplot(data = vlowcov, aes(x = Survey, y = cover), fill = "gray80", width = 3, outlier.shape = NA) +
  geom_jitter(data = vlowcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 2, pch = 21, alpha=0.5,  
              width = 0.25, 
              height = 0)+
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~ Species, nrow = 1, scales = "free") +
  scale_y_continuous("Percent Cover", limits = c(0,2)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20','July21',
                                                'May22','December22'),
                   expand = expansion(add = c(4, 4))) +   
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"))

vlowcovp

legendp <- ggplot() +
  geom_jitter(data = medcov, aes(x = TimePoint, y = cover, fill = SiteName), size = 4, pch = 21, alpha=0.5,  
              width = 0.25, 
              height = 0) +
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
legendp

library(gridExtra)
library(cowplot)

legend <- cowplot::get_legend(legendp)

covplot <- grid.arrange(highcovp, medcovp, lowcovp, ncol=1, nrow =3)

covplot2 <- cowplot::plot_grid(covplot, legend, rel_widths = c(3/4, 1/4), axis = 't', align = "v")

covplot2

covplot3 <- cowplot::plot_grid(covplot2, vlowcovp, rel_heights = c(7/10,3/10),
                               ncol = 1, axis = 't', align = "v")

covplot3

#tiff("Figures/current/spec_cover_plot.tif",width = 9, height = 8, units = "in", res = 400)
covplot3
#dev.off()

png("Figures/current/Fig 6.png", width = 9, height = 8, units = "in", res = 400)
covplot3
dev.off()
