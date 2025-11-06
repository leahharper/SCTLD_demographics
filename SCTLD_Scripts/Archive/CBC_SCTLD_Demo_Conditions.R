library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)


#To do


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




############################
#DEMO
############################
demo <- read.csv("CBC_Coral_Demographics.csv")
cond <- read.csv("CBC_Coral_Conditions.csv")

demo <- demo %>% subset(year < 2023)

dead2 <- cond %>% subset(percent_mortality == 100)
#3 EFAS and 1 PSTR at SR Central were totally dead!
#PSTR = 6cm
#EFAS at 6, 10, and 13cm

#this code removes dead colonies from demo
demo <- demo %>%
  mutate(X5_to_10_cm = ifelse(diver == "Leah Harper" &
                                month == 5 &
                                year == 2022 & 
                                location_name == "South Reef Central" &
                                scientific_name == "EFAS", 
                              X5_to_10_cm - 2, X5_to_10_cm)) %>%
  mutate(X11_to_20_cm = ifelse(diver == "Leah Harper" &
                                 month == 5 &
                                 year == 2022 & 
                                 location_name == "South Reef Central" &
                                 scientific_name == "EFAS", 
                               X11_to_20_cm - 1, X11_to_20_cm)) %>%
  mutate(X5_to_10_cm = ifelse(diver == "Leah Harper" &
                                month == 5 &
                                year == 2022 & 
                                location_name == "South Reef Central" &
                                scientific_name == "PSTR", 
                              X5_to_10_cm - 1, X5_to_10_cm))


demo <- demo %>%
  mutate(total_small = rowSums(.[17:19])) %>% #5-40cm
  mutate(total_lg = rowSums(.[20:21])) %>% #>40cm
  mutate(adult = rowSums(.[17:21])) #>5cm

levels(as.factor(demo$location_name))


sitecolors = c('CBC Central'='red3','CBC Lagoon'='darkorchid','CBC30N'='gold1',
               'Curlew Patch' = 'blue', 'House Reef' ='seagreen4', 
               'South Reef Central' = 'orange', 'SR30N' = 'pink')

demo$location_name <- recode(demo$location_name,
                             "CBC 30 North" = "CBC30N",
                             "South Reef 30 North" = "SR30N")

demo$scientific_name <- recode(demo$scientific_name, "Montrastraea cavernosa" = "Montastraea cavernosa")
#make all orbicellas one category
demo$code <- gsub("ofav|oann|ofra", "orbi", demo$code)
demo$scientific_name <- gsub("Orbicella faveolata|Orbicella annularis|Orbicella franksi", "Orbicella spp.", demo$scientific_name)


demo$year <- as.factor(demo$year)
levels(demo$year)

demo$Year_a <- demo$year
demo$Month_a <- demo$month
demo$Day_a <- demo$day
demo$location_name_a <- demo$location_name


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

write.csv(demo, "CBC_demo_curated.csv")


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

library(kableExtra)

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
library(lme4)
library(detectseparation)

############################
#species specific densities
dens <- njlong %>% rename("Species" = scientific_name)
dens$count_adult[is.na(dens$count_adult)] <- 0 


denssum <- dens %>% group_by(Species) %>% summarize(total = sum(count_adult)) %>%
  arrange(total)

dens <- dens %>% right_join(denssum, by = "Species") %>% subset(Species != "Agaricia spp.") %>% droplevels()

dens <- dens %>% subset(total > 7) %>% mutate_at(.vars = vars("time_point"),
                                                 .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
  mutate(count1 = count_adult + 1)


hist(dens$count_adult)

densmod <- glmer(count_adult ~ survey*Species + (1|location_name), family = "poisson", 
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 data = dens)

hist(resid(densmod)) 
Anova(densmod)
summary(densmod)
fixed_effects_coefficients <- summary(densmod)$coefficients
fixed_effects_coefficients

library(GLMMadaptive)
#remotes::install_github("glmmTMB/glmmTMB/glmmTMB")



nbinommod <- mixed_model(fixed = count_adult ~ survey * Species, random = ~1 | location_name, 
                         data = dens, family = negative.binomial(), 
                         max_coef_value = 30,
                         control = list(iter_EM = 0))
hist(resid(nbinommod))


zibinom <- mixed_model(fixed = count_adult ~ survey * Species, random = ~1 | location_name, 
                       data = dens, family = zi.negative.binomial(), 
                       zi_fixed = ~ Species, zi_random = NULL,
                       max_coef_value = 30,
                       control = list(iter_EM = 0))
hist(resid(zibinom))

AIC(zibinom)
AIC(nbinommod)

#densmod2 <- glmer.nb(count_adult ~ survey*Species + (1|location_name), 
#control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
#data = dens)
#hist(resid(densmod2))
#Anova(densmod2)

#AIC(densmod)
#AIC(densmod2)

sd(dens$count_adult)

library(emmeans)
library(rcompanion)


isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}

emm <- emmeans(zibinom, ~ survey*Species)
simple <- pairs(emm, simple = "survey")
pairwise <- as.data.frame(pairs(emm, simple = "survey"))

eff <- as.data.frame(eff_size(emm, sigma = 26, edf = Inf))


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

high <- dens3 %>% subset(Species == "Agaricia agaricites" | Species == "Agaricia tenuifolia"|
                           Species == "Porites astreoides" | Species == "Porites porites")

highdensp <- ggplot() +
  geom_jitter(data = high, aes(x = time_point, y = count_adult, fill = location_name), size = 2, pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  geom_text(data = high, aes(x = survey, y = max+10, label = Letter)) +
  geom_boxplot(data = high, aes(x = survey, y = count_adult), fill = "gray80", width = 7) +
  geom_jitter(data = high, aes(x = time_point, y = count_adult, fill = location_name), size = 4, alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~Species, nrow = 1, scales = "free") +
  scale_y_continuous("Colonies in 30m2 transect", limits = c(0,150), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

png("Figures/zibinom_highdens.png",width = 9, height = 4, units = "in", res = 400)
highdensp
dev.off()


med <- dens3 %>% subset(Species == "Orbicella spp." | Species == "Siderastrea radians"|
                          Species == "Siderastrea siderea")

meddensp <- ggplot() +
  geom_jitter(data = med, aes(x = time_point, y = count_adult, fill = location_name), size = 2, pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  geom_text(data = med, aes(x = survey, y = max+10, label = Letter)) +
  geom_boxplot(data = med, aes(x = survey, y = count_adult), fill = "gray80", width = 7) +
  geom_jitter(data = med, aes(x = time_point, y = count_adult, fill = location_name), size = 4, alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~Species, nrow = 1, scales = "free") +
  scale_y_continuous("Colonies in 30m2 transect", limits = c(0,105), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

png("Figures/zibinom_meddens.png",width = 9, height = 4, units = "in", res = 400)
meddensp
dev.off()


low <- dens3 %>% subset(Species == "Diploria labyinthiformis" | Species == "Montastraea cavernosa"|
                          Species == "Pseudodiploria strigosa" | Species == "Stephanocoenia intersepta" |
                          Species == "Meandrina meandrites" | Species == "Eusmilia fastigiata" |
                          Species == "Dichocoenia stokesii")

lowdensp <- ggplot() +
  geom_jitter(data = low, aes(x = time_point, y = count_adult, fill = location_name), size = 2, pch = 21,  
              width = 2, 
              height = 0, alpha = 0) +
  geom_text(data = low, aes(x = survey, y = max+2, label = Letter)) +
  geom_boxplot(data = low, aes(x = survey, y = count_adult), fill = "gray80", width = 7) +
  geom_jitter(data = low, aes(x = time_point, y = count_adult, fill = location_name), size = 4, alpha = 0.5, pch = 21,  
              width = 2, height = 0) +
  geom_vline(xintercept = "July21", 
             color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~Species, nrow = 1, scales = "free") +
  scale_y_continuous("Colonies in 30m2 transect", limits = c(0,20), expand = expansion(mult = c(0, 0.1))) +
  #scale_x_discrete("", drop = FALSE, breaks = every_nth(n=4)) +
  scale_x_discrete("", drop = FALSE, breaks = c('October19','January20', 'July21',
                                                'May22','December22'),
                   expand = expansion(add = c(5, 5))) + 
  scale_fill_manual("Site",values=c(sitecolors)) +
  theme(plot.title = element_text(size = 16,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", hjust = 1, size = 10, angle = 45),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

png("Figures/zibinom_lowdens.png",width = 9, height = 4, units = "in", res = 400)
lowdensp
dev.off()




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

#Conditions
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

demo <- demo %>% subset(year < 2023)
cond <- cond %>% subset(year < 2023)


#################
#Conditions
#################
levels(as.factor(cond$code))

aata <- cond %>% subset(code == "sdstr")

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

sum_tl <- sum_tl %>% unite("Event", c("survey_a", "location_name_a", "scientific_name_a"), sep = "_")

demo <- demo %>% 
  mutate(survey = case_when(
    survey == "October19" | survey == "January20" ~ "November19",
    TRUE ~ as.factor(survey)))

sum_demo <- demo %>%
  group_by(time_point, survey, location_name, scientific_name) %>%
  summarize(Total = sum(Total))

sum_demo$survey_a <- sum_demo$survey
sum_demo$location_name_a <- sum_demo$location_name
sum_demo$scientific_name_a <- sum_demo$scientific_name

sum_demo <- sum_demo %>% unite("Event", c("survey_a", "location_name_a", "scientific_name_a"), sep = "_")


tldf <- sum_demo %>% left_join(sum_tl, by = "Event") %>% 
  dplyr::select(-c(time_point.y:scientific_name.y)) %>%
  replace(is.na(.), 0) %>%
  mutate(prev_tl = n_tl/Total) %>%
  mutate(n_healthy = Total - n_tl) %>%
  subset(!(is.na(prev_tl)))


colnames(tldf) <- gsub("\\.x","",colnames(tldf))


tldf$survey <- recode(tldf$survey, 
                      "January20" = "November19")





tldf$scientific_name <- as.factor(tldf$scientific_name)

sum <- tldf %>% group_by(scientific_name) %>%
  summarize(total_tl = sum(n_tl)) %>%
  arrange(total_tl)

zero_tl <- sum %>% subset(total_tl < 1)

tldf_sub <- tldf %>% subset(!(scientific_name %in% zero_tl$scientific_name))

#convert to binomial (this didn't help)
sum_demo_long <- tldf_sub %>% dplyr::select(Event, survey, location_name, scientific_name, 
                                            prev_tl, n_healthy) %>%
  uncount(n_healthy) %>% mutate("tl_val"  = 0)

sum_tl_long <- tldf_sub %>% dplyr::select(Event, survey, location_name,
                                          scientific_name, prev_tl, n_tl) %>%
  uncount(n_tl) %>% mutate("tl_val" = 1)

tl_df_long <- rbind(sum_demo_long, sum_tl_long)

tldf_sub <- tldf_sub %>% arrange(prev_tl)


library(GLMMadaptive)

negbinommod <-mixed_model(fixed = tl_val ~ survey * scientific_name, random = ~1 | location_name, 
                          data = tl_df_long, family = zi.negative.binomial(), 
                          zi_fixed = ~ scientific_name, zi_random = NULL,
                          max_coef_value = 30,
                          control = list(iter_EM = 0))

AIC(negbinommod)
hist(resid(negbinommod))
Anova(negbinommod)

binommod <- mixed_model(fixed = tl_val ~ survey * scientific_name, random = ~1 | location_name, 
                        data = tl_df_long, family = binomial(), 
                        max_coef_value = 30,
                        control = list(iter_EM = 0))

AIC(binommod)
hist(resid(binommod))
Anova(binommod)

zibinom <- mixed_model(fixed = tl_val ~ survey * scientific_name, random = ~1 | location_name, 
                       data = tl_df_long, family = zi.binomial(), 
                       zi_fixed = ~ scientific_name, zi_random = NULL,
                       max_coef_value = 30,
                       control = list(iter_EM = 0))

AIC(zibinom)

hist(resid(zibinom))

Anova(zibinom)

betabinom <- mixed_model(
  fixed = cbind(n_tl, n_healthy) ~ survey * scientific_name,
  random = ~ 1 | location_name,
  data = tldf_sub,
  family = beta.binomial(),
  max_coef_value = 30,
  control = list(iter_EM = 0))

AIC(betabinom)
hist(resid(betabinom))
summary(betabinom)

summary <- tldf %>% group_by(scientific_name, survey) %>%
  summarize(mean = mean(prev_tl))

library(emmeans)
library(rcompanion)


isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}

emm <- emmeans(betabinom, ~ survey*scientific_name)
simple <- pairs(emm, simple = "survey")
pairwise <- as.data.frame(pairs(emm, simple = "survey"))

eff <- as.data.frame(eff_size(emm, sigma = 26, edf = Inf))


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

###########
#most species for supplement (those that appear every year and/or aren't high suscep)
sitecolors = c('CBC Central'='red3','CBC Lagoon'='darkorchid','CBC30N'='gold1',
               'Curlew Patch' = 'blue', 'House Reef' ='seagreen4', 
               'South Reef Central' = 'orange', 'SR30N' = 'pink')


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
  scientific_name == "Siderastrea siderea"|
    scientific_name == "Pseudodiploria strigosa"|scientific_name == "Montastraea cavernosa"|
    scientific_name == "Orbicella spp.") %>%
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



#PrevLetters <- df %>%
#unite("GroupEvent", c(Species,Group))

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


