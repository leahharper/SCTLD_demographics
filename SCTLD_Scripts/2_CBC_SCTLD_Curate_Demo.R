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

demo <- demo %>% subset(year < 2023)

dlab <- demo %>% subset(scientific_name == "Diploria labyrinthiformis"|
                          scientific_name == "Pseudodiploria strigosa")  %>%
  subset(location_name == "South Reef 30 North")

#correct Vivi's probably misidentified DLABs at SR30N
demo <- demo %>%
    mutate(scientific_name = ifelse(diver == "Viviana Bravo" &
                                 scientific_name == "Diploria labyrinthiformis" &
                               location_name == "South Reef 30 North",
                         "Pseudodiploria strigosa", scientific_name)) 


#3 EFAS and 1 PSTR at SR Central were totally dead!
#PSTR = 6cm
#EFAS at 6, 10, and 13cm

# #this code removes dead colonies from demo
# demo <- demo %>%
#  mutate(X5_to_10_cm = ifelse(diver == "Leah Harper" &
#                               month == 5 &
#                              year == 2022 & 
#                             location_name == "South Reef Central" &
#                          scientific_name == "EFAS", 
#                        X5_to_10_cm - 2, X5_to_10_cm)) %>%
#  mutate(X11_to_20_cm = ifelse(diver == "Leah Harper" &
#                               month == 5 &
#                              year == 2022 & 
#                             location_name == "South Reef Central" &
#                            scientific_name == "EFAS", 
#                         X11_to_20_cm - 1, X11_to_20_cm)) %>%
# mutate(X5_to_10_cm = ifelse(diver == "Leah Harper" &
#                              month == 5 &
#                             year == 2022 & 
#                            location_name == "South Reef Central" &
#                           scientific_name == "PSTR", 
#                        X5_to_10_cm - 1, X5_to_10_cm))


demo <- demo %>%
  mutate(total_small = rowSums(.[17:20])) %>% #5-40cm
  mutate(total_lg = rowSums(.[21:22])) %>% #>40cm
  mutate(adult = rowSums(.[17:22])) %>%
  rename("total" = Total)#>5cm

demo$scientific_name <- recode(demo$scientific_name, "Porites porities" = "Porites porites")
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

#tiff("Figures/adultNMDS_jaccard.tif",width = 7, height = 7, units = "in", res = 400)
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
#dev.off()

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

tiff("Figures/current/Fig 5.png",width = 7, height = 7, units = "in", res = 400)
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
rownames(braydf) <- c("Transect", "Survey Timepoint", "Residual", "Total")


braytable <- braydf %>%
  kbl(caption = "<span style='color: black;'> <b>Table 4.</b> PERMANOVA Results. <span>",
          digits = 4) %>%
  kable_styling()

save_kable(braytable, file = "Tables/Table 4.pdf")




