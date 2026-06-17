setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/BEL-CBC Dec 2022")

library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(ggmap)
library(sp)
library(ggrepel)


corals <- read_csv("CBC_ColonyData.csv")


corals$Direction <- dplyr::recode(corals$Direction, "L" = "left",
                                  "R" = "right")

corals <- corals %>% mutate(Meters_90 = ifelse(Direction == "left",
                                     -Meters_90, Meters_90)) %>%
                     mutate(Condition = ifelse(Date_DocumentedMortality != "Healthy"& Date_DocumentedMortality != "Diseased",
                              "Dead", Date_DocumentedMortality)) %>%
                    mutate(MaxDiameter = ifelse(is.na(MaxDiameter),
                            40, MaxDiameter))

speccolors = c('SSID'='red3','MCAV'='darkorchid','PAST'='gold1',
                            'MMEA' = 'lightblue', 'PSTR' ='seagreen4', 
                            'CNAT' = 'brown', 'OFAV' = 'pink', 'OANN' = 'blue',
                            'DLAB' = 'yellowgreen', "Mosaic" = 'black')

specalpha = c('Dead'= 1,'Diseased'= 0,'Healthy'= 0)

corals$MaxDiameter <- as.numeric(corals$MaxDiameter)

corals$Condition <- as.factor(corals$Condition)

corals$Species <- dplyr::recode(corals$Species, "OANN/OFAV?" = "OANN")


##CBC30N
CBC30N <- corals %>% subset(Transect == "CBC30N") %>%
  subset(!(is.na(NewTagNum)))


tiff("CBC30NMapTags.tif",width = 5, height = 8, units = "in", res = 300)
ggplot() +
  geom_point(data=CBC30N, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter),
             pch = 21, color = "black") +
  geom_point(data=CBC30N, aes(x = Meters_90, y = Meter, alpha = Condition),
             pch = 4, color = "black", stroke = 2) +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=CBC30N, aes(x=Meters_90, y=Meter, label=NewTagNum), color="black", size = 4, hjust=-0.2) +
  scale_y_continuous("Transect Length (m)", breaks = seq(0, 31, by = 1)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-5, 7, by = 1)) +
  scale_size_continuous(range = c(2,6.5), name = "", guide = 'none') +
  scale_fill_manual("Species",values=c(speccolors)) +
  scale_alpha_manual("",values=c(specalpha), guide = 'none') +
  ggtitle("CBC 30 N") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()

##CBC Lagoon
Lagoon <- corals %>% subset(Transect == "Lagoon") %>%
  subset(NewTagNum != "not_found")

tiff("CBCLagoonMapTags.tif",width = 5, height = 8, units = "in", res = 300)
ggplot() +
  geom_point(data=Lagoon, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter),
             pch = 21, color = "black") +
  geom_point(data=Lagoon, aes(x = Meters_90, y = Meter, alpha = Condition),
             pch = 4, color = "black", stroke = 2) +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=Lagoon, aes(x=Meters_90, y=Meter, label=NewTagNum), max.overlaps = 20, color="black", size = 4, hjust=-0.2) +
  scale_y_continuous("Transect Length (m)", breaks = seq(0, 42, by = 1)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-12, 7, by = 1)) +
  scale_size_continuous(range = c(2,6.5), name = "", guide = 'none') +
  scale_fill_manual("Species",values=c(speccolors)) +
  scale_alpha_manual("",values=c(specalpha), guide = 'none') +
  ggtitle("CBC Lagoon") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()



#SR30N
SR30N <- corals %>% subset(Transect == "SR30N") %>%
  subset(NewTagNum != "not_found")

tiff("SR30NMapTags.tif",width = 5, height = 8, units = "in", res = 300)
ggplot() +
  geom_point(data=SR30N, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter),
             pch = 21, color = "black") +
  geom_point(data=SR30N, aes(x = Meters_90, y = Meter, alpha = Condition),
             pch = 4, color = "black", stroke = 2) +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=SR30N, aes(x=Meters_90, y=Meter, label=NewTagNum), max.overlaps = 40, color="black", size = 4, hjust=-0.2) +
  scale_y_continuous("Transect Length (m)", breaks = seq(0, 42, by = 1)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-12, 7, by = 1)) +
  scale_size_continuous(range = c(2,6.5), name = "", guide = 'none') +
  scale_fill_manual("Species",values=c(speccolors)) +
  scale_alpha_manual("",values=c(specalpha), guide = 'none') +
  ggtitle("SR30N") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()

#Curlew
Curlew <- corals %>% subset(Transect == "Curlew")

tiff("CurlewTags.tif",width = 5, height = 8, units = "in", res = 300)
ggplot() +
  geom_point(data=Curlew, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter),
             pch = 21, color = "black") +
  geom_point(data=Curlew, aes(x = Meters_90, y = Meter, alpha = Condition),
             pch = 4, color = "black", stroke = 2) +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=Curlew, aes(x=Meters_90, y=Meter, label=NewTagNum), max.overlaps = 40, color="black", size = 4, hjust=-0.2) +
  scale_y_continuous("Transect Length (m)", breaks = seq(0, 42, by = 1)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-6, 10, by = 1)) +
  scale_size_continuous(range = c(2,6.5), name = "", guide = 'none') +
  scale_fill_manual("Species",values=c(speccolors)) +
  scale_alpha_manual("",values=c(specalpha), guide = 'none') +
  ggtitle("Curlew") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()


