setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/Little Cayman 2023")

library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(ggmap)
library(sp)
library(ggrepel)



corals <- read_csv("CCMI_ColonyData_LHplot.csv")


corals <- corals %>% mutate(Meters_90 = ifelse(Direction == "left" | Direction == "L",
                                                 -Meters_90, Meters_90))

corals$Species <- recode(corals$Species, "Efas" = "EFAS", "Dlab" = "DLAB", "Mlam" = "MLAM",
                         "Dsto" = "DSTO")

speccolors = c('SSID'='red3','MCAV'='darkorchid1','PAST'='gold1',
               'MMEA' = 'lightblue', 'PSTR' ='seagreen4', 'CNAT'='darkorange2', 'OFAV'='deeppink4',
               'DSTO'='goldenrod4','EFAS'='pink','MLAM'='navy', 'DLAB' = 'limegreen')


##Sailfin
sailfin <- corals %>% subset(TransectName == "Sailfin")

tiff("Sailfin_23.tif",width = 6, height = 8, units = "in", res = 300)
ggplot(sailfin, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter)) +
  geom_point(pch = 21, color = "black") +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=sailfin, aes(x=Meters_90, y=Meter, label=NewTag), color="black", size = 4, hjust=-0.2, nudge_x = -0.5, max.overlaps = 20) +
  #geom_text_repel(data=sailfin, aes(x=Meters_90, y=Meter, label=Species), color="black", size = 3, hjust=0, nudge_x = 1, max.overlaps = 20) +  scale_y_continuous("Transect Length (m)", breaks = seq(0,33,by=1)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-13,13,by=1)) +
  scale_size_continuous(range = c(3,6.5), name = "Diameter") +
  scale_fill_manual("Species",values=c(speccolors)) +
  #scale_fill_grey("Species") +
  ggtitle("Sailfin") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()

##Martha's Finyard
martha <- corals %>% subset(TransectName == "Marthas_Finyard")

tiff("Marthas_23.tif",width = 7, height = 9, units = "in", res = 100)
ggplot(martha, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter)) +
  geom_point(pch = 21, color = "black") +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=martha, aes(x=Meters_90, y=Meter, label=NewTag), color="black", size = 4, hjust=-0.2, nudge_x = -0.7, max.overlaps = 20) +
  #geom_text_repel(data=martha, aes(x=Meters_90, y=Meter, label=Species), color="black", size = 3, hjust=0, nudge_x = 1, max.overlaps = 20) +
  scale_y_continuous("Transect Length (m)", breaks = seq(0,35,by=1)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-12,5,by=1)) +
  scale_size_continuous(range = c(3,6.5), name = "Diameter") +
  scale_fill_manual("Species",values=c(speccolors)) +
  #scale_fill_grey("Species") +  
  ggtitle("Martha's Finyard") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()

##Jigsaw Puzzle
jig <- corals %>% subset(TransectName == "Jigsaw")

tiff("Jigsaw_23.tif",width = 8, height = 9, units = "in", res = 100)
ggplot(jig, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter)) +
  geom_point(pch = 21, color = "black") +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=jig, aes(x=Meters_90, y=Meter, label=NewTag), color="black", size = 4, hjust=-0.2, nudge_x = -0.5) +
  #geom_text_repel(data=jig, aes(x=Meters_90, y=Meter, label=Species), color="black", size = 3, hjust=0, nudge_x = 1) +
  scale_y_continuous("Transect Length (m)", breaks = seq(0,40,by=2)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-12,9,by=1)) +
  scale_size_continuous(range = c(3,6.5), name = "Diameter") +
  scale_fill_manual("Species",values=c(speccolors)) +
  #scale_fill_grey("Species") +
  ggtitle("Jigsaw Puzzle") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()