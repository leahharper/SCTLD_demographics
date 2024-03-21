setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/PAN-BDT Sept 2022")

library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(ggmap)
library(sp)
library(ggrepel)



corals <- read_csv("PAN-BDT_ColonyData.csv")


corals <- corals %>% mutate(Meters_90 = ifelse(Direction == "left" | Direction == "L",
                                                 -Meters_90, Meters_90))


speccolors = c('SSID'='red3','MCAV'='darkorchid1','PAST'='gold1',
               'MMEA' = 'lightblue', 'PSTR' ='seagreen4', 'CNAT'='darkorange2', 'ORBI'='deeppink4',
               'DSTO'='goldenrod4','EFAS'='pink','MLAM'='navy', 'DLAB' = 'limegreen')


##STRI
STRI <- corals %>% subset(Transect == "STRI Reef")

tiff("STRI_23.tif",width = 6, height = 10, units = "in", res = 300)
ggplot(STRI, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter)) +
  geom_point(pch = 21, color = "black") +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=STRI, aes(x=Meters_90, y=Meter, label=Current_tag_num), color="black", size = 4, hjust=-0.2, nudge_x = -0.5, max.overlaps = 20) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-16,3,by=1)) +
  scale_y_continuous("Meter", breaks = seq(-0,60,by=5)) +
  scale_size_continuous(range = c(3,6.5), name = "Diameter") +
  scale_fill_manual("Species",values=c(speccolors)) +
  ggtitle("STRI Reef") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()

##Juan Point
juan <- corals %>% subset(Transect == "Juan Point Reef")

tiff("JuanPoint_23.tif",width = 7, height = 10, units = "in", res = 100)
ggplot(juan, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter)) +
  geom_point(pch = 21, color = "black") +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=juan, aes(x=Meters_90, y=Meter, label=Current_tag_num), color="black", size = 4, hjust=-0.2, nudge_x = -0.7, max.overlaps = 20) +
  scale_y_continuous("Transect Length (m)", breaks = seq(0,50,by=5)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-6,12,by=1)) +
  scale_size_continuous(range = c(3,6.5), name = "Diameter") +
  scale_fill_manual("Species",values=c(speccolors)) +
  ggtitle("Juan Point") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()

##Crawl Cay
crawl <- corals %>% subset(Transect == "Crawl Cay")

tiff("CrawlCay_23.tif",width = 6, height = 7, units = "in", res = 100)
ggplot(crawl, aes(x = Meters_90, y = Meter, fill = Species, size = MaxDiameter)) +
  geom_point(pch = 21, color = "black") +
  geom_vline(xintercept = 0) +
  geom_text_repel(data=crawl, aes(x=Meters_90, y=Meter, label=Current_tag_num), max.overlaps = 20, color="black", size = 4, hjust=-0.2, nudge_x = -0.5) +
  scale_y_continuous("Transect Length (m)", breaks = seq(0,20,by=2)) +
  scale_x_continuous("Meters Perpendicular", breaks = seq(-6,4,by=1)) +
  scale_size_continuous(range = c(3,6.5), name = "Diameter") +
  scale_fill_manual("Species",values=c(speccolors)) +
  ggtitle("Crawl Cay") +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", hjust = 1, size = 12),
        axis.title = element_text(size = 12))
dev.off()