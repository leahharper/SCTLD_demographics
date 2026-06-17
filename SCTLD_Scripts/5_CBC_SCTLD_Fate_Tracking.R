library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)

#################
#Fate Tracking
#################

# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/GitHub/SCTLD_demographics/")

fate <- read_csv("CBC_ColonyData.csv")

# Establish time levels
TimeLevels <- c("October19","November19","December19",
                "January20", "February20", "March20" ,"April20","May20", "June20","July20",
                "August20" ,"September20", "October20","November20","December20",
                "January21","February21","March21" ,"April21","May21", "June21",
                "July21","August21" ,"September21", "October21","November21","December21",
                "January22","February22","March22" ,"April22","May22", "June22","July22",
                "August22" ,"September22", "October22","November22","December22")


# correct coding errors ----
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
  subset(Species != "CNAT") %>%
  mutate(Condition_052022 = if_else(Condition_052022 == "CLB,DC",
                                    "Healthy", Condition_052022)) %>%
  mutate(Condition_122022 = if_else(Condition_122022 == "CLB,DC"|
                                      Condition_122022 == "CLB",
                                    "Healthy", Condition_122022))

# identify recoveries----
rec <- fate4 %>% subset(Condition_052022 == "Diseased" &
                          Condition_122022 == "Healthy")


Long <-melt(fate4, id.vars =c("Species","NewTagNum",
                              "Transect"), measure.vars = 
              c("Condition_052022", "Condition_122022", "Condition_062019")) %>%
  rename("Time_Point" = variable, "Condition" = value) 

Long$Condition <- recode(Long$Condition, "Healthy?" = "Healthy")

Long$Time_Point <- recode(Long$Time_Point, "Condition_052022" = "May22",
                          "Condition_122022" = "December22",
                          "Condition_062019" = "October19")


SpecLevels = c('SSID','MCAV','PSTR','PAST','MMEA','CNAT')

Long <- Long %>% mutate_at(.vars = vars("Time_Point"), 
                           .funs = funs(factor(.,levels = TimeLevels, ordered = TRUE))) %>%
  mutate_at(.vars = vars("Species"),
            .funs = funs(factor(.,levels = SpecLevels, ordered = TRUE)))

# summaries----

sumsimple <- Long %>% group_by(Time_Point,Species) %>% count(Condition) %>% ungroup() %>%
  unite("timespp", c("Time_Point", "Species"), remove = FALSE) %>%
  rename("n_cond" = n)

sumsimpler <- Long %>% group_by(Time_Point) %>% count(Condition) 

sumtime <- Long %>% group_by(Time_Point, Species) %>% count() %>%
  unite("timespp", c("Time_Point", "Species"), remove = TRUE)

sum <- sumsimple %>% left_join(sumtime, by = "timespp") %>%
  mutate(percent = (n_cond/n)*100)


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

#remotes::install_github("davidsjoberg/ggsankey")
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

sankeydf$x <- recode(sankeydf$x, "Condition_062019" = "October19",
                     "Condition_052022" = "May22",
                     "Condition_122022" = "December22")

sankeydf$Species <- recode(sankeydf$Species, "MCAV" = "Montastraea cavernosa",
                           "PSTR" = "Pseudodiploria strigosa",
                           "SSID" = "Siderastrea siderea",
                           "PAST" = "Porites astreoides",
                           "MMEA" = "Meandrina meandrites")

sankeydf <- sankeydf %>% mutate(Species = fct_relevel(Species,
                                                      "Siderastrea siderea", "Porites astreoides", 
                                    "Meandrina meandrites", "Montastraea cavernosa", "Pseudodiploria strigosa")) 

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
        strip.text = element_text(size = 12, face = "italic"),
        axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 12),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.position = c(0.8,0.8),
        legend.title = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"),
        legend.text = element_text(size = 12))
sankey

png("Figures/current/Fig 3.png",width = 9, height = 6, units = "in", res = 300)
sankey
dev.off()
