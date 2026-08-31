library(nlme)
library(plyr)
library(car)
library(tidyverse)
library(reshape2) 
library(piecewiseSEM)
library(ggh4x)

#################
#Fate Tracking
#################

# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/GitHub/SCTLD_demographics/")

fate <- read_csv("CBC_ColonyData.csv")

# Establish time levels
TimeLevels <- c("Oct2019","Nov2019","Dec2019",
                "Jan2020", "Feb2020", "Mar2020" ,"Apr2020","May2020", "Jun2020","Jul2020",
                "Aug2020" ,"Sep2020", "Oct2020","Nov2020","Dec2020",
                "Jan2021","Feb2021","Mar2021" ,"Apr2021","May2021", "Jun2021",
                "Jul2021","Aug2021" ,"Sep2021", "Oct2021","Nov2021","Dec2021",
                "Jan2022","Feb2022","Mar2022" ,"Apr2022","May2022", "Jun2022","Jul2022",
                "Aug2022" ,"Sep2022", "Oct2022","Nov2022","Dec2022")


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

Long$Time_Point <- recode(Long$Time_Point, "Condition_052022" = "May2022",
                          "Condition_122022" = "Dec2022",
                          "Condition_062019" = "Oct2019")


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


condcolors2 = c('Dead'='#581845',                  # cool dark red
                'Diseased'='#F0B400',               # yellow
                'Healthy'='#6B8E23',                 # warm medium green
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

sankeydf$x <- recode(sankeydf$x, "Condition_062019" = "Oct2019",
                     "Condition_052022" = "May2022",
                     "Condition_122022" = "Dec2022")

sankeydf$Species <- recode(sankeydf$Species, "MCAV" = "Montastraea cavernosa",
                           "SSID" = "Siderastrea siderea",
                           "PAST" = "Porites astreoides",
                           "MMEA" = "Meandrina meandrites",
                           "PSTR" = "Pseudodiploria strigosa")

sankeydf <- sankeydf %>% mutate(x = if_else(Species != "Pseudodiploria strigosa" & x == "Oct2019",
                                                   "Jun2019", x))


sankeydf$Species <- factor(sankeydf$Species,
                           levels = c("Siderastrea siderea", "Porites astreoides", 
                                      "Montastraea cavernosa", "Meandrina meandrites", "Pseudodiploria strigosa"))

sankeydf$node <- factor(sankeydf$node, levels = c("Dead", "Diseased", "Healthy"))

sankeydf <- sankeydf %>% mutate(x = fct_relevel(x,
                                        "Jun2019", "Oct2019", "May2022", "Dec2022")) 

common_ymax <- max(sankeydf$count) * 1.5

library(patchwork)

make_sankey <- function(df, labs = NULL, show_x = TRUE, ymax = NULL, total = NULL, show_legend = FALSE) {
  pad <- (ymax - total) / 2
  p <- ggplot(df, aes(x = x, next_x = next_x, node = node,
                      next_node = next_node, fill = factor(node))) +
    geom_sankey(flow.alpha = 0.6, node.color = 'black', flow.color = 'black') +
    geom_sankey_label(aes(label = after_stat(freq)), hjust = 1.1,
                      size = 7 / .pt, color = "black", fill = "white") +
    scale_fill_manual("Condition", values = condcolors2,
                      breaks = c("Dead", "Diseased", "Healthy"),
                      drop = FALSE) +
    coord_cartesian(ylim = c(-ymax / 2, ymax / 2))  +
    ggtitle(unique(df$Species)) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "italic"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = if (show_legend) "right" else "none")
  
  
  if (!is.null(labs)) p <- p + scale_x_discrete(labels = labs)
  if (show_x) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
  } else {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  p
}


extract_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[idx]]
}

ssid <- subset(sankeydf, Species == "Siderastrea siderea")
p_legend_source <- make_sankey(ssid,
                               show_x = FALSE, show_legend = TRUE, ymax = common_ymax)

shared_legend <- extract_legend(p_legend_source)

p_ssid <- make_sankey(ssid, show_x = FALSE, ymax = common_ymax)

past <- subset(sankeydf, Species == "Porites astreoides")
p_past <- make_sankey(past, show_x = TRUE, ymax = common_ymax)

mcav <- subset(sankeydf, Species == "Montastraea cavernosa")
p_mcav <- make_sankey(mcav, show_x = FALSE, ymax = common_ymax)

mmea <- subset(sankeydf, Species == "Meandrina meandrites")
p_mmea <- make_sankey(mmea, show_x = TRUE, ymax = common_ymax)

pstr <- subset(sankeydf, Species == "Pseudodiploria strigosa")
p_pstr <- make_sankey(pstr,
                      labs = c("Oct2019","May2022","Dec2022"), show_x = TRUE, ymax = common_ymax)


sankey <- (p_mcav | p_ssid | wrap_elements(full = shared_legend) ) /
  (p_past |p_mmea | p_pstr) +
  plot_layout(guides = "collect") 

sankey

png("Figures/current/Fig 3.png",width = 8, height = 5, units = "in", res = 300)
sankey
dev.off()
