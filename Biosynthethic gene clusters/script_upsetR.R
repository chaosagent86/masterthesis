library(openxlsx)
library(tidyverse)
library(scales)
library(UpSetR)


setwd("G:/Meine Ablage/_Masterthesis/10_Biosynthethic gene clusters")

tbl_original <- read.xlsx("RNASeq_reshaped_output_full.xlsx")
glimpse(tbl_original)

tbl_filtered <- tbl_original %>% 
  filter(Cluster_type == "CDS") %>%
  select("Cluster-Name",Reads_total_crude.glycerol.pH4, Reads_total_crude.glycerol.pH7, Reads_total_pure.glycerol.pH7) %>%
  rename(Crude_Glycerol_pH4 = Reads_total_crude.glycerol.pH4,
         Crude_Glycerol_pH7 = Reads_total_crude.glycerol.pH7, 
         Pure_Glycerol_pH7 = Reads_total_pure.glycerol.pH7,
         CDS = "Cluster-Name") %>%
  group_by(CDS) %>%
    summarise(
      Crude_Glycerol_pH4 = mean(Crude_Glycerol_pH4, na.rm = TRUE),
      Crude_Glycerol_pH7 = mean(Crude_Glycerol_pH7, na.rm = TRUE),
      Pure_Glycerol_pH7 = mean(Pure_Glycerol_pH7, na.rm = TRUE)
    ) %>%
  mutate(
    Crude_Glycerol_pH4  = as.integer(Crude_Glycerol_pH4  > 0),
    Crude_Glycerol_pH7  = as.integer(Crude_Glycerol_pH7  > 0),
    Pure_Glycerol_pH7   = as.integer(Pure_Glycerol_pH7   > 0)
  ) %>%
  column_to_rownames("CDS")

upset(
  tbl_filtered,
  sets            = c("Crude_Glycerol_pH4", "Crude_Glycerol_pH7", "Pure_Glycerol_pH7"),
  nsets           = 3,            # alle drei Sets anzeigen
  nintersects     = NA,           # alle Intersektionen
  order.by        = "freq",       # nach Häufigkeit sortieren
  empty.intersections = "off"     # leere Schnittmengen ausblenden
)
