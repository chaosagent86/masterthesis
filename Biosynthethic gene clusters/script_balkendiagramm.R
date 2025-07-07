library(openxlsx)
library(tidyverse)
library(stringr)

setwd("G:/Meine Ablage/_Masterthesis/10_Biosynthethic gene clusters")

tbl_original <- read.xlsx("RNASeq_reshaped_output_full.xlsx")
glimpse(tbl_original)

tbl_filtered <- tbl_original %>% 
  filter(Cluster_type == "cand_cluster") %>%
  select(Completeness_crude.glycerol.pH4, Completeness_crude.glycerol.pH7, Completeness_pure.glycerol.pH7) %>%
  rename("Crude Glycerol, pH4" = Completeness_crude.glycerol.pH4,
         "Crude Glycerol, pH7" = Completeness_crude.glycerol.pH7, 
         "Pure Glycerol, pH7" = Completeness_pure.glycerol.pH7) %>%
  pivot_longer(
    cols      = c("Crude Glycerol, pH4", 
                  "Crude Glycerol, pH7", 
                  "Pure Glycerol, pH7"),
    names_to  = "Sample",
    values_to = "Values"
  ) %>%
  mutate(
    Values = readr::parse_number(Values),
    Status = if_else(Values > 95, "complete", "incomplete")
  )
  
counts <- tbl_filtered %>% count(Sample, Status)

ggplot(counts, aes(Sample, n, fill = Status)) +
  geom_col(position = position_dodge(0.6), width = 0.5) +
  scale_fill_manual(values = c(complete = "green", incomplete = "red")) +
  labs(
    title    = "Complete vs. Incomplete Cluster",
    subtitle = "Threshold: Completeness > 95 %",
    x        = "Probe",
    y        = "Amount of Cluster",
    fill     = "Status"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    axis.title    = element_text(size = 12)
  )
