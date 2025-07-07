library(openxlsx)
library(tidyverse)
library(scales)


setwd("G:/Meine Ablage/_Masterthesis/10_Biosynthethic gene clusters")

tbl_original <- read.xlsx("RNASeq_reshaped_output_full.xlsx")
glimpse(tbl_original)

tbl_filtered <- tbl_original %>% 
  filter(Cluster_type == "CDS") %>%
  select(Reads_total_crude.glycerol.pH4, Reads_total_crude.glycerol.pH7, Reads_total_pure.glycerol.pH7) %>%
  rename("Crude Glycerol, pH4" = Reads_total_crude.glycerol.pH4,
         "Crude Glycerol, pH7" = Reads_total_crude.glycerol.pH7, 
         "Pure Glycerol, pH7" = Reads_total_pure.glycerol.pH7)

tbl_filtered %>% summary()

tbl_filtered_pivot <- tbl_original %>% 
  filter(Cluster_type == "CDS") %>%
  select(Reads_total_crude.glycerol.pH4, Reads_total_crude.glycerol.pH7, Reads_total_pure.glycerol.pH7) %>%
  rename("Crude Glycerol, pH4" = Reads_total_crude.glycerol.pH4,
         "Crude Glycerol, pH7" = Reads_total_crude.glycerol.pH7, 
         "Pure Glycerol, pH7" = Reads_total_pure.glycerol.pH7) %>%
  pivot_longer(
    cols = c("Pure Glycerol, pH7", 
             "Crude Glycerol, pH7", 
             "Crude Glycerol, pH4"),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  mutate(Value = Value +1)

glimpse(tbl_filtered_pivot)

#Variante 2 -- bisschen chaotischer?
ggplot(tbl_filtered_pivot, aes(x = Variable, y = Value, fill = Variable)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, size = 0.7, alpha = 0.5) +
  scale_y_log10(
    breaks = 10^seq(0, 6, by = 1),
    labels = comma_format(scale = 1),
    guide  = "axis_logticks"
  ) +
  #annotation_logticks(sides = "l") + 
  scale_fill_manual(
    values = c(
      "Pure Glycerol, pH7" = "green",
      "Crude Glycerol, pH7" = "blue",
      "Crude Glycerol, pH4" = "blue"
    ) # end values
  ) +
  labs(
    title    = "Reads per coding sequence (CDS)",
    x        = "Samples",
    y        = "Reads (Log10)"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    axis.title    = element_text(size = 12)
  )


