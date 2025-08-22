# Dependancies
library(tidyverse)
library(dplyr)

##### Data input
metadata <- read.csv("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/Burkholderia_All_meta.csv",
                     header = TRUE,
                     row.names = 1)

##### Data Visualizations

# Stacked Spp. w/ Source
metadata <- metadata %>%
  mutate(ANI.Confirmed.Spp = fct_infreq(ANI.Confirmed.Spp))

# Summarize totals per species
totals <- metadata %>%
  count(ANI.Confirmed.Spp)                               # Count Spp. 

# Generate counts per source for each Spp.
label_data <- metadata %>%
  count(ANI.Confirmed.Spp, Clinical_Environmental) %>%   # Count source by Spp.
  group_by(ANI.Confirmed.Spp) %>%                        # Group by Spp. 
  mutate(
    percent = n / sum(n) * 100,                          # Generate % Values
    label = paste0(round(percent,1), "% (", n, ")")      # Generate figure labels
  )

#svg("Burk25_Source_Breakdown.svg", height = 10, width = 10)

ggplot(metadata, 
       aes(fill=Clinical_Environmental,                  # Colour by Source
           x=ANI.Confirmed.Spp)) +                       # Bars by Spp. 
  geom_bar(position="stack") +                           # Generate stacked bars. 
  
  # Assign labels
  geom_text(
    data = label_data,
    aes(
      x = ANI.Confirmed.Spp,
      y = n,
      label = label,
      group = Clinical_Environmental
    ),
    position = position_stack(vjust = 0.5),
    size = 2.5,
    color = "black"
  ) +
  
  # Customising plot labels. 
  labs(
    x = "Species",
    y = "Count",
    fill = "Source"
  ) +
  
  # Customising source fill 
  scale_fill_manual(
    values = c(
      "Clinical_(CF)" = "cornflowerblue",
      "Clinical_(Non-CF)" = "lightsteelblue",
      "Clinical_(Unknown)" = "lightskyblue4",
      "Environmental" = "seagreen2",
      "Laboratory_Experiment" = "ivory2",
      "Manufacturer" = "brown2",
      "Not_Provided" = "coral4")
  ) +
  
  # Customise label positions. 
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

#dev.off()