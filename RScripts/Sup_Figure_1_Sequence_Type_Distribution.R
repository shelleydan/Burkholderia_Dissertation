# Dependancies
library(tidyverse)
library(dplyr)

metadata <- read.csv("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/metadata.csv",
                     header = TRUE,
                     row.names = 1)


# Check unique levels for ST
unique(metadata$Sequence_Type)

# Summarise counts per ST
ST_counts <- metadata %>%
  group_by(Sequence_Type) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# Plot STs in acending order. 
#svg("Cenocepacia_Sequence_Type.svg", height = 10, width = 10)

ggplot(ST_counts, aes(x = reorder(Sequence_Type, -Count), y = Count)) +
  
  # Customise bars.
  geom_bar(stat = "identity", fill = "cornflowerblue") +
  
  # Assigning bar labels
  geom_text(aes(label = paste0("n=", Count)), hjust = -0.1, size = 3.5, angle = 0) +
  
  # Assign x-axis label
  xlab("Sequence Type") +
  
  # Assign y-axis label
  ylab("Count (n)") +
  
  # Customise Axis Labels
  theme(axis.text.x = element_text(angle = 0, hjust = 0.4, size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip()

#dev.off()