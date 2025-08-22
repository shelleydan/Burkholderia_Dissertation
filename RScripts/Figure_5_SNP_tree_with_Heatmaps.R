library(ggtree)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(tidyr)
library(ggrepel)
library(Polychrome)
library(ape)
library(tidyverse)

# Read in Metadata
metadata <- read.csv("C:/Users/Danie/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/metadata.csv",
                     header = TRUE)
#row.names = 1)
metadata$Clinical_Environmental
colnames(metadata)[1] <- "label"

# Read in SNP RAxML support tree
tree_data <- read.tree("C:/Users/Danie/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/data/snp_tree.raxml.support")

tree_data <- drop.tip(tree_data, "Reference") # Remove Ref to match metadata

## COLORUS #####
source_colours <- c(
  "Clinical_(CF)" = "cornflowerblue",
  "Clinical_(Non-CF)" = "lightsteelblue",
  "Clinical_(Unknown)" = "lightskyblue4",
  "Environmental" = "seagreen2",
  "Laboratory_Experiment" = "ivory2",
  "Manufacturer" = "brown2",
  "Not_Provided" = "coral4"
)

STs <- c("1813", "241", "Unknown_(No_DB_Match)", "32", "1387", "208",
         "248", "621", "306", "28", "227", "709", "210", "234", "1922",
         "1952", "1921", "1943", "1929", "1495", "665", "1924", "1632",
         "33", "278", "1945", "242", "1506", "1502", "258", "1499", "216",
         "824", "841", "250", "209", "807", "309", "201", "1869", "218",
         "358", "230", "1628", "2147", "1076", "964", "728", "217", "674",
         "2128", "1414", "1475", "1477", "839", "624", "1798", "1299", "1300", NA)

GroupColours %in% unique(tree1$data$Sequence_Type)

set.seed(0118999)

GroupColours <- setNames(
  createPalette(
    length(STs),
    seedcolors = c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"),
    range = c(30, 90)  # Limit luminance range to avoid dark/light overlap
  ),
  STs
)

#### PHYLOGENY WITH ANNOTATIONS
# data_matrix1 = AMR
# data_matrix2 = VF
# data_matrix3 = BCESM
# data_matrix1 = Source Information


data_matrix1 <- read.csv("C:/Users/Danie/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/data/customdb_final_.csv", 
                         header = TRUE, row.names = 1)
data_matrix2 <- read.csv("C:/Users/Danie/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/data/vf_final_.csv", 
                         header = TRUE, row.names = 1)

data_matrix3 <- read.csv("C:/Users/Danie/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/data/BCESM_final.csv", 
                         header = TRUE, row.names = 1)

data_matrix4 <- read.csv("C:/Users/Danie/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/data/metadata_clinical.csv", 
                         header = TRUE)

#### Data formatting
# Replace . with 0s
# Extract only first % value (if identified >1)
# Ensure numeric value

data_matrix1[data_matrix1 == "."] <- 0
data_matrix1[] <- lapply(data_matrix1, function(x) sub(";.*", "", x))
data_matrix1[] <- lapply(data_matrix1, function(x) as.numeric(as.character(x)))


data_matrix2[data_matrix2 == "."] <- 0
data_matrix2 <- data_matrix2 %>% select(-NUM_FOUND)
data_matrix2[] <- lapply(data_matrix2, function(x) sub(";.*", "", x))
data_matrix2[] <- lapply(data_matrix2, function(x) as.numeric(as.character(x)))

data_matrix3[data_matrix3 == "."] <- 0
data_matrix3[] <- lapply(data_matrix3, function(x) sub(";.*", "", x))
data_matrix3[] <- lapply(data_matrix3, function(x) as.numeric(as.character(x)))

data_matrix4$Model <- as.factor(data_matrix4$Model)

### Generate Tree plot
tree1 <- ggtree(tree_data, 
                layout = "rectangular",
                branch.length = 'none') %<+% metadata 

# Add Branch Colours for STs
tree1$data$branch_color <- GroupColours[tree1$data$Sequence_Type]

tree2 <- tree1 +
  geom_tree(aes(color = I(branch_color)), size = 1) +
  geom_tippoint(aes(color = Sequence_Type), size = NA) +
  scale_color_manual(
    name = "Sequence Type",
    values = GroupColours
  ) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing = unit(0.75, "cm")
  )  +
  guides(color = guide_legend(
    override.aes = list(
      shape = 15,   # solid square
      size = 5,     # square size
      linetype = 0  # no line in legend
    )
  )) +
  theme(
    panel.background = element_rect(fill =  "transparent"), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

### Add sequential heatmaps

tree3 <- gheatmap(tree2, data_matrix1, offset = 1, width = 1, 
                  font.size = 7.5, hjust = 0, colnames_angle = 90, colnames = TRUE,
                  colnames_offset_y = -30) +
  scale_fill_gradient(low = "white", high = "blue", name = "AMR") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing = unit(0.75, "cm")
  ) 

tree3 <- tree3 + ggnewscale::new_scale_fill()

tree4 <- gheatmap(tree3, data_matrix2, offset = 37, width = 1, 
                  font.size = 7.5, hjust = 0, colnames_angle = 90, colnames = TRUE,
                  colnames_offset_y = -30) +
  scale_fill_gradient(low = "white", high = "red", name = "VF")+
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing = unit(0.75, "cm")
  ) 

tree4 <- tree4 + ggnewscale::new_scale_fill()

tree5 <- gheatmap(tree4, data_matrix3, offset = 70, width = 0.05, 
                  font.size = 7.5, hjust = 0, colnames_angle = 90, colnames = TRUE,
                  colnames_offset_y = -35) +
  scale_fill_gradient(low = "white", high = "limegreen", name = "BCESM") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30),
    legend.key.size = unit(1, "cm"),
    legend.spacing = unit(1, "cm")
  ) 

tree5 <- tree5 + ggnewscale::new_scale_fill()

final_plot <- tree5 +
  geom_fruit(
    data = data_matrix4,
    geom = geom_tile,
    mapping = aes(y = Run, fill = Model),
    width = 2,
    offset = 2.3
  ) +
  scale_fill_manual(values = source_colours, name = "Source")

final_plot <- final_plot +
  theme(
    plot.margin = unit(c(2, 2, 2, 2), "cm")  # top, right, bottom, left
  )

## Save the final plot.
ggsave("SNP_ANNOTATION_RECTANGLE.svg", 
       limitsize = FALSE,
       height = 35, 
       width = 70)

