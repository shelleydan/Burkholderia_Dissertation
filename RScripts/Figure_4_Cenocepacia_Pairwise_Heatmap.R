library(ComplexHeatmap)
library(grid)
library(Polychrome)

# Read in ANI Matrix for B. cenocepacia
data <- read.table("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/identity_matrix.tsv",
                   header = TRUE,
                   row.names = 1)
str(data)

matrix_data <- as.matrix(data)       # Generate the matrix.

# Read in the Metadata
metadata <- read.csv("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/metadata.csv",
                     header = TRUE,
                     row.names = 1)

str(metadata)

# Keep only metadata rows that match rownames in the matrix
filtered_metadata <- metadata[rownames(metadata) %in% rownames(data), ]

# Extract the ST and Source columns
filtered_metadata <- filtered_metadata[,c(8,13), drop=FALSE]

data[is.na(data)] <- 0   # Replace any NA values with 0. 

# Assign Source for better naming
names(filtered_metadata)[names(filtered_metadata) == "Clinical_Environmental"] <- "Source"

# Ensure the length of the metadata and matrix match
all(length(filtered_metadata$Source) == ncol(matrix_data))

# Compare and match names and check if all are present. 
filtered_metadata <- filtered_metadata[rownames(matrix_data), , drop = FALSE]
identical(rownames(matrix_data), rownames(filtered_metadata))  # should return TRUE

## COLOURS ##########

#Assign source colours
source_colours <- c(
  "Clinical_(CF)" = "cornflowerblue",
  "Clinical_(Non-CF)" = "lightsteelblue",
  "Clinical_(Unknown)" = "lightskyblue4",
  "Environmental" = "seagreen2",
  "Laboratory_Experiment" = "ivory2",
  "Manufacturer" = "brown2",
  "Not_Provided" = "coral4"
)

#Generate a unique ST list.
STs <- c("1813", "241", "Unknown_(No_DB_Match)", "32", "1387", "208",
         "248", "621", "306", "28", "227", "709", "210", "234", "1922",
         "1952", "1921", "1943", "1929", "1495", "665", "1924", "1632",
         "33", "278", "1945", "242", "1506", "1502", "258", "1499", "216",
         "824", "841", "250", "209", "807", "309", "201", "1869", "218",
         "358", "230", "1628", "2147", "1076", "964", "728", "217", "674",
         "2128", "1414", "1475", "1477", "839", "624", "1798", "1299", "1300")

set.seed(0118999)

# Assign colours to STs
GroupColours <- setNames(
  createPalette(
    length(STs),
    seedcolors = c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"),
    range = c(30, 90)  # Limit luminance range to avoid dark/light overlap
  ),
  STs
)

# Generate ANI colours
colours <- circlize::colorRamp2(c(0,94.99,95,97.5,100), c("white","white","cornflowerblue","yellow","red"))
## ANNOTATION ##########
col_anno <- columnAnnotation(
  Source = filtered_metadata$Source,
  Sequence_Type = filtered_metadata$Sequence_Type,
  gp = gpar(fontsize = 20),
  col = list(
    Source = source_colours,
    Sequence_Type = GroupColours
  ),
  show_annotation_name = TRUE,
  show_legend = FALSE
  
)

## HEATMAP LEGEND ##########
ani_legend <- Legend(
  title = "ANI (%)",
  at = c(95, 96, 97, 98, 99, 100),
  labels = c("95", "96", "97", "98", "99", "100"),
  col_fun = colours, 
  legend_height = unit(10, "cm"),
  #title_position = "topcenter",
  labels_gp = gpar(fontsize = 18),
  title_gp = gpar(fontsize = 18, fontface = "bold")
)

## SOURCE LEGEND
labels <- names(source_colours)
fills <- as.vector(source_colours)  # force a plain character vector

# Define one fill per label
source_legend <- Legend(
  labels = labels,
  legend_gp = gpar(fill = fills),
  title = "Source",
  title_gp = gpar(fontsize = 18, fontface = "bold"),
  labels_gp = gpar(fontsize = 18)
)

## ST Legend
labels_ST <- names(GroupColours)
fills_ST <- as.vector(GroupColours)  # force a plain character vector

# Define one fill per label
ST_legend <- Legend(
  labels = labels_ST,
  legend_gp = gpar(fill = fills_ST),
  title = "Secquence Type",
  title_gp = gpar(fontsize = 18, fontface = "bold"),
  labels_gp = gpar(fontsize = 18)
)


pd = packLegend(ani_legend, source_legend, ST_legend)

HMPLOT <- Heatmap(
  matrix_data,
  top_annotation = col_anno,
  col = colours,
  width = unit(50, "cm"),
  height = unit(50, "cm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  
  #DEND
  column_dend_height = unit(50, "mm"),
  row_dend_width = unit(50, "mm"),
  show_heatmap_legend = FALSE
)


# Open SVG and save figure.
pdf("Burk25-CENOCEPACIA-Heatmap-TEST.pdf", height = 30, width = 30)

# Draw heatmap WITHOUT default legends
draw(
  HMPLOT,
  heatmap_legend_list = NULL,         # disables default ANI legend
  annotation_legend_list = NULL,      # disables annotation legends
  merge_legend = TRUE
)

# Manually draw your packed legends (pd) to the right
draw(
  pd,
  x = unit(1, "npc") - unit(8, "cm"),  # position 5 cm from right edge
  y = unit(0.47, "npc"),                # center vertically
  just = c("left", "center")
)

dev.off()
