library(ComplexHeatmap)

####### COMPLEXHEATMAP

# Read in the ANI Matrix
data <- read.table("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/data/full_identity_matrix.tsv",
                   header = TRUE,
                   row.names = 1)
str(data)

# Read in the metadata with confirmed ANI Spp in column 1. 
metadata <- read.csv("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/Burkholderia_All_meta.csv",
                     header = TRUE,
                     row.names = 1)

str(metadata)

# Ensure IDs match
filtered_metadata <- metadata[rownames(metadata) %in% rownames(data), ]

# Remove all columns apart from confirmed Spp. 
filtered_metadata <- filtered_metadata[,1, drop=FALSE]

# Assign colours for 95-100%
colours <- circlize::colorRamp2(c(0,94.99,95,97.5,100), 
                                c("white","white","cornflowerblue","yellow","red"))

# Converting original Dataframe into a matrix.
matrix_data <- as.matrix(data)

# Ensure both datasets have the same labels. 
filtered_metadata <- filtered_metadata[rownames(matrix_data), , drop = FALSE]

#Confirm all labels match
identical(rownames(matrix_data), rownames(filtered_metadata))

# Get unique species to assign colours too. 
unique(filtered_metadata$ANI.Confirmed.Spp)

species_colours <- c(
  "Burkholderia cenocepacia" = "#6495ED",  
  "Burkholderia cepacia" = "pink",   
  "Burkholderia contaminans" = "#FFFF00",    
  "Burkholderia orbicola" = "#FF0000",    
  "Burkholderia pseudomultivorans" = "#008000",     
  "Burkholderia semiarida" = "#FFA500",    
  "Burkholderia sola" = "purple",  
  "Non-Bcc" = 'gray'
)

# Generate the annotation based on species colours. 
col_anno <- columnAnnotation(
  df = filtered_metadata,
  col = list(ANI.Confirmed.Spp = species_colours),
  annotation_legend_param = list(
    ANI.Confirmed.Spp = list(
      title = "Species",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    )
  )
)

# Build the species legend. 
species_legend <- Legend(
  at = names(species_colours),              # categories
  legend_gp = gpar(fill = species_colours), # colors for categories
  title = "Species",
  title_gp = gpar(fontsize = 16, fontface = "bold"),
  labels_gp = gpar(fontsize = 14)
)

# Create the heatmap. 
heatmap_plot <- Heatmap(
  matrix_data,
  
  #Allow ANI clustering
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  
  
  col = colours,
  
  # Remove IDs to avoid clutter
  show_column_names = FALSE,
  show_row_names = FALSE,
  
  # Set Heatmap dimensions
  width = unit(50, "cm"),   
  height = unit(50, "cm"),
  
  # Add in Spp annotations
  top_annotation = col_anno,
  
  # Assign tree height (for row and column)
  column_dend_height = unit(50, "mm"),
  row_dend_width = unit(50, "mm"),
  
  # Add in customised colouring for ANI %
  heatmap_legend_param = list(
    title = "ANI (%)",
    at = c(95, 97.5, 100),        
    labels = c("95", "97.5", "100"),
    legend_height = unit(10, "cm"),
    title_position = "topcenter",   
    color_bar = "continuous",
    labels_gp = gpar(fontsize = 14),
    title_gp = gpar(fontsize = 16, fontface = "bold")
  )
)

# Save the Heatmap
svg("Burk25-Full-Heatmap.svg", height = 30, width = 30)

# Draw heatmap and legends together
draw(heatmap_plot, heatmap_legend_list = list(species_legend))

dev.off()




