library(ggtree)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(tidyr)
library(ggrepel)
library(Polychrome)
library(ape)
library(tidyverse)
library(phytools)
library(ggpubr)

# Read in Metadata
metadata <- read.csv("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP1_Genome_Collection/metadata.csv",
                     header = TRUE)

metadata$Clinical_Environmental
colnames(metadata)[1] <- "label"

# Read in the Pangenome tree from RAxML
core_tree <- read.tree("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/Data/pangenome_trees.raxml.support")

# Read in the SNP tree from RAxML
snp_tree <- read.tree("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/Data/snp_tree.raxml.support")
snp_tree <- drop.tip(snp_tree, c("Reference", "SRR25301686"))

# Check and match labels between trees. 
setdiff(core_tree$tip.label, snp_tree$tip.label)
setdiff(snp_tree$tip.label, core_tree$tip.label)

# Generate the distance matrix for each tree.
# This was used for PhyloTree and may be used again with more time. 
dist.topo(core_tree, snp_tree)

###### ST COlours
STs <- c("1813", "241", "Unknown_(No_DB_Match)", "32", "1387", "208",
         "248", "621", "306", "28", "227", "709", "210", "234", "1922",
         "1952", "1921", "1943", "1929", "1495", "665", "1924", "1632",
         "33", "278", "1945", "242", "1506", "1502", "258", "1499", "216",
         "824", "841", "250", "209", "807", "309", "201", "1869", "218",
         "358", "230", "1628", "2147", "1076", "964", "728", "217", "674",
         "2128", "1414", "1475", "1477", "839", "624", "1798", "1299", "1300", NA)

set.seed(0118999)

GroupColours <- setNames(
  createPalette(
    length(STs),
    seedcolors = c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"),
    range = c(30, 90)  # Limit luminance range to avoid dark/light overlap
  ),
  STs
)

tip_colors <- GroupColours[metadata$Sequence_Type]

tree1 <- ggtree(core_tree, branch.length = "none") %<+% metadata 

tree1$data$branch_color <- GroupColours[tree1$data$Sequence_Type]

tree2 <- tree1 +
  geom_tree(aes(color = I(branch_color)), size = 1) +
  geom_tippoint(aes(color = Sequence_Type), size = NA) +
  geom_tiplab(aes(label=""), align=TRUE) +
  scale_color_manual(
    name = "Sequence Type",
    values = GroupColours
  ) +
  theme(
    legend.position = "left",
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

tree3 <- ggtree(snp_tree, branch.length = "none") %<+% metadata 

tree3$data$branch_color <- GroupColours[tree3$data$Sequence_Type]

tree4 <- tree3 +
  geom_tree(aes(color = I(branch_color)), size = 1) +
  geom_tippoint(aes(color = Sequence_Type), size = NA) +
  geom_tiplab(aes(label=""), align=TRUE) +
  scale_color_manual(
    name = "Sequence Type",
    values = GroupColours
  ) +
  theme(
    legend.position = "left",
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

tree4 <- flip(tree4)

##### Arrange Trees

TreeArranged <- ggarrange(tree2, tree4,
                          nrow = 1,
                          ncol =2,
                          labels = c("A. Core","B. SNP"),
                          common.legend = TRUE,
                          legend = "bottom",
                          heights = c(2,2))

ggsave("Core_SNP_tree.svg",
       TreeArranged,
       height = 25,
       width = 20)
