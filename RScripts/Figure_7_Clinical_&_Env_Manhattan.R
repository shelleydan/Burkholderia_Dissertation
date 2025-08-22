library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rtracklayer)
library(tidyr)
library(ggrepel)

theme_set(theme_bw())

# Read in the Unitig data from microGWAS
mappedregions <- read.table("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/clinical/mapped_all.tsv",
                            sep = "\t",
                            header = TRUE)

unique(mappedregions$strain)
unique(mappedregions$contig)

mappedregions <- mappedregions %>% filter(strain == "J2315")

mappedregions <- mappedregions %>% filter(notes != "bad-chisq")

mappedregions$lrt_pval <- as.numeric(mappedregions[["lrt.pvalue"]])

mappedregions$pos <- (mappedregions$start + mappedregions$end) / (2 * 1000000)
mappedregions$pos <- as.numeric(mappedregions$pos)

mappedregions$p_val_log <- -log10(mappedregions$lrt_pval)
mappedregions$p_val_log <- as.numeric(mappedregions$p_val_log)

# Assign significance and which way associations lie. 
mappedregions$colour_group <- with(mappedregions, ifelse(p_val_log >= 7 & beta < 0, "Significant_Environment",
                                                         ifelse(p_val_log >= 7 & beta > 0, "Significant_Clinical",
                                                                "Non-Significant")))
unique(mappedregions$colour_group)

# Rename contigs as chromosomes
mappedregions <- mappedregions %>%
  mutate(contig = recode(contig,
                         "NC_011000.1" = "Chromosome 1",
                         "NC_011001.1" = "Chromosome 2",
                         "NC_011002.1" = "Chromosome 3",
                         "NC_011003.1" = "Plasmid"))

facet_labels <- c(
  "NC_011000.1" = "Chromosome 1",
  "NC_011001.1" = "Chromosome 2",
  "NC_011002.1" = "Chromosome 3",
  "NC_011003.1"  = "Plasmid"
)

######## Prep Annotations

# Read in significant genes
AnnotationSum <- read.delim("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/clinical/annotated_summary.tsv", 
                            sep = "\t",
                            header = TRUE,
                            quote = )

# Read in conversion table to match annotations and unitigs
SNP_Loci <- read.table("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/cysticfibrosis/regions.tsv", 
                       sep = "\t",
                       header = FALSE)

colnames(SNP_Loci)[colnames(SNP_Loci) == "V1"] <- "J2315"
colnames(SNP_Loci)[colnames(SNP_Loci) == "V2"] <- "Start_End"

Annotation_Loci <- merge(AnnotationSum, SNP_Loci, by = "J2315")

Annotation_Loci <- Annotation_Loci %>%
  separate(Start_End, into = c("contig", "coords"), sep = ":") %>%
  separate(coords, into = c("Start", "End"), sep = "-")

Annotation_Loci$Start <- as.numeric(Annotation_Loci$Start)
Annotation_Loci$End <- as.numeric(Annotation_Loci$End)

Annotation_Loci$gene_edit <- sub("~~~.*", "", Annotation_Loci$gene)

Annotation_Loci$pos <- (Annotation_Loci$Start + Annotation_Loci$End) / (2 * 1000000)

Annotation_Loci$p_val_log <- -log10(Annotation_Loci$avg.lrt.pvalue)
Annotation_Loci$p_val_log <- as.numeric(Annotation_Loci$p_val_log)

Annotation_Loci <- Annotation_Loci %>%
  arrange(desc(p_val_log)) %>%
  mutate(label = ifelse(row_number() <= 20, gene_edit, NA))

Annotation_Loci <- Annotation_Loci %>%
  mutate(contig = recode(contig,
                         "NC_011000" = "Chromosome 1",
                         "NC_011001" = "Chromosome 2",
                         "NC_011002" = "Chromosome 3",
                         "NC_011003" = "Plasmid"))

##########
# Plot the full manhattan and divide by chromosome
FullManhatten <- ggplot(mappedregions, aes(x = pos, y = p_val_log, colour = colour_group)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_colour_manual(values = c("Non-Significant" = "grey", 
                                 "Significant_Clinical" = "cornflowerblue",
                                 "Significant_Environment" = "darkseagreen2"),
                      name = "Significance") +
  geom_hline(yintercept = 7, 
             linetype = "dashed",
             color = "brown2",
             linewidth = 1) +
  labs( y = expression(italic("-log"[10]*"(pvalue)")),
        x = "Genome Position (Mb)") +
  #scale_x_continuous(breaks = seq(0,4, by = 0.5),
  #                  expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,30, by = 5),
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 30)) +
  facet_grid(~ contig, scales = "free_x",
             labeller = labeller(contig = facet_labels))


FullAnnMan <- FullManhatten +
  geom_point(
    data = Annotation_Loci,
    aes(x = pos, y = p_val_log, fill = pangenome.category),
    inherit.aes = FALSE,
    shape = 21,  # shape that supports fill and outline color
    colour = "black",  # outline color for points
    size = 3,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("core" = "grey", 
                               "soft-core" = "brown2",
                               "shell" = "orange",
                               "cloud" = "black"),
                    name = "Pangenome Category")

### Create Zoom

subregion <- mappedregions[mappedregions$pos >= 1 & mappedregions$pos <= 2 & mappedregions$contig == "Chromosome 1",]
subregion$colour_group <- ifelse(mappedregions$p_val_log > 7, "Significant", "Non-Significant")
unique(subregion$colour_group)

Annotation_subregion <- Annotation_Loci[Annotation_Loci$pos >= 1 & Annotation_Loci$pos <= 2 & Annotation_Loci$contig == "Chromosome 1",]

region1 <- ggplot(subregion, aes(x = pos, y = p_val_log, colour = colour_group)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_colour_manual(values = c("Non-Significant" = "grey", 
                                 "Significant_Clinical" = "cornflowerblue",
                                 "Significant_Environment" = "darkseagreen2"),
                      name = "Significance") +
  geom_hline(yintercept = 7, 
             linetype = "dashed",
             color = "brown2",
             linewidth = 1) +
  labs( y = expression(italic("-log"[10]*"(pvalue)")),
        x = "Genome Position (Mb)") +
  coord_cartesian(ylim = c(0, 30),
                  xlim = c(1,2)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,30, by = 5),
                     expand = c(0, 0))+
  ggtitle("Chromosome 1")+
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

chr1_region <- region1 +
  geom_point(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, fill = pangenome.category),
    inherit.aes = FALSE,
    shape = 21,  # shape that supports fill and outline color
    colour = "black",  # outline color for points
    size = 3,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("core" = "grey", 
                               "soft-core" = "brown2",
                               "shell" = "orange",
                               "cloud" = "black"),
                    name = "Pangenome Category") +
  geom_text_repel(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, label = gene_edit),
    inherit.aes = FALSE,
    size = 3,
    colour = "black",
    max.overlaps = 10  # adjust or remove to control label hiding
  )

### region 2
subregion <- mappedregions[mappedregions$pos >= 1.1 & mappedregions$pos <= 1.3 & mappedregions$contig == "Chromosome 2",]
subregion$colour_group <- ifelse(mappedregions$p_val_log > 7, "Significant", "Non-Significant")
unique(subregion$colour_group)

Annotation_subregion <- Annotation_Loci[Annotation_Loci$pos >= 1.1 & Annotation_Loci$pos <= 1.3 & Annotation_Loci$contig == "Chromosome 2",]

region1 <- ggplot(subregion, aes(x = pos, y = p_val_log, colour = colour_group)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_colour_manual(values = c("Non-Significant" = "grey", 
                                 "Significant_Clinical" = "cornflowerblue",
                                 "Significant_Environment" = "darkseagreen2"),
                      name = "Significance") +
  geom_hline(yintercept = 7, 
             linetype = "dashed",
             color = "brown2",
             linewidth = 1) +
  labs( y = expression(italic("-log"[10]*"(pvalue)")),
        x = "Genome Position (Mb)") +
  coord_cartesian(ylim = c(0, 30),
                  xlim = c(1.1,1.3)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,30, by = 5),
                     expand = c(0, 0))+
  ggtitle("Chromosome 2")+
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

chr2_region1 <- region1 +
  geom_point(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, fill = pangenome.category),
    inherit.aes = FALSE,
    shape = 21,  # shape that supports fill and outline color
    colour = "black",  # outline color for points
    size = 3,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("core" = "grey", 
                               "soft-core" = "brown2",
                               "shell" = "orange",
                               "cloud" = "black"),
                    name = "Pangenome Category") +
  geom_text_repel(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, label = gene_edit),
    inherit.aes = FALSE,
    size = 3,
    colour = "black",
    max.overlaps = 10  # adjust or remove to control label hiding
  )

## Chromosome 2, region 2

subregion <- mappedregions[mappedregions$pos >= 1.5 & mappedregions$pos <= 2.5 & mappedregions$contig == "Chromosome 2",]
subregion$colour_group <- ifelse(mappedregions$p_val_log > 7, "Significant", "Non-Significant")
unique(subregion$colour_group)

Annotation_subregion <- Annotation_Loci[Annotation_Loci$pos >= 1.5 & Annotation_Loci$pos <= 2.5 & Annotation_Loci$contig == "Chromosome 2",]

region1 <- ggplot(subregion, aes(x = pos, y = p_val_log, colour = colour_group)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_colour_manual(values = c("Non-Significant" = "grey", 
                                 "Significant_Clinical" = "cornflowerblue",
                                 "Significant_Environment" = "darkseagreen2"),
                      name = "Significance") +
  geom_hline(yintercept = 7, 
             linetype = "dashed",
             color = "brown2",
             linewidth = 1) +
  labs( y = expression(italic("-log"[10]*"(pvalue)")),
        x = "Genome Position (Mb)") +
  coord_cartesian(ylim = c(0, 30),
                  xlim = c(1.5,2.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,30, by = 5),
                     expand = c(0, 0))+
  ggtitle("Chromosome 2")+
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

chr2_region2 <- region1 +
  geom_point(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, fill = pangenome.category),
    inherit.aes = FALSE,
    shape = 21,  # shape that supports fill and outline color
    colour = "black",  # outline color for points
    size = 3,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("core" = "grey", 
                               "soft-core" = "brown2",
                               "shell" = "orange",
                               "cloud" = "black"),
                    name = "Pangenome Category") +
  geom_text_repel(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, label = gene_edit),
    inherit.aes = FALSE,
    size = 3,
    colour = "black",
    max.overlaps = 10  # adjust or remove to control label hiding
  )

## Chromosome 2, region 3

subregion <- mappedregions[mappedregions$pos >= 2.5 & mappedregions$pos <= 3 & mappedregions$contig == "Chromosome 2",]
subregion$colour_group <- ifelse(mappedregions$p_val_log > 7, "Significant", "Non-Significant")
unique(subregion$colour_group)

Annotation_subregion <- Annotation_Loci[Annotation_Loci$pos >= 2.5 & Annotation_Loci$pos <= 3 & Annotation_Loci$contig == "Chromosome 2",]

region1 <- ggplot(subregion, aes(x = pos, y = p_val_log, colour = colour_group)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_colour_manual(values = c("Non-Significant" = "grey", 
                                 "Significant_Clinical" = "cornflowerblue",
                                 "Significant_Environment" = "darkseagreen2"),
                      name = "Significance") +
  geom_hline(yintercept = 7, 
             linetype = "dashed",
             color = "brown2",
             linewidth = 1) +
  labs( y = expression(italic("-log"[10]*"(pvalue)")),
        x = "Genome Position (Mb)") +
  coord_cartesian(ylim = c(0, 30),
                  xlim = c(2.5,3)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,30, by = 5),
                     expand = c(0, 0))+
  ggtitle("Chromosome 2")+
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

chr2_region3 <- region1 +
  geom_point(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, fill = pangenome.category),
    inherit.aes = FALSE,
    shape = 21,  # shape that supports fill and outline color
    colour = "black",  # outline color for points
    size = 3,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("core" = "grey", 
                               "soft-core" = "brown2",
                               "shell" = "orange",
                               "cloud" = "black"),
                    name = "Pangenome Category") +
  geom_text_repel(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, label = gene_edit),
    inherit.aes = FALSE,
    size = 3,
    colour = "black",
    max.overlaps = 10  # adjust or remove to control label hiding
  )

## Chromosome 3, region 1

subregion <- mappedregions[mappedregions$pos >= 0.2 & mappedregions$pos <= 0.4 & mappedregions$contig == "Chromosome 3",]
subregion$colour_group <- ifelse(mappedregions$p_val_log > 7, "Significant", "Non-Significant")
unique(subregion$colour_group)

Annotation_subregion <- Annotation_Loci[Annotation_Loci$pos >= 0.2 & Annotation_Loci$pos <= 0.4 & Annotation_Loci$contig == "Chromosome 3",]

region1 <- ggplot(subregion, aes(x = pos, y = p_val_log, colour = colour_group)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_colour_manual(values = c("Non-Significant" = "grey", 
                                 "Significant_Clinical" = "cornflowerblue",
                                 "Significant_Environment" = "darkseagreen2"),
                      name = "Significance") +
  geom_hline(yintercept = 7, 
             linetype = "dashed",
             color = "brown2",
             linewidth = 1) +
  labs( y = expression(italic("-log"[10]*"(pvalue)")),
        x = "Genome Position (Mb)") +
  coord_cartesian(ylim = c(0, 30),
                  xlim = c(0.2,0.4)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,30, by = 5),
                     expand = c(0, 0))+
  ggtitle("Chromosome 3") +
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

chr3_region1 <- region1 +
  geom_point(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, fill = pangenome.category),
    inherit.aes = FALSE,
    shape = 21,  # shape that supports fill and outline color
    colour = "black",  # outline color for points
    size = 3,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("core" = "grey", 
                               "soft-core" = "brown2",
                               "shell" = "orange",
                               "cloud" = "black"),
                    name = "Pangenome Category") +
  geom_text_repel(
    data = Annotation_subregion,
    aes(x = pos, y = p_val_log, label = gene_edit),
    inherit.aes = FALSE,
    size = 3,
    colour = "black",
    max.overlaps = 10  # adjust or remove to control label hiding
  )

#### ggarrange

# Arrange Chromosome 2 Zoom ins
middle_row <- ggarrange(chr2_region1, chr2_region2, chr2_region3, 
                        labels = c("B", "C", "D"), 
                        ncol = 3, 
                        common.legend = FALSE, 
                        legend = "none")

# Arrange chromosome 1 and 3 zoom ins
bottom_row <- ggarrange(chr1_region, chr3_region1, 
                        labels = c("E", "F"), 
                        ncol = 2, 
                        common.legend = FALSE, 
                        legend = "none")

# Arrange all plots together
final_plot <- ggarrange(FullAnnMan, middle_row, bottom_row,
                        nrow = 3,
                        labels = "A",
                        common.legend = TRUE,
                        legend = "bottom",
                        heights = c(2,2,2))


ggsave("ManhattenClinical.png",
       plot = final_plot,
       height = 15,
       width = 15)
