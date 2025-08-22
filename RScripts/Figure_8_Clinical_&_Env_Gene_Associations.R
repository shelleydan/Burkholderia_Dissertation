### Library
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rtracklayer)
library(tidyr)
library(ggrepel)

theme_set(theme_bw())

### Data

### Cleaning Annotation File
AnnotationSum <- read.delim("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/clinical/annotated_summary.tsv", 
                            sep = "\t",
                            header = TRUE,
                            quote = "")


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

top30 <- Annotation_Loci %>%
  arrange(desc(p_val_log)) %>%
  slice_head(n = 30)

### Volcano Plots 

ggplot(Annotation_Loci, aes(x = avg.beta, 
                            y = p_val_log, 
                            size = unitigs,
                            color = pangenome.category)) +
  geom_point(alpha = 0.6) +
  scale_size_continuous(range = c(1, 8)) +   # adjust dot size scaling
  labs(
    x = "Average Beta",
    y = "-log10(p-value)",
    size = "Unitigs",
    color = "Pangenome Category"
  ) +
  geom_text_repel(data = top30,
                  aes(label = label),   # change 'gene_name' to your column
                  size = 3,
                  max.overlaps = 50,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color = "grey50",
                  show.legend = FALSE) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             colour = "grey50") +
  geom_hline(yintercept = 7, 
             linetype = "dashed", 
             colour = "brown2") +
  scale_color_manual(values = c("core" = "grey", 
                                "soft-core" = "brown2",
                                "shell" = "orange",
                                "cloud" = "black"),
                     name = "Pangenome Category")


ggsave("ClinicalVolcano.svg",
       height = 5,
       width = 7.5)
