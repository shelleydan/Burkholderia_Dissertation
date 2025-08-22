library(dplyr)
library(tidyr)
library(tidyverse)

theme_set(theme_bw())

### Cleaning Annotation File
AnnotationSum <- read.delim("C:/Users/Dan_PC/OneDrive - Cardiff University/002 - Cardiff University - Postgraduate/MSc Big Data Biology/BIT104_Dissertation/CHP2_GWAS/cysticfibrosis/annotated_summaryCF.tsv", 
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




## Subsetting the data frame to find the genes
# Regions Used

### Clinical GWAS
# Chromosome ID, Start, End
# Chromosome 1, 1.3, 1.34
# Chromosome 1, 1.72, 1.78
# Chromosome 2, 1.14. 1.18
# Chromosome 2, 1.75, 2.05
# Chromosome 2, 2.675, 2.725
# Chromosome 2, 2.780, 2.825
# Chromosome 3, 0.255, 0.270
# Chromosome 3, 0.360, 0.380

### CF GWAS
# Chromosome ID, Start, End
# Chromosome 1, 1.75, 1.92
# Chromosome 2, 0.425, 0.435
# Chromosome 2, 2.650, 2.710
# Chromosome 2, 2.790, 2.810
# Chromosome 2, 2.910, 3.100
# Chromosome 3, 0.350, 0.400

Annotation_subregion <- Annotation_Loci[Annotation_Loci$pos >= 0.425 & Annotation_Loci$pos <= 0.435 & Annotation_Loci$contig == "Chromosome 2" & Annotation_Loci$p_val_log >=7,]

Annotation_subregion$gene_edit

Annotation_subregion %>%
  count(COG_name) %>%
  arrange(desc(n))

# Confirming plot patterns to find the correct region. 
ggplot(Annotation_subregion, aes(x = pos, y = p_val_log)) +
  geom_point(size = 3, alpha = 0.5)

