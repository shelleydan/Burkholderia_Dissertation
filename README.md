# BIT104: Burkholderia Dissertation
MSc Big Data Biology @ Cardiff University
<br>
#### This repository supports the analysis which contributes to a dissertation submission for BIT104 of MSc Big Data Biology at Cardiff University. This dissertation is titled: A comparative genomic interrogation of _Burkholderia cenocepacia_ as a cystic fibrosis lung pathogen  and rarely encountered environmental bacterium. Here, we aimed to indentify key associations of _B. cenocepacia_ with the environment and CF-derived strains. We reclassified the _B. cenocepacia_ databse to identify true _B. cenocepacia_ isolates. We identified that 37.6% of the database was misclassified. Approches used to analyse these genomes are outlines below with links to the relevant scripts used. Any scripts not self-curated are appropriately referenced in their relevant methods sections. 

# Overview of Genome Curation


<img width="8865" height="4757" alt="Flowchart" src="https://github.com/user-attachments/assets/a4799682-0fe4-41c7-bad0-1547f214db9e" />


# Script Summary

<table>
  <thead>
    <tr>
      <th>Script</th>
      <th>Description</th>
      <th>Tool Version</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/SRA-RawReads-Download.sbatch" target="_blank">SRA.sbatch</a></td>
      <td>The SRA toolkit was used to download genomes using accession IDs, specific tools included prefetch and fasterq-dump.</td>
      <td>3.2.1</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/QC_TRIM.sbatch" target="_blank">QC_TRIM.sbatch</a></td>
      <td>Two tools; FastQC was used to generate quality reports prior and following trimming, FastP was used to perform trimming.</td>
      <td>FastQC; 0.12.1 <br> FastP; 0.23.4</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/unicyler.sbatch" target="_blank">Unicycler.sbatch</a></td>
      <td>Unicycler was used to assemble genomes from trimmed reads.</td>
      <td>0.5.0</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/quast.sbatch" target="_blank">QUAST.sbatch</a></td>
      <td>Generates a table including contig information.</td>
      <td>5.3.0</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/checkM2.sbatch" target="_blank">CheckM2.sbatch</a></td>
      <td>Generates a table of genome information, including completeness, contamination etc.</td>
      <td>1.10</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/tree/main/nextflow/modules/FastANI" target="_blank">FastANI.nf</a></td>
      <td>A self-curated FastANI module for nextflow was generated to generate an average nucleotide identity matrix.</td>
      <td>1.34</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/mlst.sbatch" target="_blank">MLST.sbatch</a></td>
      <td>Performs MLST analysis.</td>
      <td>2.13.0</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/abricate.sbatch" target="_blank">ABRicate.sbatch</a></td>
      <td>Provides annotation of antimicrobial resistance genes and virulence factors.</td>
      <td>1.0.0</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/bakta.sbatch" target="_blank">Bakta.sbatch</a>**</td>
      <td>Provide high quality annotations - used to confirm some ABRicate results</td>
      <td>1.11.3</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/prokka.sbatch" target="_blank">Prokka.sbatch</a></td>
      <td>Provide whole genome annotations.</td>
      <td>1.14.6</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/core_gene_tree.sbatch" target="_blank">Panaroo.sbatch</a><br><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/SNP.sbatch" target="_blank">SNP.sbatch</a></td>
      <td>Provide phylogenetic analysis, based on core genes and single nucleotide polymorphisms.</td>
      <td>Panaroo; 1.2.10 <br> Snippy; 4.6.0 <br> SNP-Sites; 2.5.1 <br> RAxML; 1.2.0</td>
    </tr>
    <tr>
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/microGWAS.sbatch" target="_blank">MicroGWAS.sbatch</a>**</td>
      <td>Perform a Genome-Wide Association Study.</td>
      <td>0.6.0</td>
    </tr>
  </tbody>
</table>
** A Conda install of MicroGWAS and Bakta was used, information about the installation method can be found here: <a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/conda/conda_module_HPC.txt" target="_blank">Conda HPC Module Install</a>

# Figure & Table Curation

<table>
  <thead>
    <tr>
      <th>Figure</th>
      <th>Script</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Figure 1</td>
      <td>Figure 1 was generated using Microsoft Powerpoint following quality filtering. <br> See Overview of Genome Curation</td>
    </tr>
    <tr>
      <td>Figure 2</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_2_Metadata_Visualisations.R" target="_blank">Figure_2_Metadata_Visualisations.R</a></td>
    </tr>
    <tr>
      <td>Figure 3</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_3_Full_Pairwise_Heatmap.R" target="_blank">Figure_3_Full_Pairwise_Heatmap.R</a></td>
    </tr>
    <tr>
      <td>Figure 4</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_4_Cenocepacia_Pairwise_Heatmap.R" target="_blank">Figure_4_Cenocepacia_Pairwise_Heatmap.R</a></td>
    </tr>
    <tr>
      <td>Figure 5</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_5_SNP_tree_with_Heatmaps.R" target="_blank">Figure_5_SNP_tree_with_Heatmaps.R</a></td>
    </tr>
    <tr>
      <td>Figure 6</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_6_Pangenome_SNP_Trees.R" target="_blank">Figure_6_Pangenome_SNP_Trees.R</a></td>
    </tr>
    <tr>
      <td>Figure 7</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_7_Clinical_%26_Env_Manhattan.R" target="_blank">Figure_7_Clinical_&_Env_Manhattan.R</a></td>
    </tr>
    <tr>
      <td>Figure 8</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_8_Clinical_%26_Env_Gene_Associations.R" target="_blank">Figure_8_Clinical_&_Env_Gene_Associations.R</a></td>
    </tr>
    <tr>
      <td>Figure 9</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_9_CF_Manhattan.R" target="_blank">Figure_9_CF_Manhattan.R</a></td>
    </tr>
    <tr>
      <td>Figure 10</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Figure_10_CF_Gene_Association.R" target="_blank">Figure_10_CF_Gene_Association.R</a></td>
    </tr>
    <tr>
      <td>Supplementary Figure 1</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Sup_Figure_1_Sequence_Type_Distribution.R" target="_blank">Sup_Figure_1_Sequence_Type_Distribution.R</a></td>
    </tr>
    <tr>
      <td>Supplementary Tables 1 & 2</td>
      <td><a href="https://github.com/shelleydan/Burkholderia_Dissertation/blob/main/RScripts/Sup_Tables_Annotations.R" target="_blank">Sup_Tables_Annotations.R</a></td>
    </tr>
  </tbody>
</table>

# RStudio Session.info()
```
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8   
[3] LC_MONETARY=English_United Kingdom.utf8 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.utf8    

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] phytools_2.4-4        maps_3.4.2.1          ape_5.8-1             Polychrome_1.3.1     
 [5] cowplot_1.1.3         ggnewscale_0.5.0      ggtree_3.10.1         ComplexHeatmap_2.18.0
 [9] ggrepel_0.9.6         ggpubr_0.6.0          data.table_1.16.4     lubridate_1.9.4      
[13] forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
[17] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1        
[21] tidyverse_2.0.0      

loaded via a namespace (and not attached):
 [1] mnormt_2.1.1            phangorn_2.12.1         rlang_1.1.4            
 [4] magrittr_2.0.3          clue_0.3-66             GetoptLong_1.0.5       
 [7] matrixStats_1.5.0       compiler_4.3.2          png_0.1-8              
[10] vctrs_0.6.5             combinat_0.0-8          quadprog_1.5-8         
[13] pkgconfig_2.0.3         shape_1.4.6.1           crayon_1.5.3           
[16] backports_1.5.0         tzdb_0.4.0              aplot_0.2.4            
[19] clusterGeneration_1.3.8 jsonlite_2.0.0          broom_1.0.7            
[22] parallel_4.3.2          cluster_2.1.4           R6_2.6.1               
[25] stringi_1.8.7           RColorBrewer_1.1-3      car_3.1-3              
[28] numDeriv_2016.8-1.1     Rcpp_1.0.13-1           iterators_1.0.14       
[31] optimParallel_1.0-2     IRanges_2.36.0          Matrix_1.6-1.1         
[34] igraph_2.1.4            timechange_0.3.0        tidyselect_1.2.1       
[37] rstudioapi_0.17.1       abind_1.4-8             doParallel_1.0.17      
[40] codetools_0.2-19        lattice_0.21-9          treeio_1.26.0          
[43] withr_3.0.2             coda_0.19-4.1           gridGraphics_0.5-1     
[46] circlize_0.4.16         pillar_1.10.1           carData_3.0-5          
[49] foreach_1.5.2           stats4_4.3.2            ggfun_0.1.8            
[52] generics_0.1.3          S4Vectors_0.40.2        hms_1.1.3              
[55] munsell_0.5.1           scales_1.3.0            tidytree_0.4.6         
[58] glue_1.8.0              scatterplot3d_0.3-44    lazyeval_0.2.2         
[61] tools_4.3.2             ggsignif_0.6.4          fs_1.6.5               
[64] fastmatch_1.1-6         colorspace_2.1-1        nlme_3.1-163           
[67] patchwork_1.3.0         Formula_1.2-5           cli_3.6.3              
[70] DEoptim_2.2-8           expm_1.0-0              gtable_0.3.6           
[73] rstatix_0.7.2           yulab.utils_0.1.9       digest_0.6.37          
[76] BiocGenerics_0.48.1     ggplotify_0.1.2         rjson_0.2.23           
[79] farver_2.1.2            lifecycle_1.0.4         GlobalOptions_0.1.2    
[82] MASS_7.3-60  
```
