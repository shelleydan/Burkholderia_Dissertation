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
      <th>Dependencies</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Figure 1</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 2</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 3</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 4</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 5</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 6</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 7</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 8</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 9</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Figure 10</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Supplementary Fingure 1</td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>
# Package and Tool Versions
