# BIT104: Burkholderia Dissertation
MSc Big Data Biology @ Cardiff University

# Script Function

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
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/bakta.sbatch" target="_blank">Bakta.sbatch</a></td>
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
      <td><a href="https://github.com/shelleydan/Llyfrgell_Bersonol/blob/main/bin/microGWAS.sbatch" target="_blank">MicroGWAS.sbatch</a></td>
      <td>Perform a Genome-Wide Association Study.</td>
      <td>0.6.0</td>
    </tr>
  </tbody>
</table>


# Figure & Table Curation

# Package and Tool Versions
