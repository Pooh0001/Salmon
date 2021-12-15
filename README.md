# Mosquito viral infection

## Setup
```{bash}
ln -s /mnt/picea/projects/labs/datafolder/ data
```

Describe the R scripts used to analyse RNA-seq data, that have been published DOI:

Short description of study: RNA-seq data of Aedes albopictus cells infected with either West Nile Virus (WNV), insect-specific Lammi virus (LamV) or an infection scheme whereby cells were pre-infected with LamV prior to challenge with WNV. Our aim is to improve the understanding of the insect-specific flavivirus interaction and its effects on the host. Furthermore, we investigated the possible interfering effects of a prior LamV infection on WNV replication.

Raw sequencing data have been made publicly available at NCBI Sequence Read Archive (SRA) under accession number Bioproject ID PRJNA785413.

We strongly recommend to use the RNA sequencing analysis pipeline nf-core/rnaseq to process the raw data (https://nf-co.re/rnaseq). The pipeline included quality control (FastQC v0.11.9), adapter and quality trimming (Trim Galore! V0.6.6), removal of ribosomal RNA (SortMeRNA v4.2.0) and pseudo-alignment and quantification (Salmon v1.4.0) using the reference Ae. albopictus Foshan genome sequence and annotation (version Aalo1.2) retrieved from Vectorbase https://vectorbase.org/vectorbase/app/record/dataset/DS_58c436b555#category:data-identity-and-mapping.

The resulting Salmon quant.sf files should first be analysed with the R script named QC_Biological_QA.R, which includes a quality control of the data to see that the sequencing depth is homogeneous, that the samples clustered as expected when analyzed with a Principal-Component Analysis (PCA) and if it is enough variance in the data to expect differential expression. For this analysis you need to give a CSV file with samples information, an example of such CSV file can be found here. For later analysis you need full rank (data files in all fields of the matrix), you may remove or add "fake" data points to achieve this.

For the differential expression analysis using DESeq2, the R script DE_all_in_one.R should be used. This script compare the altered gene expression upon infection at different time points compared to non-infected cells at the same time points. The DE expression is visualized in different ways such as heatmaps and volcano plots. It also generates two tables with log2fold changes and p-values, those with only the DE transcripts (named genes) and and a table with expression data of all transcripts (named results).
