# Harmful.Algae.eDNA
Code needed to reproduce the results of Jacobs-Palmer et at, 2020

## Raw Sequencing Files

Fastq files can be found in these Google Drive Links

  * For the WOAC Cruises samples: https://drive.google.com/drive/folders/1RWxxaM9jMAD_F2Ca61z_nSYccbpJ-RW6?usp=sharing
  
  * For the Intertidal eDNA samples: https://drive.google.com/drive/folders/1XtuDMdlnk9V2acNiPJo5T6O1Tp0KjaCz?usp=sharing
  
## FASTQ file processing

Raw Fastq reads were demultiplexed, quality filtered and denoised using a bioinformatics pipeline that reorients reads for DADA2 denoising. The pipeline can be downloaded from https://github.com/ramongallego/demultiplexer_for_DADA2, and each Miseq run processed with the metadata and parameters file included in each folder.

## From Post processing to the Manuscript data and figures

Follow the main script under `Manuscript/Submission_FrontiersEcolEvol_2020/Jacobs-Palmer_Frontiers_2020.Rmd`

## Supplementary Figures and Tables

The main script generates the datasets included on each of the Tables

