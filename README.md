# Harmful.Algae.eDNA
Code needed to reproduce the results of Jacobs-Palmer et at, 2020

## Raw Sequencing Files

Fastq files can be found in these Google Drive Links

  * For the WOAC Cruises samples: https://drive.google.com/drive/folders/1RWxxaM9jMAD_F2Ca61z_nSYccbpJ-RW6?usp=sharing
  
  * For the Intertidal eDNA samples: https://drive.google.com/drive/folders/1XtuDMdlnk9V2acNiPJo5T6O1Tp0KjaCz?usp=sharing
  
## FASTQ file processing

Raw Fastq reads were demultiplexed, quality filtered and denoised using a bioinformatics pipeline that reorients reads for DADA2 denoising. The pipeline can be downloaded from https://github.com/ramongallego/demultiplexer_for_DADA2, and each Miseq run processed with the metadata and parameters file included in each folder.

## From Post processing to the Manuscript

We removed artifactal sequences, tag-jumping, and other contamination by using a custom script. The Rmarkdown that replicatets the process for the WOAC cruises samples is  `Manuscript/data/WOAC_Run_20190201` and produces the output `Copy_of_WOAC_ASV.csv` and `WOAC_Hash_Key_all_together.csv`. An identical cleaning process for the Intertidal eDNA samples is detailed in the repository for that original publication https://github.com/ramongallego/eDNA.and.Ocean.Acidification.Gallego.et.al.2020.

## Secondary clustering with Swarm

DADA2 and clean-up got rid of sequencing inaccuracies - but the PCR cycles previous to that could have created artifacts - we dealt with that with secondary clustering using a snow-balling algorithm (swarm, https://github.com/torognes/swarm). The code, input and output is in `Manuscript/data/swarm`. You need a local installation of  `swarm` to run this.

## Taxonomical Annotation

The starting point of the taxonomical annotation of these sequences is the output of the CRUX pipeline with the Leray Primers. The Classification file is `Copy_of_hash.annotated.csv`. 

## Creating the manuscript, analysis, Tables and Tigures

Follow the main script under `Manuscript/Submission_FrontiersEcolEvol_2020/Jacobs-Palmer_Frontiers_2020.Rmd`.
To knit your own copy of the manuscript, figures and results, you need to clone this repository and run the above Rmarkdown. You need these Rpackages installed:
   
   *"tidyverse",
   *"lubridate",
   * "vegan",
   * "seacarb",
   * "modelr",
   *"broom",
   *"pander",
   *"gridExtra",
   * "here",
   * "tidymodels",
   *"rstanarm",
   *"loo",
   *"tidybayes",
   *"bayesplot"

## Supplementary Figures and Tables

The main script generates the datasets included on each of the Tables, and the final tables are rendered in the File `/Manuscript/Supplement/Jacobs-Palmer_Frontiers_2020_Supp.Rmd`

