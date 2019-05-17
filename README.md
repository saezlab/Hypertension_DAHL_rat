# DAHL rat analysis

The following repository contains the scripts used for the analysis of the
phosphorylation data from DAHL rats using
[PHONEMeS-ILP](https://github.com/saezlab/PHONEMeS-ILP). Note that in order to run the ILP algorithm
you need to obtain a copy of `cplex`'s executable (see PHONEMeS-ILP documentation for more information).

It is assumed that the contents of the [`PHONEMeS-ILP/Public` folder](https://github.com/saezlab/PHONEMeS-ILP/tree/master/Public) are inside a subfolder on the working directory.

## Dependencies

All the scripts are written in R, the packages needed (apart from the scripts from PHONEMeS-ILP mentioned above) in order to run the complete analysis are:

- readr
- dplyr
- limma
- XML
- UniProt.ws
- igraph
- BioNet
- PHONEMeS

## Kinase-substrate network

The kinase-substrate background network was extracted from OmniPath as of 15/2/2019. The PTM database was queried with the following parameters: `organisms=10116` and `types=phosphorylation,dephosphorylation` [LINK](http://omnipathdb.org/ptms?organisms=10116&types=phosphorylation,dephosphorylation).

## Pipeline

1. `1_preprocessing.R`
2. `2_diff_exp.R`
3. `3_phonemes.R`
