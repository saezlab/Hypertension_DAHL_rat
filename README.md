# DAHL rat analysis

The following repository contains the scripts used for the analysis of the
phosphorylation data from DAHL rats using
[PHONEMeS-ILP](https://github.com/saezlab/PHONEMeS-ILP).
Here you can find the corresponding article:
> Rinschen MM, Palygin O, Guijas C, Palermo A, Palacio-Escat N, et al.
[Metabolic rewiring of the hypertensive kidney
](https://stke.sciencemag.org/content/12/611/eaax9760). *Science Signaling* 2019. **12**; 611: eaax9760

Note that in order to run the ILP algorithm
you need to obtain a copy of `cplex`'s executable (see PHONEMeS-ILP documentation for more information).

It is assumed that a copy of the [`PHONEMeS-ILP/Public` folder](https://github.com/saezlab/PHONEMeS-ILP/tree/master/Public) is inside the working directory along with `cplex` executable.

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

## Pipeline description

- `1_preprocessing.R`: Cleans the metadata from the data table and generates the p
hosphosite identifiers.
- `2_diff_exp.R`: Performs the differential expression (phosphorylation) analysis.
- `3_phonemes.R`: Runs the PHONEMeS-ILP analysis for multiple time-points. The tar
get nodes used for the analysis are defined in lines 131 and 132.

User can either run each script independently (using Rstudio or another GUI, or from the command line using `$ Rscript <name_of_the_script>`) or execute the whole pipeline with the makefile (`$ make`).
