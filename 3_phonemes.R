# Copyright (C) 2019 Nicolàs Palacio
#
# Contact: palacio@combine.rwth-aachen.de
#
# GNU-GLPv3:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# A full copy of the GNU General Public License can be found on
# http://www.gnu.org/licenses/.
#
# PHONEMeS-ILP analysis
# =====================

library(PHONEMeS)
library(readr)
library(dplyr)
library(UniProt.ws)
library(igraph)
library(BioNet)

setwd('/mnt/DB/Dropbox/JRC_COMBINE/DAHL_rat_Markus')

# Loading PHONEMeS-ILP scripts
sdir = 'PHONEMeS_ILP'
scripts = list.files(sdir)
scripts = scripts[endsWith(scripts, '.R')]

for(s in scripts){
  source(paste(sdir, s, sep='/'))
}

in_dir = 'data'
out_dir = 'results'

# Loading UniProt web services
up = UniProt.ws(taxId=10116)

# Functions
convert = function(upids, up){
    res = select(up, upids, c('UNIPROTKB', 'GENES'), 'UNIPROTKB')
    res$GENES = stringr::str_to_upper(gsub('\\s.*', '', res$GENES))
    
    return(unique(res))
}

###############################################################################
# Building the background network table
ksn = read_tsv(paste(in_dir,
                     'omnipath_webservice_ptms.csv',
                     #'ks_net.txt',
                     sep='/'))

bn = as.data.frame(ksn[-ncol(ksn)])
colnames(bn) = c('K.AC', 'S.AC', 'res', 'pos')

# Converting UniProtIDs to GeneSymbols
for(i in c('K', 'S')){
    ac = paste(i, 'AC', sep='.')
    aux = convert(unique(bn[, ac]), up)
    
    colnames(aux) = c(ac, paste(i, 'ID', sep='.'))
    bn = left_join(bn, aux, by=ac)
}

bn['SID'] = paste0(rep('e', times=nrow(bn)), 1:nrow(bn))
bn['S.cc'] = paste(bn[, 'S.ID'], paste0(bn[, 'res'], bn[, 'pos']), sep='_')
bn = bn[, c("S.AC", "S.ID", "K.AC", "K.ID", "res",  "pos",  "SID",  "S.cc")]

bgn = new("KPSbg", interactions=bn, species=unique(c(bn$K.ID, bn$S.cc)))

###############################################################################
# Loading differential expression data
files = list.files(out_dir)
files = files[startsWith(files, 'ttop')]
ttops = list()

for(f in files){
    ttop = as.data.frame(read_csv(paste(out_dir, f, sep='/')))
    ttop = ttop[, c(colnames(ttop)[ncol(ttop)],
                    colnames(ttop)[1:ncol(ttop) - 1])]
    rownames(ttop) = ttop$ids
    
    name = substr(f, 6, nchar(f) - 4)
    ttops[[name]] = ttop
}

bn <- bn[complete.cases(bn), ]
proteins <- c(bn$S.AC, bn$K.AC)
genes <- c(bn$S.ID, bn$K.ID)

for(ii in 1:length(ttops)){
  
  for(jj in 1:nrow(ttops[[ii]])){
    
    prot <- strsplit(x = ttops[[ii]][jj, 1], split = "_", fixed = TRUE)[[1]][1]
    residue <- strsplit(x = ttops[[ii]][jj, 1], split = "_", fixed = TRUE)[[1]][2]
    
    idx <- which(proteins==prot)
    if(length(idx)>0){
      # print(idx)
      ttops[[ii]][jj, 1] <- paste0(genes[idx[1]], "_", residue)
      
    }
    
  }
  
}

for(ii in 1:length(ttops)){
  
  rownames(ttops[[ii]]) <- ttops[[ii]][, 1]
  
}

bn$S.AC <- bn$S.ID
bn$K.AC <- bn$K.ID

###############################################################################
# Building GMM object
#ggm_list = buildDataObject(ttList=ttops, pThresh=0.1, pValIdx=6, measIdx=1,
#                          fcIdx=2,
#                          organism="RAT")
# bn$dataID = paste(bn[, 'S.AC'], paste0(bn[, 'res'], bn[, 'pos']), sep='_')
gmm_list = buildInputs(tableTopList=ttops, fcThresh=NULL, pThresh=c(0.1, 0.1),
                       idxID=1, idxFC=2, idxPval=5, mappingTable=NULL,
                       namesConditions=names(ttops))

#load('dataGMM.RData')
dataGMM = new("GMMres", res=gmm_list$res, IDmap=gmm_list$IDmap,
              resFC=gmm_list$resFC)
###############################################################################
# Setting up and running PHONEMeS
conditions = list(c("d21_con"), c("d7_con"))

names(conditions) = c("d21", "d7")

targets.P = list(d21=c("PRKAA1", "PRKAA2", "MTOR"),
                 d7=c("PRKAA1", "PRKAA2", "MTOR"))

experiments = list(tp1=c(2), tp2=c(1))

#Generating the networks

############# FIXME v: error in next line


tpSIF = runPHONEMeS_dt(targets.P=targets.P, conditions=conditions,
                       dataGMM=dataGMM, experiments=experiments, bg=bgn,
                       nIter=100)

write.table(x=tpSIF, file="pdgfr_tp_analysis.txt", quote=F, sep="\t",
            row.names=F)

# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=5), ], file = "pdgfr_tp_analysis_cutoff_5.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=10), ], file = "pdgfr_tp_analysis_cutoff_10.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=20), ], file = "pdgfr_tp_analysis_cutoff_20.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=50), ], file = "pdgfr_tp_analysis_cutoff_50.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=60), ], file = "pdgfr_tp_analysis_cutoff_60.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=70), ], file = "pdgfr_tp_analysis_cutoff_70.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(x = tpSIF[which(as.numeric(tpSIF[, 2])>=80), ], file = "pdgfr_tp_analysis_cutoff_80.txt", quote = FALSE, sep = "\t", row.names = FALSE)


nodesAttributes = assignAttributes_tp(sif=tpSIF, dataGMM=dataGMM,
                                      targets=targets.P[[1]], writeAttr=T)
