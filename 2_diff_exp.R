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
# Differential expression analysis
# ================================

library(limma)
library(readr)
library(dplyr)

setwd('/mnt/db/Dropbox/JRC_COMBINE/DAHL_rat_Markus')

in_dir = 'data'
out_dir = 'results'

# Load raw data
data = as.data.frame(read_delim(paste(in_dir, 'data.csv', sep='/'), ',', escape_double=F,
                                trim_ws=T))
rownames(data) = data$ids
data = data[-ncol(data)]

targets = read_csv(paste(in_dir, 'targets.csv', sep='/'))

# Design and contrast matrices
f = factor(targets$condition, levels=unique(targets$condition))
design = model.matrix(~0 + f)
cont.matrix = makeContrasts(d7_con=f7d - fcon,
                            d21_con=f21d - fcon,
                            levels=design)

# Fitting the linear models
pfit = lmFit(data, design)
pfit = contrasts.fit(pfit, cont.matrix)
pfit = eBayes(pfit)

ttops = list(topTable(pfit, coef=1, adjust='fdr', n=nrow(pfit)),
             topTable(pfit, coef=2, adjust='fdr', n=nrow(pfit)))
names(ttops) = colnames(cont.matrix)

for(t in names(ttops)){
    ttop = ttops[t][[1]]
    ttop$ids = rownames(ttop)
    write_csv(ttop, paste(out_dir, paste0('ttop_', t, '.csv'),
                                   sep='/'))
}
