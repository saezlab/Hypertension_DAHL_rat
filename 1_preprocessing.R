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
# Data pre-processing
# ===================
#
# Data cleaning and plot generation

library(readr)
library(dplyr)

#library(omicToolsTest)

in_dir = 'data'
out_dir = 'results'
ifelse(!dir.exists(out_dir), dir.create(out_dir, recursive=T), F)

# Load raw data
raw_data = read_delim(paste(in_dir, 'P-siteexpression_fornico.txt', sep='/'),
                      '\t', escape_double=F, trim_ws=T)

usecols = colnames(raw_data[(startsWith(colnames(raw_data), 'Intensity') &
                             grepl('.*\\d', colnames(raw_data)))])
cond = as.character(raw_data[2, usecols])
cond[1] = substr(cond[1], regexpr('\\}', cond[1])[1] + 1, nchar(cond[1]))

raw_data = raw_data[3:nrow(raw_data), ]

data = as.data.frame(lapply(raw_data[, usecols], as.numeric))
ids = paste(raw_data$Protein, paste0(raw_data$`Amino acid`, raw_data$Position),
            sep='_')

rownames(data) = ids

# Defining sample conditions
targets = as.data.frame(matrix(NA, nrow=ncol(data), ncol=2))
colnames(targets) = c('sample', 'condition')
targets[, 'sample'] = colnames(data)
targets[, 'condition'] = cond

write.csv(targets, paste(in_dir, 'targets.csv', sep='/'), row.names=F, quote=F)

# Save the results
#magicPlotMaker(data, outpath=out_dir, targets=targets)
data$ids = rownames(data)
write_csv(data, paste(in_dir, 'data.csv', sep='/'))
