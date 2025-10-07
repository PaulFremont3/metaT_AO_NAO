#!/bin/env/usr/env Rscript

frac = commandArgs(trailingOnly = T)
v <- read.table(paste('metat_analysis_',frac,'.txt', sep=''))
saveRDS(v, paste('metat_analysis_',frac,'.rds', sep=''))
