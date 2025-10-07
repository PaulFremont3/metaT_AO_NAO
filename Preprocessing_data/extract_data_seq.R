#!/bin/env/usr/env Rscript
library("gbm")
library("dismo")
# library("FactoMineR")
# library("factoextra")
#library("readxl")
# library("ggplot2")
# library("reshape2")
library("gplots")
# library("plotly")
library("stringr")
# library("caret")
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library('scales')
library('parallel')
library('bestglm')
library('data.table')
library('tidyr')
library('dplyr')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_2")
# env_a <- read_excel('Env_arctic.xlsx')
env_a <- read.table('Env_arctic.txt', header = T)
fraction = commandArgs(trailingOnly = T)[1]
fraction1 = commandArgs(trailingOnly = T)[2]
stations <- paste(env_a$Station, 'SUR', sep='')
headers_G <- read.delim('/env/cns/proj/TaraOcean/scratch/GeneSet-v2/Occurrences/metaG/MATOU-v2.metaG.matrix36.h.long', header=FALSE, stringsAsFactors=FALSE)
headers_G <- append('UID', headers_G)
headers_G <- as.character(headers_G)
headers_T <- read.delim('/env/cns/proj/Taranalyse/TaraEukData/GeneSet-v2/MATOU-v2.metaT.matrix35.h.long', header=FALSE, stringsAsFactors=FALSE)
headers_T <- append('UID', headers_T)
headers_T <- as.character(headers_T)
stations_filter_G <- rep(0, length(headers_G))
for (st in stations){
  stations_filter_G <- stations_filter_G+grepl(st, headers_G)
}
stations_filter_G <- as.logical(stations_filter_G)
if (fraction1!='0'){
  mG <- which(grepl(fraction, headers_G) & stations_filter_G | grepl('UID',headers_G) | grepl(fraction1, headers_G) & stations_filter_G) 
} else{
  mG <- which(grepl(fraction, headers_G) & stations_filter_G | grepl('UID',headers_G))
}
stations_filter_T <- rep(0, length(headers_T))
for (st in stations){
  stations_filter_T <- stations_filter_T+grepl(st, headers_T)
}
stations_filter_T <- as.logical(stations_filter_T)
if (fraction1!='0'){
  mT <- which(grepl(fraction, headers_T) & stations_filter_T | grepl('UID',headers_T) | grepl(fraction1, headers_T) & stations_filter_T) 
} else{
  mT <- which(grepl(fraction, headers_T) & stations_filter_T | grepl('UID',headers_T))
}
# 
st_T <- substr(headers_T[mT], start = 1, stop=6)
st_G <- substr(headers_G[mG], start = 1, stop=6)

if (length(st_T)<length(st_G)){
  mG<-mG[!is.na(match(st_G, st_T))]
  st_T <- substr(headers_T[mT], start = 1, stop=6)
  st_G <- substr(headers_G[mG], start = 1, stop=6)
} else if (length(st_T)<length(st_G)){
  mT<-mT[!is.na(match(st_T, st_G))]
  st_T <- substr(headers_T[mT], start = 1, stop=6)
  st_G <- substr(headers_G[mG], start = 1, stop=6)
}
metaT_ <- fread('/env/cns/proj/Taranalyse/TaraEukData/GeneSet-v2/MATOU-v2.metaT.matrix35', select = mT, header = F)
metaG_ <- fread('/env/cns/proj/TaraOcean/scratch/GeneSet-v2/Occurrences/metaG/mg.mode36.matrix',select = mG, header = F)
# metaT_ <- fread('/env/cns/proj/Taranalyse/TaraEukData/GeneSet-v2/testT.txt', select = mT, header = F)
# metaG_ <- fread('/env/cns/proj/Taranalyse/TaraEukData/GeneSet-v2/testG.txt',select = mG, header = F)

colnames(metaT_) <- st_T
colnames(metaG_) <- st_G

if (fraction1!='0'){
  sel <- !duplicated(colnames(metaG_))
  metaG_ <- as.matrix(metaG_)
  metaG_ <- metaG_[,sel]
  metaG_ <- as.data.frame(metaG_)
}

uids <- unique(metaT_$UID, metaG_$UID)
num <- length(uids)

subseting <- function(i){
  j=which(sequence==i)
  if (i!= sequence[length(sequence)]){
    uids0 <- uids[i:(i+199999)]
  } else {
    uids0 <- uids[i:num]
  }
  tableT <- metaT_[metaT_$UID %in% uids0, ]
  tableG <- metaG_[metaG_$UID %in% uids0, ]
  tableT_ <- tableT %>% gather(key='Station', value='MetaT', -UID)
  tableG_ <- tableG %>% gather(key='Station', value='MetaG', -UID) 
  data <- merge(tableG_, tableT_, by=c('Station', 'UID'), all=T)
  data$MetaG[is.na(data$MetaG)] <- 0
  data$MetaT[is.na(data$MetaT)] <- 0
  saveRDS(data,paste(fraction,"metaTnMetaG_",j,".rds", sep=''))
}
write(num, paste('count_uni',fraction,'.txt', sep=''))

sequence=seq(1,num, 200000)
#no_cores <- detectCores()-1
# Initiate cluster
#cl <- makeCluster(no_cores)
#clusterExport(cl=cl, varlist=c("metaT_","metaG_", 'sequence', 'num', 'fraction', 'uids'))
#cor_list <- parLapply(cl = cl,sequence, fun = subseting)
#stopCluster(cl)
for (i in sequence){
    subseting(i)
}
# metaT_ <- metaT_ %>% gather(key='Station', value='MetaT', -UID)
# metaG_ <- metaG_ %>% gather(key='Station', value='MetaG', -UID)
# 
# #data <- inner_join(metaG_GGZZ, metaT_GGZZ,by=c('Station', 'UID'))
# data <- merge(metaG_, metaT_, by=c('Station', 'UID'), all=T)
# data$MetaG[is.na(data$MetaG)] <- 0
# data$MetaT[is.na(data$MetaT)] <- 0
# saveRDS(data, paste(fraction, 'metaTnMetaG.rds', sep=''))



