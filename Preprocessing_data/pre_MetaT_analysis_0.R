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
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")

fraction = commandArgs(trailingOnly = T)[1]  #'KKQQ' ,'GGZZ', 'SSUU', 'QQSS'
arg = commandArgs(trailingOnly = T)[2] 
df <- readRDS(paste("subset_metat_",fraction,"/",fraction,"metaTnMetaG_",arg,".rds", sep=''))
# uids <- readRDS(paste('subset_metat_QQSS/unilist_QQSS_',arg,'.rds', sep=''))
# df <- df[df$UID %in% uids,]

df$MetaT[is.na(df$MetaT)]<-0
df$MetaG[is.na(df$MetaG)]<-0

stations <- unique(df$Station)
if (fraction=='GGZZ'){
  good_stat <- stations[stations != '191SUR']
} else{
  good_stat <- stations
}
good_stat <- good_stat[!(good_stat %in% c('142SUR','201SUR','205SUR', '206SUR', '208SUR', '209SUR', '210SUR'))]
df <- df[df$Station %in% good_stat,]
uids <- unique(df$UID)

func0 <- function(i){
  
  class_genes <-function(uni, df0){
    su_T <- sum(df0$MetaT[df0$UID==uni] , na.rm=T)
    if (su_T==0){
      return(F)
    } else{  
      return(T)
    }
  }
  if (arg != 792){
    uds <- uids[((i-1)*20000+1):((i-1)*20000+20000)]
  } else if (arg==792 & i==5){
    uds <- uids[((i-1)*20000+1):length(uids)]
  } else if (arg==792 & i!=5){
    uds <- uids[((i-1)*20000+1):((i-1)*20000+20000)]
  }
  df0 <- df[df$UID %in% uds,]
  uds_split <- sapply(uds, class_genes, df0=df0)
  return(uds_split)
}
start_time <- Sys.time()
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c("df", 'uids', 'arg'))
if (arg != 792){
  uni_list <- parLapply(cl = cl, c(1:10), fun =func0)
} else{
  uni_list <- parLapply(cl = cl, c(1:5), fun =func0)
}
stopCluster(cl)
end_time <- Sys.time()

t <- unlist(uni_list)
saveRDS(t, paste('expressed_unilist_',fraction,'_',arg,'.rds', sep=''))
