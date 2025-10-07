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

fraction = commandArgs(trailingOnly = T)[1] # 'GGZZ', 'SSUU', 'QQSS'
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
uids <- unique(df$UID)
good_stat <- good_stat[!(good_stat %in% c('142SUR','201SUR','205SUR', '206SUR', '208SUR', '209SUR', '210SUR'))]
df <- df[df$Station %in% good_stat,]

func <- function(i){
  correlations1 <- function(uni, df0){
    cond <- df0$UID==uni  
    metat1 <- df0$MetaT[cond]
    if (length(metat1[metat1!=0])>4){
      return(T)
    } else{
      return(F)
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
  uds_split <- NULL
  for (u in uds){
    ud <- correlations1(u, df0=df0)
    uds_split <- append(uds_split, ud)
  }
  return(uds_split)
}

no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c("df", 'uids', 'arg'))
if (arg != 792){
  uni_list <- parLapply(cl = cl, c(1:10), fun =func)
} else{
  uni_list <- parLapply(cl = cl, c(1:5), fun =func)
}
stopCluster(cl)


t <- unlist(uni_list)
unilist_all0 <- uids[t==T]
saveRDS(unilist_all0, paste('pre_unilist_',fraction,'_',arg,'.rds', sep=''))
