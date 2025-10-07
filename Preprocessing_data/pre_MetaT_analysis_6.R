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

fraction = commandArgs(trailingOnly = T)[1]  # 'GGZZ', 'SSUU', 'QQSS'
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

arctic_stations <- c('158SUR','163SUR',"168SUR","173SUR", "175SUR", "178SUR","180SUR", "188SUR", "189SUR", "191SUR",
                       "193SUR", "194SUR", "196SUR")
atlantic_stations <- c("142SUR", "143SUR", "144SUR","145SUR","146SUR", "147SUR", "148SUR", "149SUR", "150SUR" ,
                         "151SUR" ,"152SUR", '155SUR')
outlier_stations <-c('')


func0 <- function(i){

  class_genes <-function(uni, df0){
    su_T <- sum(df0$MetaT[df0$UID==uni] , na.rm=T)
    if (su_T==0){
      return(4)
    } else{    
      su_T_arctic <- sum(df0$MetaT[df0$UID==uni & df0$Station %in% arctic_stations], na.rm = T)
      su_T_outlier <- sum(df0$MetaT[df0$UID==uni & df0$Station %in% outlier_stations], na.rm = T)
      su_T_atlantic <- sum(df0$MetaT[df0$UID==uni & df0$Station %in% atlantic_stations], na.rm = T)
      if ( su_T_atlantic>0 & su_T_arctic==0 ){
        return(2)
      } else if (su_T_atlantic==0 & su_T_arctic>0){
        return(3)
      } else if (su_T_atlantic>0 & su_T_arctic>0){
        return(1)
      } else{
        return(4)
      }
    }
  }
  expr_uds <- readRDS(paste("subset_metat_",fraction,'/expressed_unilist_',fraction,'_',arg,'.rds', sep=''))
  if (arg != 792){
    uds <- uids[((i-1)*20000+1):((i-1)*20000+20000)]
    uds <- uds[expr_uds[((i-1)*20000+1):((i-1)*20000+20000)]]
  } else if (arg==792 & i==5){
    uds <- uids[((i-1)*20000+1):length(uids)]
    uds <- uds[expr_uds[((i-1)*20000+1):length(uids)]]
  } else if (arg==792 & i!=5){
    uds <- uids[((i-1)*20000+1):((i-1)*20000+20000)]
    uds <- uds[expr_uds[((i-1)*20000+1):((i-1)*20000+20000)]]
  }
  df0 <- df[df$UID %in% uds,]
  uds_split <- sapply(uds, class_genes, df0=df0)
  return(uds_split)
}
#start_time <- Sys.time()
#no_cores <- detectCores()
#cl <- makeCluster(no_cores)
#clusterExport(cl=cl, varlist=c("df", 'uids', 'arg', 'fraction','arctic_stations', 'atlantic_stations', 'outlier_stations'))
if (arg != 792){
  uni_list <- lapply(c(1:10),FUN=func0)
} else{
  uni_list <- lapply(c(1:5),FUN=func0)
}
#stopCluster(cl)
#end_time <- Sys.time()

t <- unlist(uni_list)
saveRDS(t, paste('class_unilist_6_',fraction,'_',arg,'.rds', sep=''))
