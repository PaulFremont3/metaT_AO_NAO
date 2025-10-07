#!/bin/env/usr/env Rscript
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library('scales')
library('parallel')
library('bestglm')
library('reshape2')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")

fraction = commandArgs(trailingOnly = T)[1]# 'GGZZ'
n <- 792
extract_pfam_metaT <- function(arg){
  expr_uds <- readRDS(paste("subset_metat_",fraction,'/expressed_unilist_',fraction,'_',arg,'.rds', sep=''))
  df <- readRDS(paste("subset_metat_",fraction,"/",fraction,"metaTnMetaG_",arg,".rds", sep=''))
  uids <- unique(df$UID)
  uids <- uids[expr_uds]
  df <- df[df$UID %in% uids,]
  df$MetaT[is.na(df$MetaT)]<-0
  df$MetaG[is.na(df$MetaG)]<-0
  if (fraction=='GGZZ'){
    df <- df[df$Station != '191SUR',]
  }
  matou <- readRDS(paste("subset_metat_",fraction,"/",'matouUK', fraction, '_', arg,'.rds', sep=''))
  colnames(matou)<- c('geneID','pfam','Evalue')
  matou <-as.data.frame(matou)
  matou <- matou[matou$geneID %in% uids,]
  counts <- as.data.frame(table(matou$geneID))
  
  v<-merge(df, matou, by.x='UID', by.y='geneID', all=T)

  v$uni <- 1
  v$uni[v$MetaT==0]<-0
  v$uniG <- 1
  v$uniG[v$MetaG==0]<-0
  
  v1 <- v
  v1$id <- paste(v$UID, v$Station)
  v1 <- v1[!duplicated(v1$id), ]  
  
  test <- aggregate(v$MetaT , by = list(v$pfam,v$Station), FUN = sum)
  test1 <- reshape2::acast(test, test$Group.1~test$Group.2)
  
  test2 <- aggregate(v$uni , by = list(v$pfam,v$Station), FUN = sum)
  test3 <- reshape2::acast(test2, test2$Group.1~test2$Group.2)
  
  test4 <- aggregate(v$MetaG , by = list(v$pfam,v$Station), FUN = sum)
  test5 <- reshape2::acast(test4, test4$Group.1~test4$Group.2)
  
  test6 <- aggregate(v$uniG , by = list(v$pfam,v$Station), FUN = sum)
  test7 <- reshape2::acast(test6, test6$Group.1~test6$Group.2)
  
  test8 <- aggregate(v1$uni , by = list(v1$pfam,v1$Station), FUN = sum)
  test9 <- reshape2::acast(test8, test8$Group.1~test8$Group.2)
  
  test10 <- aggregate(v1$uniG , by = list(v1$pfam,v1$Station), FUN = sum)
  test11 <- reshape2::acast(test8, test8$Group.1~test8$Group.2)
  return(list(test1, test3, test5, test7, test9, test11))
}

no_cores <- detectCores()
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c('fraction'))
whole_data <- parLapply(cl=cl, c(1:n), fun=extract_pfam_metaT)
stopCluster(cl)
saveRDS(whole_data, paste('MetaTG_pfams_station_table_allpf_',fraction,'_bis.rds', sep=''))
