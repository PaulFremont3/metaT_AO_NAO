#!/bin/env/usr/env Rscript
#library("readxl")
# library("ggplot2")
# library("reshape2")
library("gplots")
# library("plotly")
library("stringr")
# library("caret")
library('mapproj')
library('mapplots')
#library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library('scales')
library('parallel')
#library('bestglm')
library('reshape2')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")

fraction = commandArgs(trailingOnly = T)[1]# 'GGZZ', 'SSUU', 'QQSS'
taxo = commandArgs(trailingOnly = T)[2] ## '0', 'MGT', 'MGT-v2', 'class', 'groups', 'groups2'
if (taxo==0){
  taxo <- ''
}
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
  
  if (taxo == ''){
    taxi <- readRDS(paste("subset_metat_",fraction,"/TaxID_",fraction, '_',arg,'.rds',sep=''))
  } else{
    taxi <- readRDS(paste("subset_metat_",fraction,"/TaxID_",taxo, '_',fraction, '_',arg,'.rds',sep=''))
  }
  
  v<-merge(df, matou, by.x='UID', by.y='geneID', all=T)
  v$pfam <-as.character(levels(v$pfam))[v$pfam]
  #v$pfam[is.na(v$pfam)]<-'Unknown'
  v$count <- counts$Freq[match(v$UID,counts$Var1)]
  v$count[is.na(v$count)]<-1
  v$MetaT <- v$MetaT/v$count
  v$MetaG <- v$MetaG/v$count
  v$uni <- 1
  v$uni[v$MetaT==0]<-0
  v$uniG <- 1
  v$uniG[v$MetaG==0]<-0
  v$tax <- taxi[match(v$UID, uids)]
  
  t1 <- aggregate(v$MetaT , by = list(v$tax,v$Station, v$pfam), FUN = sum)
  t2 <- aggregate(v$uni , by = list(v$tax,v$Station, v$pfam), FUN = sum)
  t3 <- aggregate(v$MetaG , by = list(v$tax,v$Station, v$pfam), FUN = sum)
  t4 <- aggregate(v$uniG , by = list(v$tax,v$Station, v$pfam), FUN = sum)
  
  extract_tax_pfam <- function(pfam, t){
    u <- t[t$Group.3==pfam,]
    t2 <- reshape2::acast(u, u$Group.1~u$Group.2)
    return(t2)
  }
  
  data_pf1 <- lapply(unique(v$pfam),FUN=extract_tax_pfam, t=t1)
  names(data_pf1) <- unique(v$pfam)
  data_pf2 <- lapply(unique(v$pfam),FUN=extract_tax_pfam, t=t2)
  names(data_pf2) <- unique(v$pfam)
  data_pf3 <- lapply(unique(v$pfam),FUN=extract_tax_pfam, t=t3)
  names(data_pf3) <- unique(v$pfam)
  data_pf4<- lapply(unique(v$pfam),FUN=extract_tax_pfam, t=t4)
  names(data_pf4) <- unique(v$pfam)
  
  return(list(data_pf1, data_pf2, data_pf3, data_pf4))
}

no_cores <- 30
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c('fraction', 'taxo'))
whole_data <- parLapply(cl=cl, c(1:n), fun=extract_pfam_metaT)
stopCluster(cl)
if (taxo==''){
  saveRDS(whole_data, paste('MetaTG_pfams_station_table_allpf_',fraction,'_taxo.rds', sep=''))
} else{
  saveRDS(whole_data, paste('MetaTG_pfams_station_table_allpf_',fraction,'_taxo_',taxo,'.rds', sep=''))
}
