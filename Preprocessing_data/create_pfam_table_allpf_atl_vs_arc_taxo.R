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
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library('scales')
library('parallel')
library('bestglm')
library('reshape2')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")

fraction = commandArgs(trailingOnly = T)[1]# 'GGZZ', 'SSUU', 'QQSS', 'KKQQ'
n_clusts =  commandArgs(trailingOnly = T)[2]
n_clusts = as.numeric(n_clusts)
n=792
taxo = commandArgs(trailingOnly = T)[3]
if (taxo==0){
  taxo <- ''
}
extract_pfam_metaT <- function(arg){
  expr_uds <- readRDS(paste("subset_metat_",fraction,'/expressed_unilist_',fraction,'_',arg,'.rds', sep=''))
  df <- readRDS(paste("subset_metat_",fraction,"/",fraction,"metaTnMetaG_",arg,".rds", sep=''))
  uids <- unique(df$UID)
  uids <- uids[expr_uds]
  df <- df[df$UID %in% uids,]
  if (fraction=='GGZZ'){
    df <- df[df$Station != '191SUR',]
  }
  df$MetaT[is.na(df$MetaT)]<-0
  df$MetaG[is.na(df$MetaG)]<-0
  matou <- readRDS(paste("subset_metat_",fraction,"/",'matouUK', fraction, '_', arg,'.rds', sep=''))
  colnames(matou)<- c('geneID','pfam','Evalue')
  matou <-as.data.frame(matou)
  matou$geneID <- as.numeric(levels(matou$geneID ))[matou$geneID ]
  matou <- matou[matou$geneID %in% uids,]
  
  counts <- as.data.frame(table(matou$geneID))
  counts$Var1 <- as.numeric(levels(counts$Var1))[counts$Var1]
  counts <- counts[order(counts$Var1, decreasing = F),]
  #counts <- counts[counts$Var1 %in% uids,]
  rownames(counts)<-NULL
  
  class_pf <- readRDS(paste("subset_metat_",fraction,"/class_unilist_",n_clusts,"_",fraction, '_',arg,'.rds',sep=''))
  if (taxo == ''){
    taxi <- readRDS(paste("subset_metat_",fraction,"/TaxID_",fraction, '_',arg,'.rds',sep=''))
  } else{
    taxi <- readRDS(paste("subset_metat_",fraction,"/TaxID_",taxo, '_',fraction, '_',arg,'.rds',sep=''))
  }
  if (n_clusts ==4){
    n_clusts0=3
  } else if (n_clusts %in% c(5,6)){
    n_clusts0=2
  }  else{
    n_clusts0=n_clusts
  }
  output_list <- rep(list(NA), (n_clusts0+1)*4)
  for (cl in 1:(n_clusts0+1)){
    ok_uid <- counts$Var1[class_pf==cl]
    tax <- taxi[class_pf==cl]
    if (length(ok_uid)>0){
      df0 <- df[df$UID %in% ok_uid,]
      matou0 <- matou[matou$geneID %in% ok_uid,]
      v<-merge(df0, matou0, by.x='UID', by.y='geneID', all=T)
      v$pfam <-as.character(levels(v$pfam))[v$pfam]
      v$pfam[is.na(v$pfam)]<-'Unknown'
      v$count <- counts$Freq[match(v$UID,counts$Var1)]
      v$count[is.na(v$count)]<-1
      v$MetaT <- v$MetaT/v$count
      v$MetaG <- v$MetaG/v$count
      v$uni <- 1
      v$uni[v$MetaT==0]<-0
      v$uniG <- 1
      v$uniG[v$MetaG==0]<-0
      v$tax <- tax[match(v$UID, ok_uid)]
      
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
      
      output_list[[(cl-1)*4+1]]=data_pf1
      output_list[[(cl-1)*4+2]]=data_pf2
      output_list[[(cl-1)*4+3]]=data_pf3
      output_list[[(cl-1)*4+4]]=data_pf4
    }
  }
  
  return(output_list)
}

no_cores <- detectCores()
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c('fraction', 'n_clusts', 'taxo'))
whole_data <- parLapply(cl=cl, c(1:n), fun=extract_pfam_metaT)
stopCluster(cl)
if (n_clusts==2){
  id <- ''
} else{
  id <- paste('_', n_clusts, sep='')
}
if (taxo==''){
  saveRDS(whole_data, paste('MetaTG_pfams_station_table_allpf_',fraction,'_atl_vs_arc',id,'_taxo.rds', sep=''))
} else{
  saveRDS(whole_data, paste('MetaTG_pfams_station_table_allpf_',fraction,'_atl_vs_arc',id,'_taxo_',taxo,'.rds', sep=''))
}
