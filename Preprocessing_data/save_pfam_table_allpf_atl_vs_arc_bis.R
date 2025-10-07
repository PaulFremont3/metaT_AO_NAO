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
library('dplyr')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")

fraction =commandArgs(trailingOnly = T)[1]# c('GGZZ', 'SSUU', 'QQSS', 'GGZZ', 'QQSS')
n_c = commandArgs(trailingOnly = T)[2]# c(2,2,2,3,3)
taxo = commandArgs(trailingOnly = T)[3] # c('taxo_groups2', 'taxo_MGTv2', 'taxo_groups3')

n_clusts <- as.numeric(n_c)
if (n_clusts==2){
  id <-''
} else{
  id <-paste('_', n_clusts, sep='')
}
if (taxo != "0" ){
  tax_id = taxo
}else{
  tax_id = ''
}
n <- 792
if (taxo != "0" ){
  file <- readRDS(paste('MetaTG_pfams_station_table_allpf_',fraction,'_atl_vs_arc',id,'_',tax_id,'_bis.rds', sep=''))
} else{
  file <- readRDS(paste('MetaTG_pfams_station_table_allpf_',fraction,'_atl_vs_arc',id,'_bis.rds', sep=''))
}
matou_baseline <- readRDS(paste('matou_baseline_allpf_', fraction, '.rds', sep=''))
table_pfam <- NULL
if (fraction=='GGZZ'){
  nstat <- 27
} else if (fraction=='SSUU'){
  nstat <- 31
} else if (fraction=='QQSS'){
  nstat <- 30
} else if (fraction=='KKQQ'){
  nstat <- 28
}
pfams_all <- unique(matou_baseline$Var1)
pfams_all<- as.character(levels(pfams_all))[pfams_all]
pfams_all[is.na(pfams_all)]<-'Unknown'
nrows <-  length(pfams_all)
if (n_clusts %in% c(2,5, 6) & taxo != "0"){
  types <- paste(c('T_com', 'uni_T_com', 'G_com','uni_G_com' , 'uni_T1_com' , 'uni_G1_com',
             'T_atl', 'uni_T_atl', 'G_atl','uni_G_atl', 'uni_T1_atl' , 'uni_G1_atl',
             'T_arc', 'uni_T_arc', 'G_arc','uni_G_arc', 'uni_T1_arc' , 'uni_G1_arc'), tax_id, sep='_')
} else if (n_clusts %in% c(3,4) & taxo != "0"){
  types <- paste(c('T_trans', 'uni_T_trans', 'G_trans','uni_G_trans' ,
             'T_atl', 'uni_T_atl', 'G_atl','uni_G_atl',
             'T_arc', 'uni_T_arc', 'G_arc','uni_G_arc',
             'T_com', 'uni_T_com', 'G_com','uni_G_com'), tax_id, sep='_')
} else if (n_clusts %in% c(2,5,6) & taxo == "0"){
  types <- c('T_com', 'uni_T_com', 'G_com','uni_G_com' ,'uni_T1_com' , 'uni_G1_com',
             'T_atl', 'uni_T_atl', 'G_atl','uni_G_atl', 'uni_T1_atl' , 'uni_G1_atl',
             'T_arc', 'uni_T_arc', 'G_arc','uni_G_arc', 'uni_T1_arc' , 'uni_G1_arc')
} else if (n_clusts %in% c(3,4) & taxo == "0"){
  types <- c('T_trans', 'uni_T_trans', 'G_trans','uni_G_trans' ,
             'T_atl', 'uni_T_atl', 'G_atl','uni_G_atl',
             'T_arc', 'uni_T_arc', 'G_arc','uni_G_arc',
             'T_com', 'uni_T_com', 'G_com','uni_G_com')
}
if (n_clusts==4){
  n_clusts0 <- 3
} else if (n_clusts %in% c(5,6)){
  n_clusts0 <- 2
} else{
  n_clusts0 <- n_clusts
}
su <- readRDS(paste('pfams_sum_T_', fraction,'.rds', sep=''))
suG <- readRDS(paste('pfams_sum_G_', fraction,'.rds', sep=''))
seq0 <- 1:((n_clusts0+1)*4)
if (n_clusts %in% c(2,5,6)){
  seq0 <- 1:((n_clusts0+1)*6)
}
#seq0 <- seq0[seq0%%5 != 0]
for (j in 1:length(types)){
  table_pfam <- matrix(0, ncol=nstat, nrow=nrows)
  list_pfam_taxo <- rep(list(NULL), nrows)
  names(list_pfam_taxo) <- pfams_all
  
  ##concatenating info of whole MATOU (1 to 792)
  for (i in 1:n){
    if (taxo =="0"){
      if (!is.na(file[[i]][[j]])){
        mapping <- match(rownames(file[[i]][[j]]),  pfams_all)
        table_pfam[mapping,] = table_pfam[mapping,]+file[[i]][[j]]
      }
    } else{
      for (pf in names(file[[i]][[j]])){
        list_pfam_taxo[[pf]] <- rbind(list_pfam_taxo[[pf]],file[[i]][[j]][[pf]])
      }
    }
  }
  
  if (taxo =="0"){
    rownames(table_pfam)<-pfams_all
    colnames(table_pfam)<-colnames(file[[1]][[1]])
    if (substr(types[j],1,1) == 'T'){
      table_pfam <- t(table_pfam)/su
      table_pfam <- t(table_pfam)
    }
    if (substr(types[j],1,1) == 'G'){
      table_pfam <- t(table_pfam)/suG
      table_pfam <- t(table_pfam)
    }
    saveRDS(table_pfam, paste('pfams_station_table_',types[j],'_allpf_',fraction,id,'_bis.rds', sep=''))
  } else{
    for (pfam in pfams_all){
      data <- list_pfam_taxo[[pfam]]
      if (!is.null(data)){
        na <- rownames(data)
        data <- as.data.frame(data)
        data$taxo <- na
        data <- data %>% group_by(taxo) %>% summarise_all(funs(sum))
        data <- as.data.frame(data)
        row.names(data) <- data[,1]
        data <- data[,-1]
        if (substr(types[j],1,1) == 'T'){
          data <- t(data)/su
          data <- t(data)
        }
        if (substr(types[j],1,1) == 'G'){
          data <- t(data)/suG
          data <- t(data)
        }
        list_pfam_taxo[[pfam]] <- data
      }
    }
    print(types[j])
    saveRDS(list_pfam_taxo, paste('pfams_station_table_',types[j],'_allpf_',fraction,id,'_bis.rds', sep=''))  
  }
}

