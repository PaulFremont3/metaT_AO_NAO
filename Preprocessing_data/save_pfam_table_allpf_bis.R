#!/bin/env/usr/env Rscript
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('RColorBrewer')
library('ncdf4')
library('scales')
library('parallel')
library('dplyr')

fraction = commandArgs(trailingOnly = T)[1]
taxo = commandArgs(trailingOnly = T)[2]
if (taxo != "0" ){
  tax_id = taxo
}else{
  tax_id = ''
}

n <- 792
if (taxo==0){
  file <- readRDS(paste('MetaTG_pfams_station_table_allpf_',fraction,'_bis.rds', sep=''))
} else{
  file <- readRDS(paste('MetaTG_pfams_station_table_allpf_',fraction,'_',tax_id,'_bis.rds', sep=''))
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
table_pfam <- matrix(0, ncol=nstat, nrow=nrows)
table_pfam_uni <- matrix(0, ncol=nstat, nrow=nrows)
table_pfamG <- matrix(0, ncol=nstat, nrow=nrows)
table_pfam_uniG <- matrix(0, ncol=nstat, nrow=nrows)
table_pfam_uni1 <- matrix(0, ncol=nstat, nrow=nrows)
table_pfam_uniG1 <- matrix(0, ncol=nstat, nrow=nrows)
list_pfam_taxo <- rep(list(NULL), nrows)
list_pfam_taxo_uni <- rep(list(NULL), nrows)
list_pfam_taxoG <- rep(list(NULL), nrows)
list_pfam_taxo_uniG <- rep(list(NULL), nrows)
list_pfam_taxo_uni1 <- rep(list(NULL), nrows)
list_pfam_taxo_uniG1 <- rep(list(NULL), nrows)
names(list_pfam_taxo) <- pfams_all
for (i in 1:n){
  if (taxo==0){
    mapping <- match(rownames(file[[i]][[1]]),  pfams_all)
    table_pfam[mapping,] = table_pfam[mapping,]+file[[i]][[1]]
    table_pfam_uni[mapping,] = table_pfam_uni[mapping,]+file[[i]][[2]]
    table_pfamG[mapping,] = table_pfamG[mapping,]+file[[i]][[3]]
    table_pfam_uniG[mapping,] = table_pfam_uniG[mapping,]+file[[i]][[4]]
    if (!is.null(file[[i]][[5]])){
      mapping_bis <- match(rownames(file[[i]][[5]]),  pfams_all)
      mapping_ter <- match(rownames(file[[i]][[6]]),  pfams_all)
      table_pfam_uni1[mapping_bis,] = table_pfam_uni1[mapping_bis,]+file[[i]][[5]]
      table_pfam_uniG1[mapping_ter,] = table_pfam_uniG1[mapping_ter,]+file[[i]][[6]]
    }
  } else{
    for (pf in names(file[[i]][[1]])){
      list_pfam_taxo[[pf]] <- rbind(list_pfam_taxo[[pf]],file[[i]][[1]][[pf]])
      list_pfam_taxo_uni[[pf]] <- rbind(list_pfam_taxo_uni[[pf]],file[[i]][[2]][[pf]])
      list_pfam_taxoG[[pf]] <- rbind(list_pfam_taxoG[[pf]],file[[i]][[3]][[pf]])
      list_pfam_taxo_uniG[[pf]] <- rbind(list_pfam_taxo_uniG[[pf]],file[[i]][[4]][[pf]])
      if (!is.null(file[[i]][[5]][[pf]])){
        list_pfam_taxo_uni1[[pf]] <- rbind(list_pfam_taxo_uni1[[pf]],file[[i]][[5]][[pf]])
        list_pfam_taxo_uniG1[[pf]] <- rbind(list_pfam_taxo_uniG1[[pf]],file[[i]][[6]][[pf]])
      }
    }
  }
}
if (taxo==0){
  rownames(table_pfam)<-pfams_all
  colnames(table_pfam)<-colnames(file[[1]][[1]])
  rownames(table_pfam_uni)<-pfams_all
  colnames(table_pfam_uni)<-colnames(file[[1]][[1]])
  rownames(table_pfamG)<-pfams_all
  colnames(table_pfamG)<-colnames(file[[1]][[1]])
  rownames(table_pfam_uniG)<-pfams_all
  colnames(table_pfam_uniG)<-colnames(file[[1]][[1]])
  rownames(table_pfam_uni1)<-pfams_all
  colnames(table_pfam_uni1)<-colnames(file[[1]][[1]])
  rownames(table_pfam_uniG1)<-pfams_all
  colnames(table_pfam_uniG1)<-colnames(file[[1]][[1]])
  su <- apply(table_pfam, 2, sum)
  saveRDS(su, paste('pfams_sum_T_', fraction,'.rds', sep=''))
  table_pfam <- t(table_pfam)/su
  table_pfam <- t(table_pfam)
  suG <- apply(table_pfamG, 2, sum)
  saveRDS(su, paste('pfams_sum_G_', fraction,'.rds', sep=''))
  table_pfamG <- t(table_pfamG)/suG
  table_pfamG <- t(table_pfamG)
  saveRDS(table_pfam, paste('pfams_station_table_T_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(table_pfam_uni, paste('pfams_station_table_uni_T_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(table_pfamG, paste('pfams_station_table_G_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(table_pfam_uniG, paste('pfams_station_table_uni_G_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(table_pfam_uni1, paste('pfams_station_table_uni_T1_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(table_pfam_uniG1, paste('pfams_station_table_uni_G1_allpf_',fraction,'_bis.rds', sep=''))
} else{
  su <- readRDS(paste('pfams_sum_T_', fraction,'.rds', sep=''))
  suG <- readRDS(paste('pfams_sum_G_', fraction,'.rds', sep=''))
  taxi_sum <- function(data0, sup){
    data_o <- data0
    for (pfam in pfams_all){
      data <- data0[[pfam]]
      if (!is.null(data)){
        na <- rownames(data)
        data <- as.data.frame(data)

        data$taxo <- na
        data <- data %>% group_by(taxo) %>% summarise_all(funs(sum))
        data <- as.data.frame(data)
        row.names(data) <- data[,1]
        data <- data[,-1]
        if (sup != 0){
          data <- t(data)/sup
          data <- t(data)
        }
        data_o[[pfam]] <- data
      }
    }
    return(data_o)
  }
  pfam_taxo <- taxi_sum(list_pfam_taxo, su)
  pfam_taxo_uni <- taxi_sum(list_pfam_taxo_uni, 0)
  pfam_taxoG <- taxi_sum(list_pfam_taxoG, suG)
  pfam_taxo_uniG <- taxi_sum(list_pfam_taxo_uniG, 0)
  pfam_taxo_uni1 <- taxi_sum(list_pfam_taxo_uni1, 0)
  pfam_taxo_uniG1 <- taxi_sum(list_pfam_taxo_uniG1, 0)
  saveRDS(pfam_taxo, paste('pfams_station_table_T_',tax_id,'_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(pfam_taxo_uni, paste('pfams_station_table_uni_T_',tax_id,'_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(pfam_taxoG, paste('pfams_station_table_G_',tax_id,'_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(pfam_taxo_uniG, paste('pfams_station_table_uni_G_',tax_id,'_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(pfam_taxo_uni1, paste('pfams_station_table_uni_T1_',tax_id,'_allpf_',fraction,'_bis.rds', sep=''))
  saveRDS(pfam_taxo_uniG1, paste('pfams_station_table_uni_G1_',tax_id,'_allpf_',fraction,'_bis.rds', sep=''))
}



