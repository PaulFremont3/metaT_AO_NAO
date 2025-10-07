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

fraction = commandArgs(trailingOnly = T)[1]# 'GGZZ', 'SSUU', 'QQSS'
n <- 792
extract_values <- function(arg){
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
  diff <- function(x, dt){
    apply(dt, 2,FUN = function(j){sum(abs(j-x))})
  }
  su <- function(x, dt){
    apply(dt, 2,FUN = function(j){sum(abs(j+x))})
  }
  
  t <- reshape2::acast(df,df$UID~df$Station, value.var = 'MetaT')
  g <- reshape2::acast(df,df$UID~df$Station, value.var = 'MetaG' )
  
  diffs_t <- apply(t, MARGIN = 2,FUN =  diff, dt=t)
  sums_t <- apply(t, MARGIN = 2,FUN =  su, dt=t)
  
  diffs_g <- apply(g, MARGIN = 2,FUN =  diff, dt=g)
  sums_g <- apply(g, MARGIN = 2,FUN =  su, dt=g)
  
  return(list(diffs_t,sums_t, diffs_g, sums_g))
}

no_cores <- detectCores()
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c('fraction'))
whole_data <- parLapply(cl=cl, c(1:n), fun=extract_values)
stopCluster(cl)
bc_t_d <- matrix(0, ncol=dim(whole_data[[1]][[1]])[1], nrow=dim(whole_data[[1]][[1]])[1])
bc_t_s <- matrix(0, ncol=dim(whole_data[[1]][[1]])[1], nrow=dim(whole_data[[1]][[1]])[1])
bc_g_d <- matrix(0, ncol=dim(whole_data[[1]][[1]])[1], nrow=dim(whole_data[[1]][[1]])[1])
bc_g_s <- matrix(0, ncol=dim(whole_data[[1]][[1]])[1], nrow=dim(whole_data[[1]][[1]])[1])
for (j in 1:length(whole_data)){
  bc_t_d <- bc_t_d+whole_data[[j]][[1]]
  bc_t_s <- bc_t_s+whole_data[[j]][[2]]
  bc_g_d <- bc_g_d+whole_data[[j]][[3]]
  bc_g_s <- bc_g_s+whole_data[[j]][[4]]
}
bc_t <- bc_t_d/bc_t_s
bc_g <- bc_g_d/bc_g_s
saveRDS(list(bc_t, bc_g), paste('MetaTG_bray-curtis_unigenes_',fraction,'.rds', sep=''))
