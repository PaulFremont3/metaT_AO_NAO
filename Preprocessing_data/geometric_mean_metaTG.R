#!bin/usr/bin/env Rscript
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

fraction <- commandArgs(trailingOnly = T)[1]

prod_arg <- function(arg){
  df <- readRDS(paste("subset_metat_",fraction,"/",fraction,"metaTnMetaG_",arg,".rds", sep=''))
  good_stats <- unique(df$Station)
  #if (fraction=='GGZZ'){
  #  good_stat <- stations[stations != '191SUR']
  #} else{
  #  good_stat <- stations
  #}
  #good_stat <- good_stat[!(good_stat %in% c('142SUR','201SUR','205SUR', '206SUR', '208SUR', '209SUR', '210SUR'))]
  #df <- df[df$Station %in% good_stat,] 
  df$MetaT[is.na(df$MetaT)]<-0
  df$MetaG[is.na(df$MetaG)]<-0
  n_stats <- length(good_stats)
  prods_mt <- rep(0, n_stats)
  prods_mg <- rep(0, n_stats)
  n_Ts <- rep(0, n_stats)
  n_Gs <- rep(0, n_stats)
  j=1
  for (u in good_stats){
    prodT <- sum(log(df$MetaT[df$Station==u & df$MetaT != 0]), na.rm=T)
    nT <- sum(df$Station==u & df$MetaT != 0, na.rm=T)
    n_Ts[j] =nT
    prods_mt[j] <- prodT
    prodG <- sum(log(df$MetaG[df$Station==u & df$MetaG != 0]), na.rm=T)
    nG <- sum(df$Station==u & df$MetaG != 0, na.rm=T)
    n_Gs[j] =nG
    prods_mg[j] <- prodG
    j=j+1
  }
  return(list(prods_mt, n_Ts, prods_mg, n_Gs))
}

n=792
no_cores <- detectCores()
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c("fraction"))
prods_args <- parLapply(cl = cl, c(1:n), fun = prod_arg)
stopCluster(cl)
n_stats <- length(prods_args[[1]][[1]])
prod_finalT <- rep(0, n_stats)
nT <- rep(0, n_stats)
prod_finalG <- rep(0, n_stats)
nG <- rep(0, n_stats)
for (pr in prods_args){
  pT=pr[[1]]
  n_T=pr[[2]]
  prod_finalT <- prod_finalT+pT
  nT <- nT+n_T
  pG=pr[[3]]
  n_G=pr[[4]]
  prod_finalG <- prod_finalG+pG
  nG <- nG+n_G
}
prod_finalT <- exp((1/nT)*prod_finalT)
prod_finalG <- exp((1/nG)*prod_finalG)
saveRDS(prod_finalT, paste('geom_mean_metaT_',fraction,'.rds', sep=''))
saveRDS(prod_finalG, paste('geom_mean_metaG_',fraction,'.rds', sep=''))
  

