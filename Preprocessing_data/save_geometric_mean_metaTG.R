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
prod_finalT<-readRDS(paste('geom_mean_metaT_',fraction,'.rds', sep=''))
prod_finalG <- readRDS(paste('geom_mean_metaG_',fraction,'.rds', sep=''))
prod_arg <- function(arg){
  df <- readRDS(paste("subset_metat_",fraction,"/",fraction,"metaTnMetaG_",arg,".rds", sep=''))
  good_stats <- unique(df$Station)
  df$MetaT[df$MetaT==0]<-NA
  df$MetaG[df$MetaG==0]<-NA
  j=1
  for (u in good_stats){
    df$MetaT[df$Station==u & !is.na(df$MetaT)] <- log(df$MetaT[df$Station==u & !is.na(df$MetaT)]/prod_finalT[j])
    df$MetaG[df$Station==u & !is.na(df$MetaG)] <- log(df$MetaG[df$Station==u & !is.na(df$MetaG)]/prod_finalG[j])
    j=j+1
  }
  saveRDS(df,paste(fraction,"metaTnMetaG_",arg,"_clr.rds", sep='') )
  return(NULL)
}

no_cores <- detectCores()
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c("fraction", "prod_finalT", "prod_finalG"))
prods_args <- parLapply(cl = cl, c(1:792), fun = prod_arg)
stopCluster(cl)

  

