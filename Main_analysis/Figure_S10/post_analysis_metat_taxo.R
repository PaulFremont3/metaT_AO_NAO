#!/bin/env/usr/env Rscript
library("gbm")
library("dismo")
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('RColorBrewer')
library('ncdf4')
library('scales')
library('parallel')
library('bestglm')


frac = commandArgs(trailingOnly = T)[1]
taxo = commandArgs(trailingOnly = T)[2]
data <- readRDS(paste('../data/metat_analysis_', frac, '.rds', sep=''))
data <- data[data[,111]>4,]
sel <- !(data[,105] %in% c('arctic', 'atlantic'))
data <- data[sel,]
env_a <- read.table('../data/env_arctic_3.txt', header = T)
env_a <- env_a[!(env_a$Station %in% c(142,201,205, 206, 208, 209, 210)),]
variables1 <- colnames(env_a)[4:13]
variables0 <- c('Metagenome', 'Environment', 'Travel_time',variables1)
if (taxo=='taxo'){
  data_uni_tax <- readRDS('taxID_uni.rds')
} else {
  data_uni_tax <- readRDS(paste('data_uni_',taxo,'.rds', sep=''))
  colnames(data_uni_tax)<-c('taxName', 'geneID')
  data_uni_tax<- data_uni_tax[!duplicated(data_uni_tax$geneID),]
} 

thresholds <- readRDS(paste('thresholds_signif_env_vs_expresssion_0_05_',frac,'.rds', sep=''))

data$taxo <- data_uni_tax$taxName[match(data$V1, data_uni_tax$geneID)]
data$taxo[is.na(data$taxo)]<-'unknown'

count_taxo <- as.data.frame(table(data$taxo))
count_taxo<-count_taxo[order(count_taxo$Freq, decreasing = T),]

make_density_correlation_plots <- function(seq1, seq2, data, name, colo, title, title1,tax){
  titles <- c('Metagenome', 'Environment','Travel time')
  inds <- c(1:3)
  labels <- c(title,'Pearson correlation coefficent (MetaT/MetaG)', 'Pearson correlation coefficent (MetaT/MetaG)')
  pdf(paste('correlation_distributions_metaT_vs_metag_env_lagr_',frac,name,'.pdf', sep=''), width = 20, height = 10)
  par(mfrow=c(2,2),mar=c(5.1, 5.1, 4.1, 2.1) )
  count = 1
  data0 <- data[data$taxo==tax,]
  for (i in seq1){
    d <- density(data0[!is.na(data0[,i]),i])
    d_r <- density(data[!is.na(data[,i+1]),i+1])
    u=inds[count]
    plot(c(5,5), main=titles[i/2], , xlim=c(-1,1), ylim=c(0,max(d$y, d_r$y)), xlab=labels[count],ylab='Density', col='white',cex.main=2, cex.lab=2 ,cex.axis=2)
    polygon(d, col=alpha(colo, 0.5), border=alpha('black', 0.5))
    polygon(d_r, col=alpha("blue", 0.5), border=alpha('black', 0.5))

     n <- length(data0[!is.na(data0[,i]),i])
     legend('topleft', paste('n=',n, sep=''), bty='n', cex=2)
     count=count+1
  }
  dev.off()
  
  
  variables <- variables1#c('SSD','T', 'Salinity', 'Chla', 'O2', 'NO3', 'NO2', 'NH4', 'Fe', 'Phos', 'Si')
  inds <- c(4:13)
  pdf(paste('correlation_distributions_metaT_vs_env_variables_',frac,name,'.pdf', sep=''), width = 15, height = 20)
  par(mfrow=c(4,3), mar=c(6.1, 6.6, 4.1, 2.1) ,mgp=c(5,1,0))
  c<-1
  for (i in seq2){
    d <- density(data0[!is.na(data0[,i]),i])
    d_r <- density(data[!is.na(data[,i+1]),i+1])
    plot(c(5,5), main=variables[c], xlim=c(-1,1), ylim=c(0,max(d$y, d_r$y)), xlab=title1,ylab='Density', col='white',cex.main=2, cex.lab=2 ,cex.axis=2)
    polygon(d, col=alpha(colo, 0.5), border=alpha('black', 0.5))
    polygon(d_r, col=alpha("blue", 0.5), border=alpha('black', 0.5))
    
    #u=inds[c]
    for (k in 1:2){
       abline(v=thresholds[(c-1)*2+k], col='red',lwd=2)
    }
    n=length(data0[!is.na(data0[,i]),i])
    legend('topleft', paste('n=',n, sep=''), bty='n', cex=2)
    c <- c+1
  }
  dev.off()
}
if (taxo=='MGT'){
  colors_taxon <- readRDS('color_table_MGT.rds')
  M=5
} else if (taxo=='groups2'){
  colors_taxon <- readRDS(paste('color_table_',taxo,'.rds',sep=''))
  M=19
} else if (taxo=='MGT-v2'){
  colors_taxon <- readRDS(paste('color_table_',taxo,'.rds',sep=''))
  M=20
} else if (taxo=='groups3'){
  colors_taxon <- readRDS(paste('color_table_',taxo,'.rds',sep=''))
  M=29
}

for (taxon in count_taxo$Var1[1:M]){
  #data0 <- data[data$taxo==taxon,]
  make_density_correlation_plots(seq1 = c(83, 55, 57), seq(59-2, 79-4,2), data = data, 
                                 name = paste('_',taxo, '_',taxon,'_log', sep=''),
                                 colo=as.character(colors_taxon$col[match(taxon, colors_taxon$taxon)]), 
                                 title='Erb rho coefficent', tax=taxon, title1='Pearson correlation coefficent\nlog((MetaT+1)/(MetaG+1))')
  make_density_correlation_plots(seq1 = c(83, 55, 57), seq(85, 103,2), data = data, 
                                 name = paste('_',taxo, '_',taxon,'_clr', sep=''),
                                 colo=as.character(colors_taxon$col[match(taxon, colors_taxon$taxon)]), 
                                 title='Erb rho coefficent', tax=taxon, title1='Pearson correlation coefficent')
  
}








