#!/bin/env/usr/env Rscript
library("gbm")
library("stringr")
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library('scales')
library('parallel')
library('bestglm')
library('gplots')


frac = commandArgs(trailingOnly = T)
data <- readRDS(paste('../data/metat_analysis_', frac, '.rds', sep=''))
data <- data[data[,111]>4,]
env_a <- read.table('../data/env_arctic_3.txt', header = T)
env_a <- env_a[!(env_a$Station %in% c(142,201,205, 206, 208, 209, 210)),]
variables1 <- colnames(env_a)[4:13]
variables0 <- c('Metagenome', 'Environment', 'Travel_time',variables1)

fdrs<-  c(0.001,0.01, 0.05)
fdrs0<- c('0_001', '0_01', '0_05')
signif_uids_p <- function(i, fdr, sel){
  datasel <- data[sel,]
  sorted_p <- sort(datasel[,i+1], decreasing = T)
  fdr_p <- sorted_p[fdr*length(sorted_p)/2]
  print(fdr_p)
  id_sig_p <- datasel[datasel[,i]>fdr_p,1]
  cor_sig_p <- datasel[datasel[,i]>fdr_p,i]
  pfam_sig_p <- datasel[datasel[,i]>fdr_p,39-2]
  pfams_sig_p <- datasel[datasel[,i]>fdr_p,51-2]
  signif_uids_p <- list(id_sig_p, cor_sig_p, pfam_sig_p, fdr_p, pfams_sig_p)
}

signif_uids_n <- function(i, fdr, sel){
  datasel <- data[sel,]
  sorted_n <- sort(datasel[,i+1], decreasing = F)
  fdr_n <- sorted_n[fdr*length(sorted_n)/2]
  print(fdr_n)
  id_sig_n <- datasel[datasel[,i]<fdr_n,1]
  cor_sig_n <- datasel[datasel[,i]<fdr_n,i]
  pfam_sig_n <- datasel[datasel[,i]<fdr_n,39-2]
  pfams_sig_n <- datasel[datasel[,i]<fdr_n,51-2]
  signif_uids_n <- list(id_sig_n, cor_sig_n, pfam_sig_n, fdr_n, pfams_sig_n)
}

uids_list <- rep(list(), length(fdrs))
uids_list_log <- rep(list(), length(fdrs))
uids_list_pca <- rep(list(), length(fdrs))
uids_list_clr <- rep(list(), length(fdrs))
c=1
sel <- !(data[,105] %in% c('arctic', 'atlantic'))
for (fdr in fdrs){
  significant_uids_p <- lapply(seq(2,28-2,2), FUN = signif_uids_p, fdr=fdr, sel=sel)
  significant_uids_n <- lapply(seq(2,28-2,2), FUN = signif_uids_n, fdr=fdr, sel=sel)
  uids_list[[c]] <- significant_uids_p
  uids_list[[c+1]] <- significant_uids_n
  significant_uids_p_log <- lapply(seq(53-2,79-4,2), FUN = signif_uids_p, fdr=fdr, sel=sel)
  significant_uids_n_log <- lapply(seq(53-2,79-4,2), FUN = signif_uids_n, fdr=fdr, sel=sel)
  uids_list_log[[c]] <- significant_uids_p_log
  uids_list_log[[c+1]] <- significant_uids_n_log
  significant_uids_p_pca <- lapply(seq(81-4,85-4,2), FUN = signif_uids_p, fdr=fdr, sel=sel)
  significant_uids_n_pca <- lapply(seq(81-4,85-4,2), FUN = signif_uids_n, fdr=fdr, sel=sel)
  uids_list_pca[[c]] <- significant_uids_p_pca
  uids_list_pca[[c+1]] <- significant_uids_n_pca
  significant_uids_p_clr <- lapply(c(c(83,4,6),seq(85, 103,2)), FUN = signif_uids_p, fdr=fdr, sel=sel)
  significant_uids_n_clr <- lapply(c(c(83,4,6),seq(85, 103,2)), FUN = signif_uids_n, fdr=fdr, sel=sel)
  uids_list_clr[[c]] <- significant_uids_p_clr
  uids_list_clr[[c+1]] <- significant_uids_n_clr
  c=c+2
}

make_density_correlation_plots <- function(seq1, seq2, data, uids_list, name, title, sel, labx){
  titles <- c('Metagenome', 'Environment','Travel time')
  inds <- c(1:3)
  labels <- c(title, 'Pearson correlation coefficent (MetaT/MetaG)', 'Pearson correlation coefficent (MetaT/MetaG)')
  pdf(paste('correlation_distributions_metaT_vs_metag_env_lagr_',frac,name,'.pdf', sep=''), width = 20, height = 10)
  par(mfrow=c(2,2), mar=c(5.1, 5.1, 4.1, 2.1))
  count = 1
  for (i in seq1){
    d <- density(data[sel,i])
    d_r <- density(data[sel,i+1])
    u=inds[count]
    plot(c(5,5), main=titles[i/2], , xlim=c(-1,1), ylim=c(0,max(d$y, d_r$y)), xlab=labels[count],ylab='Density', col='white', cex.main=2, cex.lab=2 ,cex.axis=2)
    polygon(d, col=alpha("red", 0.5), border=alpha('red', 0.5))
    polygon(d_r, col=alpha("blue", 0.5), border=alpha('blue', 0.5))
    for (k in 1:length(uids_list)){
      if (count>1){
        abline(v=uids_list[[k]][[u]][[4]], col='red',lwd=2)
      }
    }
    n <- length(data[sel,i])
    legend('topleft', legend=paste('n=', n, sep=''), bty='n', cex=2)
    count=count+1
  }
  dev.off()
  
  
  variables <- variables1
  inds <- c(4:13)
  pdf(paste('correlation_distributions_metaT_vs_env_variables_',frac,name,'.pdf', sep=''), width = 15, height = 20)
  par(mfrow=c(4,3), mar=c(6.1, 6.6, 4.1, 2.1) ,mgp=c(5,1,0))
  c<-1
  for (i in seq2){
    datasel <- data[sel,]
    d <- density(datasel[!is.na(datasel[,i]),i])
    d_r <- density(datasel[!is.na(datasel[,i+1]),i+1])
    plot(c(5,5), main=variables[c], xlim=c(-1,1), ylim=c(0,max(d$y, d_r$y)), xlab='Pearson correlation coefficent',ylab='Density', col='white',cex.main=2, cex.lab=2 ,cex.axis=2)
    polygon(d, col=alpha("red", 0.5), border=alpha('red', 0.5))
    polygon(d_r, col=alpha("blue", 0.5), border=alpha('blue', 0.5))
    u=inds[c]
    for (k in c(5,6)){
      abline(v=uids_list[[k]][[u]][[4]], col='red',lwd=2)
    }
    n <- length(datasel[!is.na(datasel[,i]),i])
    legend('topleft', legend=paste('n=', n, sep=''), bty='n', cex=2)
    c <- c+1
  }
  dev.off()
}
#make_density_correlation_plots(seq(2,6,2), seq(8, 28-2,2), data = data, sel=sel,labx='Pearson correlation coefficent\n(MetaT/MetaG)', uids_list = uids_list, name = '', title='Pearson correlation coefficent (MetaT)')
#make_density_correlation_plots(seq(53-2,57-2,2), seq(59-2, 79-4,2), data = data,sel=sel, labx='Pearson correlation coefficent\n(MetaT/MetaG)', uids_list = uids_list_log, name = '_log', title= 'Pearson correlation coefficent (MetaT)')
make_density_correlation_plots(c(83, 4, 6), seq(85, 103,2), data = data[!is.na(data[,83]),],sel=sel,labx='Pearson correlation coefficent\nclr(MetaT)-clr(MetaG)', uids_list = uids_list_clr, name = '_clr', title='Erb rho coefficent (MetaT)')

thresholds <- NULL
for (j in 4:13){
  thresholds <- append(thresholds, uids_list_clr[[5]][[j]][[4]])  
  thresholds <- append(thresholds, uids_list_clr[[6]][[j]][[4]])
}
saveRDS(thresholds, paste('thresholds_signif_env_vs_expresssion_0_05_', frac,'.rds', sep=''))

variables_pca <- c('pca1','pca2', 'pca3')
inds <- c(1:3)
pdf(paste('correlation_distributions_metaT_vs_pca_variables_',frac,'.pdf', sep=''), width = 15, height = 20)
par(mfrow=c(2,2))
c<-1
for (i in seq(81-4,85-4,2)){
  d <- density(data[sel,i])
  d_r <- density(data[sel,i+1])
  plot(c(5,5), main=variables_pca[c], xlim=c(-1,1), ylim=c(0,max(d$y, d_r$y)), xlab='Pearson correlation coefficent (MetaT/MetaG)',ylab='Density', col='white' ,cex.main=2, cex.lab=2 ,cex.axis=2)
  polygon(d, col=alpha("red", 0.5), border=alpha('red', 0.5))
  polygon(d_r, col=alpha("blue", 0.5), border=alpha('blue', 0.5))
  u=inds[c]
  for (k in 1:length(uids_list_pca)){
    abline(v=uids_list_pca[[k]][[u]][[4]], col='darkviolet',lwd=2)
  }
  c <- c+1
}
dev.off()


write_signif <- function(significant_uids_p,significant_uids_n, fdr, name, seq, vars){
  data_to_write <- NULL
  for (i in seq){
    uni_vec_p <- significant_uids_p[i][[1]][[1]]
    cor_vec_p <- significant_uids_p[i][[1]][[2]]
    pfm_vec_p <- significant_uids_p[i][[1]][[3]]
    pfms_vec_p <- significant_uids_p[i][[1]][[5]]
    var_vec_p <- rep(vars[i], length(uni_vec_p))
    
    uni_vec_n <- significant_uids_n[i][[1]][[1]]
    cor_vec_n <- significant_uids_n[i][[1]][[2]]
    pfm_vec_n <- significant_uids_n[i][[1]][[3]]
    pfms_vec_n <- significant_uids_n[i][[1]][[5]]
    var_vec_n <- rep(vars[i], length(uni_vec_n))
    
    data_new <- data.frame(uni=c(uni_vec_p, uni_vec_n),vari=c(var_vec_p, var_vec_n), cor=c(cor_vec_p, cor_vec_n), 
                           pfam=c(as.character(pfm_vec_p), as.character(pfm_vec_n)), pfams=c(as.character(pfms_vec_p), as.character(pfms_vec_n)))
    rownames(data_new) <- NULL
    saveRDS(data_new, file = paste('Significant_uids_all_',frac,'_',fdr,'_',vars[i],name,'.rds', sep=''))
    data_to_write <- rbind(data_to_write, data_new)
  }
  rownames(data_to_write)<- NULL
  saveRDS(data_to_write, file = paste('Significant_uids_all_',frac,'_',fdr,name,'.rds', sep=''))
}

i=1
for (fdr0 in fdrs0){
  #write_signif(uids_list[[2*i-1]], uids_list[[2*i]], fdr0, '', c(1:13), variables0)
  #write_signif(uids_list_log[[2*i-1]], uids_list_log[[2*i]], fdr0, '_log', c(1:13), variables0)
  #write_signif(uids_list_pca[[2*i-1]], uids_list_pca[[2*i]], fdr0, '_pca', c(1:3), variables_pca)
  write_signif(uids_list_clr[[2*i-1]], uids_list_clr[[2*i]], fdr0, '_clr', c(1:13), variables0)
  i = i+1
}


