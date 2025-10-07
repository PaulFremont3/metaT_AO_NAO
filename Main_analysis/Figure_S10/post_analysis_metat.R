#!/bin/env/usr/env Rscript
library("gbm")
#library("dismo")
#library("readxl")
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
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")

frac = commandArgs(trailingOnly = T)
data <- readRDS(paste('metat_analysis_', frac, '.rds', sep=''))
data <- data[data[,111]>4,]
env_a <- read.table('env_arctic_v2.txt', header = T)
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
make_density_correlation_plots(seq(2,6,2), seq(8, 28-2,2), data = data, sel=sel,labx='Pearson correlation coefficent\n(MetaT/MetaG)', uids_list = uids_list, name = '', title='Pearson correlation coefficent (MetaT)')
make_density_correlation_plots(seq(53-2,57-2,2), seq(59-2, 79-4,2), data = data,sel=sel, labx='Pearson correlation coefficent\n(MetaT/MetaG)', uids_list = uids_list_log, name = '_log', title= 'Pearson correlation coefficent (MetaT)')
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
  write_signif(uids_list[[2*i-1]], uids_list[[2*i]], fdr0, '', c(1:13), variables0)
  write_signif(uids_list_log[[2*i-1]], uids_list_log[[2*i]], fdr0, '_log', c(1:13), variables0)
  write_signif(uids_list_pca[[2*i-1]], uids_list_pca[[2*i]], fdr0, '_pca', c(1:3), variables_pca)
  write_signif(uids_list_clr[[2*i-1]], uids_list_clr[[2*i]], fdr0, '_clr', c(1:13), variables0)
  i = i+1
}


overexpressed_func <- function(data,i,j,k,l, m,name){
  p_vals <- as.numeric(levels(data[data[,i]!='arctic' & data[,i]!='atlantic',i])[data[data[,i]!='arctic' & data[,i]!='atlantic',i]])
  p_vals_r <- data[data[,i]!='arctic' & data[,i]!='atlantic',k]
  
  locs <- data[data[,i]!='arctic' & data[,i]!='atlantic',j]
  uni_ids <- data[data[,i]!='arctic' & data[,i]!='atlantic',1]
  pfams <- data[data[,i]!='arctic' & data[,i]!='atlantic',l]
  pfams_all =  data[data[,i]!='arctic' & data[,i]!='atlantic',m]
  # ratio_t <- data[data[,37]!='arctic' & data[,37]!='atlantic',40]
  # ratio_g <- data[data[,37]!='arctic' & data[,37]!='atlantic',41]
  p_vals_fdr <- p.adjust(p_vals, method='fdr')
  cond_signif <- p_vals_fdr<0.05
  uids_arctic <- uni_ids[cond_signif & locs=='arctic']
  pfams_arctic <- pfams_all[cond_signif & locs=='arctic']
  pfams_all_arctic <- pfams_all[cond_signif & locs=='arctic']
  uids_atlantic <- uni_ids[cond_signif & locs=='atlantic']
  pfams_atlantic <- pfams_all[cond_signif & locs=='atlantic']
  pfams_all_atlantic <- pfams_all[cond_signif & locs=='atlantic']

  arctic_overexpressed <- data.frame(uni=uids_arctic, pfam=pfams_arctic, pfams=pfams_all_arctic)
  atlantic_overexpressed <- data.frame(uni=uids_atlantic, pfam=pfams_atlantic, pfams=pfams_all_atlantic)
  saveRDS(arctic_overexpressed,file=paste('Significant_uids_all_',frac,'_fdr_arctic',name,'.rds', sep=''))
  saveRDS(atlantic_overexpressed,file=paste('Significant_uids_all_',frac,'_fdr_atlantic',name,'.rds', sep=''))  

  p_vals_sorted <- sort(p_vals)
  p_vals_sorted_r <- sort(p_vals_r)
  locs_sorted <- locs[order(p_vals, decreasing = F)]
  uni_ids_sorted <- uni_ids[order(p_vals, decreasing = F)]
  pfams_sorted <- pfams[order(p_vals, decreasing = F)]
  pfams_all_sorted <- pfams_all[order(p_vals, decreasing = F)]
  bins <- exp(seq(log(10^-5), log(10^0), length.out = 1000))
  cum_pval <- NULL
  cum_pval_r <- NULL
  for (b in bins){
    cum_pval <- append(cum_pval, sum(p_vals_sorted<b)/length(p_vals_sorted))
    cum_pval_r <- append(cum_pval_r, sum(p_vals_sorted_r<b)/length(p_vals_sorted_r))
  }
  fdrs <- c(0.001, 0.01, 0.05, 0.1, 0.15)
  fdrs0 <- c('0_001', '0_01', '0_05', '0_1', '0_15')
  i=1
  for (fdr in fdrs){
    fdr_coff <- p_vals_sorted_r[fdr*length(p_vals_sorted_r)]
    cond_signif <- p_vals_sorted < fdr_coff
    
    pdf(paste('cumulative_distribution_pval_arctic_atlantic_overexpressed_', frac,'_',fdrs0[i], name,'.pdf', sep=''), width = 10, height = 10)
    plot(bins, cum_pval, log='x', col='darkgreen', type='l', lwd=3, xaxt = 'n' ,cex.lab = 1.5, cex.axis = 1.5, 
         xlab='p-value', ylab='Cumulative distribution', yaxs='i', xaxs='i', ylim=c(0,1))
    points(bins, cum_pval_r, col='blue',type='l', lwd=3)
    abline(v=fdr_coff, col='red', lwd=3)
    axis(side = 1, at = c(2:10 %o% 10^(-5:1)), labels = F, tck=-0.005)
    axis(side = 1, at = c(10 %o% 10^(-5:1)), labels = c(10 %o% 10^(-5:1)), cex.axis = 1.5)
    legend("topleft", legend= c('Data', 'Random'),fill=c('darkgreen', 'blue'), cex=1.5, bty= 'n')
    dev.off()
    
    uids_arctic <- uni_ids_sorted[cond_signif & locs_sorted=='arctic']
    pfams_arctic <- pfams_sorted[cond_signif & locs_sorted=='arctic']
    pfams_all_arctic <- pfams_all_sorted[cond_signif & locs_sorted=='arctic']
    uids_atlantic <- uni_ids_sorted[cond_signif & locs_sorted=='atlantic']
    pfams_atlantic <- pfams_sorted[cond_signif & locs_sorted=='atlantic']
    pfams_all_atlantic <- pfams_all_sorted[cond_signif & locs_sorted=='atlantic']
    
    arctic_overexpressed <- data.frame(uni=uids_arctic, pfam=pfams_arctic, pfams=pfams_all_arctic)
    atlantic_overexpressed <- data.frame(uni=uids_atlantic, pfam=pfams_atlantic, pfams=pfams_all_atlantic)
    saveRDS(arctic_overexpressed,file=paste('Significant_uids_all_',frac,'_',fdrs0[i],'_arctic',name,'.rds', sep=''))
    saveRDS(atlantic_overexpressed,file=paste('Significant_uids_all_',frac,'_',fdrs0[i],'_atlantic',name,'.rds', sep=''))
    i=i+1
  }
}
overexpressed_func(data, 37-2, 38-2, 46-2, 39-2,51-2, name='' )
overexpressed_func(data, 105, 106, 109, 39-2,51-2, name='_clr' )
# ## Holm correction
# correction_holm <- length(p_vals)- c(1:length(p_vals))+1
# test_holm <- p_vals_sorted*correction_holm ## nothing is significant
# ## FDR adjustment
# correction_fdr <- c(1:length(p_vals))/length(p_vals)
# test_fdr <- p_vals_sorted/correction_fdr
# test_fdr[1:which.min(test_fdr)]<- min(test_fdr)
# methods <- c("holm", "hochberg",  "bonferroni", "BH", "BY","fdr")
# for (m in methods){
#   pval <- p.adjust(p_vals_sorted, m)
#   print(min(pval))
# }
# fdr_test <- p.adjust(p_vals_sorted, 'fdr')
# # uids_arctic <- uni_ids_sorted[p_vals_sorted<p_vals_sorted[length(adjustment_fdr)+1] & locs_sorted=='arctic']
# # pfams_arctic <- pfams_sorted[p_vals_sorted<p_vals_sorted[length(adjustment_fdr)+1] & locs_sorted=='arctic']
# uids_arctic <- uni_ids_sorted[fdr_test<0.05 & locs_sorted=='arctic']
# pfams_arctic <- pfams_sorted[fdr_test<0.05 & locs_sorted=='arctic']
# # uids_atlantic <- uni_ids_sorted[locs_sorted=='atlantic' & p_vals_sorted<p_vals_sorted[length(adjustment_fdr)+1]]
# # pfams_atlantic <- pfams_sorted[locs_sorted=='atlantic' & p_vals_sorted<p_vals_sorted[length(adjustment_fdr)+1]]
# uids_atlantic <- uni_ids_sorted[fdr_test<0.05 & locs_sorted=='atlantic']
# pfams_atlantic <- pfams_sorted[fdr_test<0.05 & locs_sorted=='atlantic']

arctic_only <- data.frame(uni=data[data[,37-2]=='arctic',1], pfam=data[data[,37-2]=='arctic',39-2], pfams=data[data[,37-2]=='arctic',51-2])
atlantic_only <- data.frame(uni=data[data[,37-2]=='atlantic',1], pfam=data[data[,37-2]=='atlantic',39-2], pfams=data[data[,37-2]=='atlantic',51-2])
saveRDS(arctic_only,file=paste('Significant_uids_all_',frac,'_arctic_only.rds', sep=''))
saveRDS(atlantic_only,file=paste('Significant_uids_all_',frac,'_atlantic_only.rds', sep=''))

pf_list = c('PF05971', 'PF01180', 'PF05351', 'PF05035', 'PF01346', 'PF03104', 'PF11999',
            'PF00249', 'PF03259', 'PF02209', 'PF00241', 'PF00626', 'PF01119','PF10551' )
arctic = readRDS(paste('Significant_uids_all_',frac,'_0_001_arctic.rds', sep=''))
atl = readRDS(paste('Significant_uids_all_',frac,'_0_001_atlantic.rds', sep=''))
#atl = atl[atl$cor<0,]

pdf(paste('plot_ratios_t_g_overlap_no_overlap_', frac,'.pdf', sep=''), height = 13, width = 13)
#cond= as.numeric(data[,1] %in% unique(atl[,1], na.rm=T))+1
par(mfrow = c(2,2))
#cond= as.numeric(data[,1] %in% unique(atl[,1], na.rm=T))+1
cond = as.numeric(!is.na(data[,39-2]))+1
cond1 = !is.na(data[,39-2])
#cond1=data[,1] %in% unique(arctic[,1], na.rm=T)
cols = c(alpha('darkgreen', 0.5), alpha('red', 0.5))
col_vec= cols[cond]
u = data[,43-2]
v = data[,42-2]
condu = u!=0 & !is.na(u) & v!=0 & !is.na(v)
u0 = u[condu]
v0 = v[condu]
u1 = u[condu & !is.na(data[,39-2])]
#maxi = max(abs(log10(u0)), abs(log10(v0)))
maxi=3
plot(log10(data[condu & cond1,43-2]), log10(data[condu& cond1,42-2]), xlab='log(MetaG_Arctic/MetaG_Atlantic)', 
     ylab='log(MetaT_Arctic/MetaT_Atlantic)', col=densCols(log10(data[condu & cond1,43-2]), 
                                                           log10(data[condu& cond1,42-2]),
                                                           colramp=colorRampPalette(RColorBrewer::brewer.pal(9, 'Reds')[-(1:4)])), 
     xlim=c(-maxi,maxi), 
     ylim=c(-maxi,maxi), pch=19)
abline(h=0, col='blue',lwd=2 )
abline(v=0, col='lightblue',lwd=2)
abline(0,1, col='purple', lwd=2)
legend('topleft', col = c('red'), 
       legend = c(paste('Known function',length(log10(data[condu & cond1,43-2])),sep=' ')), pch=1)
title('Overlap metaT/G')

plot(log10(data[condu & !cond1,43-2]), log10(data[condu& !cond1,42-2]), xlab='log(MetaG_Arctic/MetaG_Atlantic)', 
     ylab='log(MetaT_Arctic/MetaT_Atlantic)', col=densCols(log10(data[condu & !cond1,43-2]), 
                                                           log10(data[condu& !cond1,42-2]),
                                                           colramp=colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[-(1:4)])), 
     xlim=c(-maxi,maxi), 
     ylim=c(-maxi,maxi), pch=19)
abline(h=0, col='blue',lwd=2 )
abline(v=0, col='lightblue',lwd=2)
abline(0,1, col='purple', lwd=2)
legend('topleft', col = c('darkgreen'), 
       legend = c(paste('Unknown function',length(log10(data[condu & !cond1,43-2])),sep=' ')), pch=1)
title('Overlap metaT/G')


u = data[,41-2]
v = data[,40-2]
condu = u!=0 & !is.na(u) & v!=0 & !is.na(v)
cond1 = !is.na(data[,39-2])
#cond1=data[,1] %in% unique(arctic[,1], na.rm=T)
u0 = u[condu]
v0 = v[condu]
u1 = u[condu & !is.na(data[,39-2])]
u_x = data[,43-2]
v_x = data[,42-2]
condu_x = !is.na(u_x) & !is.na(v_x)
#maxi = max(abs(log10(u0)), abs(log10(v0)))
plot(log10(data[condu & cond1 & !condu_x,41-2]), log10(data[condu& cond1 & !condu_x,40-2]), xlab='log(MetaG_Arctic/MetaG_Atlantic)', 
     ylab='log(MetaT_Arctic/MetaT_Atlantic)', col=densCols(log10(data[condu & cond1 & !condu_x,41-2]), 
                                                           log10(data[condu& cond1 & !condu_x,40-2]),
                                                           colramp=colorRampPalette(RColorBrewer::brewer.pal(9, 'Reds')[-(1:4)])), 
     xlim=c(-maxi,maxi), 
     ylim=c(-maxi,maxi), pch=19)
abline(h=0, col='blue',lwd=2 )
abline(v=0, col='lightblue',lwd=2)
abline(0,1, col='purple', lwd=2)
legend('topleft', col = c('red'), 
       legend = c(paste('Known function',length(log10(data[condu & cond1 & !condu_x,41-2])),sep=' ')), pch=1)
title('No Overlap metaT/G')

plot(log10(data[condu & !cond1 & !condu_x,41-2]), log10(data[condu& !cond1 & !condu_x,40-2]), xlab='log(MetaG_Arctic/MetaG_Atlantic)', 
     ylab='log(MetaT_Arctic/MetaT_Atlantic)', col=densCols(log10(data[condu & !cond1 & !condu_x,41-2]), 
                                                           log10(data[condu& !cond1 & !condu_x,40-2]),
                                                           colramp=colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[-(1:4)])), 
     xlim=c(-maxi,maxi), 
     ylim=c(-maxi,maxi), pch=19)
abline(h=0, col='blue',lwd=2 )
abline(v=0, col='lightblue',lwd=2)
abline(0,1, col='purple', lwd=2)
legend('topleft', col = c('darkgreen'), 
       legend = c(paste('Unknown function',length(log10(data[condu & !cond1,41-2])),sep=' ')), pch=1)
title('No Overlap metaT/G')
dev.off()

pdf(paste('plot_ratios_t_g_overlap_no_overlap_', frac,'_rand.pdf', sep=''), height = 13, width = 13)
#cond= as.numeric(data[,1] %in% unique(atl[,1], na.rm=T))+1
par(mfrow = c(2,2))
#cond= as.numeric(data[,1] %in% unique(atl[,1], na.rm=T))+1
cond = as.numeric(!is.na(data[,39-2]))+1
cond1 = !is.na(data[,39-2])
#cond1=data[,1] %in% unique(arctic[,1], na.rm=T)
cols = c(alpha('darkgreen', 0.5), alpha('red', 0.5))
col_vec= cols[cond]
u = data[,50-2]
v = data[,49-2]
condu = u!=0 & !is.na(u) & v!=0 & !is.na(v)
u0 = u[condu]
v0 = v[condu]
u1 = u[condu & !is.na(data[,39-2])]
#maxi = max(abs(log10(u0)), abs(log10(v0)))
maxi=3
plot(log10(data[condu & cond1,50-2]), log10(data[condu& cond1,49-2]), xlab='log(MetaG_Arctic/MetaG_Atlantic)', 
     ylab='log(MetaT_Arctic/MetaT_Atlantic)', col=densCols(log10(data[condu & cond1,50-2]), 
                                                           log10(data[condu& cond1,49-2]),
                                                           colramp=colorRampPalette(RColorBrewer::brewer.pal(9, 'Reds')[-(1:4)])), 
     xlim=c(-maxi,maxi), 
     ylim=c(-maxi,maxi), pch=19)
abline(h=0, col='blue',lwd=2 )
abline(v=0, col='lightblue',lwd=2)
abline(0,1, col='purple', lwd=2)
legend('topleft', col = c('red'), 
       legend = c(paste('Known function',length(log10(data[condu & cond1,50-2])),sep=' ')), pch=1)
title('Overlap metaT/G')

plot(log10(data[condu & !cond1,50-2]), log10(data[condu& !cond1,49-2]), xlab='log(MetaG_Arctic/MetaG_Atlantic)', 
     ylab='log(MetaT_Arctic/MetaT_Atlantic)', col=densCols(log10(data[condu & !cond1,50-2]), 
                                                           log10(data[condu& !cond1,49-2]),
                                                           colramp=colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[-(1:4)])), 
     xlim=c(-maxi,maxi), 
     ylim=c(-maxi,maxi), pch=19)
abline(h=0, col='blue',lwd=2 )
abline(v=0, col='lightblue',lwd=2)
abline(0,1, col='purple', lwd=2)
legend('topleft', col = c('darkgreen'), 
       legend = c(paste('Unknown function',length(log10(data[condu & !cond1,43])),sep=' ')), pch=1)
title('Overlap metaT/G')


u = data[,48-2]
v = data[,47-2]
condu = u!=0 & !is.na(u) & v!=0 & !is.na(v)
cond1 = !is.na(data[,39-2])
#cond1=data[,1] %in% unique(arctic[,1], na.rm=T)
u0 = u[condu]
v0 = v[condu]
u1 = u[condu & !is.na(data[,39-2])]
u_x = data[,50-2]
v_x = data[,49-2]
condu_x = !is.na(u_x) & !is.na(v_x)
#maxi = max(abs(log10(u0)), abs(log10(v0)))
plot(log10(data[condu & cond1 & !condu_x,48-2]), log10(data[condu& cond1 & !condu_x,47-2]), xlab='log(MetaG_Arctic/MetaG_Atlantic)', 
     ylab='log(MetaT_Arctic/MetaT_Atlantic)', col=densCols(log10(data[condu & cond1 & !condu_x,48-2]), 
                                                           log10(data[condu& cond1 & !condu_x,47-2]),
                                                           colramp=colorRampPalette(RColorBrewer::brewer.pal(9, 'Reds')[-(1:4)])), 
     xlim=c(-maxi,maxi), 
     ylim=c(-maxi,maxi), pch=19)
abline(h=0, col='blue',lwd=2 )
abline(v=0, col='lightblue',lwd=2)
abline(0,1, col='purple', lwd=2)
legend('topleft', col = c('red'), 
       legend = c(paste('Known function',length(log10(data[condu & cond1 & !condu_x,48-2])),sep=' ')), pch=1)
title('No Overlap metaT/G')

plot(log10(data[condu & !cond1 & !condu_x,48-2]), log10(data[condu& !cond1 & !condu_x,47-2]), xlab='log(MetaG_Arctic/MetaG_Atlantic)', 
     ylab='log(MetaT_Arctic/MetaT_Atlantic)', col=densCols(log10(data[condu & !cond1 & !condu_x,48-2]), 
                                                           log10(data[condu& !cond1 & !condu_x,47-2]),
                                                           colramp=colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[-(1:4)])), 
     xlim=c(-maxi,maxi), 
     ylim=c(-maxi,maxi), pch=19)
abline(h=0, col='blue',lwd=2 )
abline(v=0, col='lightblue',lwd=2)
abline(0,1, col='purple', lwd=2)
legend('topleft', col = c('darkgreen'), 
       legend = c(paste('Unknown function',length(log10(data[condu & !cond1,41-2])),sep=' ')), pch=1)
title('No Overlap metaT/G')
dev.off()

pf_list = c('PF05971', 'PF01180', 'PF05351', 'PF05035', 'PF01346', 'PF03104', 'PF11999',
            'PF00249', 'PF03259', 'PF02209', 'PF00241', 'PF00626', 'PF01119','PF10551' )
Correlation_matrix <- NULL
uni_list <- NULL
cols_uni <- NULL
pfam_list <- NULL
pfams_all_list <- list()
col_vec = colorRampPalette(brewer.pal(11,"Spectral"))(11)
col_vec <- append(col_vec, c('red', 'darkgreen', 'darkblue'))
c=1
for (pf in pf_list){
  unis <- data[data[,39-2]==pf & !is.na(data[,39-2]),1]
  pfams_all <- data[data[,39-2]==pf & !is.na(data[,39-2]),51-2]
  print(pfams_all)
  pfams_all_list[[pf]] <- pfams_all
  uni_list <- append(uni_list, unis)
  cols_uni <- append(cols_uni, rep(col_vec[c], length(unis)))
  pfam_list <- append(pfam_list, rep(pf, length(unis)))
  cors <- data[data[,39-2]==pf & !is.na(data[,39-2]), seq(2,28-2,2)]
  Correlation_matrix <- rbind(Correlation_matrix, cors)
  print(pf)
  print(length(unis))
  c=c+1
}
Correlation_matrix <- as.matrix(Correlation_matrix)
pdf(paste('pfams_of_interest_', frac,'.pdf', sep=''), width=10, height =10)
heatmap.2(Correlation_matrix, trace="none",Rowv = NA, 
          dendrogram = "none", keysize=1,margins=c(10,22), labCol=variables0,
          labRow=pfam_list, col= redgreen(15),  symkey = F, colRow =  cols_uni)
hist(Correlation_matrix[,5], breaks=100,main='Correlation with temperature', xlab='Pearson correlation coefficient', xlim=c(-1,1))
heatmap.2(Correlation_matrix, trace="none",
          dendrogram = "none", keysize=1,margins=c(10,22), labCol=variables0,
          labRow=pfam_list, col= redgreen(15),  symkey = F, colRow =  cols_uni)
for (pf in pf_list){
  if (sum(pfam_list==pf, na.rm = T)>=2){
    heatmap.2(Correlation_matrix[pfam_list==pf,], trace="none",
            dendrogram = "none", keysize=1,margins=c(10,22), labCol=variables0, main = pf,
            labRow=uni_list[pfam_list==pf], col= redgreen(15),  symkey = F, colRow =  cols_uni[pfam_list==pf])
    hist(Correlation_matrix[pfam_list == pf,5], breaks=100, main='Correlation with temperature', xlab='Pearson correlation coefficient', 
         col=cols_uni[pfam_list==pf][1], border=scales::alpha('white', 0), xlim=c(-1,1))
  }
}
dev.off()
