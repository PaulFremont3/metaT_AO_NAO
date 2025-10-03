# library("rgdal")                                                                                                      
# library("raster")
library('stringr')
library("ggplot2")
#library('readxl')
library('gridExtra')
library('tidyr')
library('vegan')
#library('DT')
library('dplyr')
library('RColorBrewer')
library('viridis')
#setwd("~/Groups_metaT/Groups_plots/")
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4/Groups_metaT/Groups_plots")
source('functions_scatter_plots.R')

type <- commandArgs(trailingOnly=TRUE)[1]
#Tmin1 <- readr::read_tsv("minAijji_tarrive_min_surface_1000.csv") %>% 
Tmin1 <- read.table('minAijji_tarrive_min_surface_1000.csv') 
  #rename(Station_1 = X1) %>% 
Station_1 <- rownames(Tmin1)
Tmin1 <- cbind(Station_1, Tmin1)
new_names <- lapply(colnames(Tmin1)[2:dim(Tmin1)[2]], function(x){strsplit(x, split = 'X')[[1]][2]})
colnames(Tmin1)[2:dim(Tmin1)[2]] <- new_names
Tmin1 <- Tmin1 %>% 
  dplyr::select(Station_1, ends_with("_SUR")) %>% 
  filter(grepl("_SUR", Station_1)) %>%     
  gather(-Station_1, key = "Station_2", value = "Tmin") %>% 
  mutate(type = "Tmin")

Tmin1$Station_1 <- as.character(sapply(Tmin1$Station_1, 
                                       FUN = function(x){paste(strsplit(x, split = '_')[[1]],
                                                               collapse = '')}))
Tmin1$Station_2 <- as.character(sapply(Tmin1$Station_2, 
                                       FUN = function(x){paste(strsplit(x, split = '_')[[1]],
                                                               collapse = '')}))

simka <- readRDS('simkaTG_long_arctic.rds')
env_arctic <- read.table('../env_arctic_3.txt', header = T)
env_arctic$Station <- paste(env_arctic$Station, 'SUR', sep='')
stations <- env_arctic$Station
good_stat <- stations[stations>'142SUR' & stations <'201SUR']
latis <- env_arctic$Latitude
lonis <- env_arctic$Longitude

fractions <- c('GGZZ')
bc_unigenes <- list()
for (f in fractions){
  if (f!='KKQQ_MMQQ'){
    bc_file <- readRDS(paste('MetaTG_bray-curtis_unigenes_',f,'.rds', sep=''))
  } else {
    bc_file <- readRDS('MetaTG_bray-curtis_unigenes_KKQQ.rds')
  }
  names(bc_file) <- c('T', 'G')
  bc_unigenes[[f]] <- bc_file
}

my_round = function(x, n=2) {
  max(abs(round(x, n)), abs(signif(x, 2))) * sign(x)
}

env_dist <- readRDS('../environmental_distances_arctic_3.rds')
colnames(env_dist) <- paste(colnames(env_dist), 'SUR', sep='')
rownames(env_dist) <- colnames(env_dist)


colos <- list('SSUU'='darkorange','GGZZ'= 'dodgerblue2','QQSS'= 'darkviolet', 'KKQQ_MMQQ'='darkgreen')
basins <- c('all', 'green', 'blue')
basins_an <- c( 'all', 'Arctic-Arctic', 'Atlantic-Atlantic') 
data_slopes <- NULL
for (frac in fractions){

  test <- extract_data(type, frac)
  condit=test$tmin <1/3
  test$sk[condit]=NA
  test$sk_g[condit]=NA
  test$bc[condit]=NA
  test$bc_g[condit]=NA
  test$tmin[condit]=NA
  print(frac)
  print(cor(test$sk, test$bc, use = 'pairwise.complete.obs'))
  print(cor(test$sk_g, test$bc_g, use = 'pairwise.complete.obs'))
  print(cor(test$sk, test$sk_g, use = 'pairwise.complete.obs'))
  print(cor(test$bc, test$bc_g, use = 'pairwise.complete.obs'))
  
  aty <- c(0.01,0.1, 1, 10)
  aty0 <- c(seq(2, 9, by= 1) %o% 10^(-3:1))
  labels <- sapply(aty,function(i)
    as.expression(bquote(i))
  )
  n_bins <-50
  if (frac=='GGZZ'){
    varis <- c('ed', 'sk', 'bc')
  } else{
    varis <- c('sk', 'bc')
  }
  for (vari in varis){
    for (j in 1:3){
      b <- basins_an[j]
      colob <- basins[j]
      if (colob=='all'){
        test0 <- test
      } else{
        test0 <- test[test$col==colob,]
      }
      #trans_stats <-c('155SUR', '158SUR', '163SUR')
      test0$col[test0$col=='blue']<- 'darkorange'
      test0$col[test0$col=='green']<- 'darkblue'
      test0$col[test0$col=='red']<- 'forestgreen'
      #test0$col[test0$st1 %in% trans_stats]<-'red'
      #test0$col[test0$st2 %in% trans_stats]<-'red'
      
      #print(sqrt(dim(test0)[1]))
      cum_cors <- cum_cor(test0, var=vari)
      cum_cors_ed <- cum_cor_ed(test0, var=vari)
      condu <- which(cum_cors$mt <0.01)
      cond0 <- which(cum_cors$mt >0.01)
      lines_1 <- find_lines(condu, 1)
      lines_2 <- find_lines(cond0, 2)
      
      condu_e <- which(cum_cors_ed$mt <0.01)
      cond0_e <- which(cum_cors_ed$mt >0.01)
      lines_3 <- find_lines(condu_e, 1)
      lines_4 <- find_lines(cond0_e, 2)
      
      mx <- max(test0[[vari]], na.rm=T)
      if (vari=='ed'){
        mi <- 0
      } else{
        mi <- min(test0[[vari]][test0[[vari]]!=0], na.rm=T)
      }
      oldmin <- 0
      cum_cors$cor2 <- (((cum_cors$cor - oldmin) * (mx- mi)) / (1 - oldmin)) + mi
      cum_cors_ed$cor2 <- (((cum_cors_ed$cor - oldmin) * (mx- mi)) / (1 - oldmin)) + mi
      
      to_plot <- test0 %>% arrange(st1, st2) %>% filter(st1 < st2)
      if (vari!='ed'){
        pdf(paste('lagrangian_',type,'_',vari,'_',frac,'_',b,'.pdf', sep=''), width=9)
      } else{
        pdf(paste('lagrangian_',type,'_',vari,'_',b,'.pdf', sep=''), width=9)
      }
      if (vari=='ed'){
        yl = 'Environmental distance'
        colo = 'orchid1'
      } else if (vari=='sk'){
        yl = 'Metatranscriptomic dissimilarity'
        colo=colos[[frac]]
      } else if (vari=='bc'){
        yl = 'Metatranscriptomic dissimilarity (unigenes)'
        colo=colos[[frac]]
      }
      par(mar=c(5.1, 5.1, 4.1, 4.6))
      
      if (vari %in% c('sk', 'bc')){
        print(vari)
        print(frac)
        print(b)
        sls <- NULL
        for (t in c(100, 8, 4)){
          y=1-to_plot[[vari]]
          x=to_plot$tmin
          y=y[x<t]
          x=x[x<t]
          lin_mod <- lm(log(y)~x)
          lin_mod1 <- lm(y~x)
          slope <- round(-log(2)/lin_mod$coefficients[2],1)
          pv <- summary(lin_mod)
          pva <- pv$coefficients[8]
	  #print(as.numeric(slope))
          #print(pva)
          #print(summary(lin_mod))
          sls <- append(sls, as.numeric(slope))
          sls <- append(sls, my_round(pva))
        }
        data_slopes <- rbind(data_slopes, c(vari, frac, b, sls))
      }
      plot(to_plot$tmin, to_plot[[vari]], log='x', col=to_plot$col, xlab='Tmin (years)',
           ylab=yl, pch=19, xaxt='n', cex.lab=2, cex.axis=2, cex=1.3)
      for (cond in lines_1){
        points(cum_cors$thresh[cond], cum_cors$cor2[cond], col=colo, type='l',lwd=5)
      }
      for (cond in lines_2){
        points(cum_cors$thresh[cond], cum_cors$cor2[cond], col=colo, type='l',lwd=5, lty=2)
      }
      axis(1,at=aty,labels=aty,  cex.lab=2, cex.axis=2)
      axis(1,at=aty0,  cex.lab=2, cex.axis=2, tck=-0.005, labels=rep(NA, length(aty0)))
      ats <- seq(mi, mx, (mx-mi)/5)
      labs <- round(seq(oldmin, 1, (1-oldmin)/5), 2)
      axis(4, at=ats, labels = labs,  cex.lab=2, cex.axis=2)
      mtext("Spearman correlation", side = 4, line=3, cex=2)
      if (b=="all"){
        legend('topleft', legend = c('Atlantic-Atlantic', 'Arctic-Arctic', 'Atlantic-Arctic'), 
               col = c('darkorange', 'darkblue', 'forestgreen'), pch = 19, cex=1.3)
      } else{
        legend('topleft', legend = b, 
               col = colob, pch = 19, cex=1.3)
      }
      legend('bottomright', legend='Spearman correlation', col=colo, lty = 1, lwd=5, cex=1.3)
      plot(to_plot$tmin, to_plot[[vari]], log='xy', col=to_plot$col, xlab='Tmin (years)',
           ylab=yl, pch=19, xaxt='n', cex.lab=2, cex.axis=2, cex=1.3)
      axis(1,at=aty,labels=aty,  cex.lab=2, cex.axis=2)
      axis(1,at=aty0,  cex.lab=2, cex.axis=2, tck=-0.005, labels=rep(NA, length(aty0)))
      for (cond in lines_1){
        points(cum_cors$thresh[cond], cum_cors$cor2[cond], col=colo, type='l',lwd=5)
      }
      for (cond in lines_2){
        points(cum_cors$thresh[cond], cum_cors$cor2[cond], col=colo, type='l',lwd=5, lty=2)
      }
      axis(1,at=aty,labels=aty,  cex.lab=2, cex.axis=2)
      axis(1,at=aty0,  cex.lab=2, cex.axis=2, tck=-0.005, labels=rep(NA, length(aty0)))
      axis(4, at=ats, labels = labs,  cex.lab=2, cex.axis=2)
      mtext("Spearman correlation", side = 4, line=3, cex=2)
      if (b=="all"){
        legend('topleft', legend = c('Atlantic-Atlantic', 'Arctic-Arctic', 'Atlantic-Arctic'), 
               col = c('darkorange', 'darkblue', 'forestgreen'), pch = 19, cex=1.3)
      } else{
        legend('topleft', legend = b, 
               col = colob, pch = 19, cex=1.3)
      }
      legend('bottomright', legend='Spearman correlation', col=colo, lty = 1, lwd=5, cex=1.75)
      
      if (vari %in% c('sk', 'bc')){
        vari0 <- 'ed'
        
        colos_p <- round(100*(to_plot[[vari0]]- min(test[[vari0]], na.rm=T))/(max(test[[vari0]], na.rm=T)- min(test[[vari0]], na.rm=T)))
        palette_col <- viridis(100, alpha=0.7)
        plot(to_plot$tmin, to_plot[[vari]], log='x', col=palette_col[colos_p], xlab='Tmin (years)',
             ylab=yl, pch=19, xaxt='n', cex.lab=2, cex.axis=2, cex=1.5)
        for (cond in lines_3){
          points(cum_cors_ed$thresh[cond], cum_cors_ed$cor2[cond], col=colo, type='l',lwd=5)
        }
        for (cond in lines_4){
          points(cum_cors_ed$thresh[cond], cum_cors_ed$cor2[cond], col=colo, type='l',lwd=5, lty=2)
        }
        axis(1,at=aty,labels=aty,  cex.lab=2, cex.axis=2)
        axis(1,at=aty0,  cex.lab=2, cex.axis=2, tck=-0.005, labels=rep(NA, length(aty0)))
        ats <- seq(mi, mx, (mx-mi)/5)
        labs <- round(seq(oldmin, 1, (1-oldmin)/5), 2)
        axis(4, at=ats, labels = labs,  cex.lab=2, cex.axis=2)
        mtext("Spearman correlation", side = 4, line=3, cex=2)
        legend('bottomright', legend='Spearman correlation', col=colo, lty = 1, lwd=5, cex=1.75)
        pnt <- cbind(x =c(0,2,2,0), y =c(0,50,50,0))
        plot(0,0, col='white', xlim=c(0,17), ylim=c(0, 60), axes=FALSE, frame.plot=F, xlab = '', ylab='')
        mi <- min(test[[vari0]], na.rm=T)
        ma <- max(test[[vari0]], na.rm=T)
#        SDMTools::legend.gradient(pnt, palette_col, limits=c(mi,ma), title ='Environmental distance' , cex=1) 
      }
      dev.off()
    }
  }
}
write.table(data_slopes, paste('decay_rates_',type,'.txt', sep=''))
saveRDS(data_slopes, paste('decay_rates_',type,'.rds', sep=''))
