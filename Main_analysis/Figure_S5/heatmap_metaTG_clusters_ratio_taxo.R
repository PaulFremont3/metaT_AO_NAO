#!/bin/env/usr/env Rscript
library('gplots')
library('viridis')
library('FactoMineR')
library('RColorBrewer')
library('viridis')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")
source('vioplot.R')

arctic_stations <- c("158SUR","163SUR","168SUR" ,"173SUR", "175SUR", "178SUR","180SUR", "188SUR", "189SUR","191SUR",
                     "193SUR", "194SUR", "196SUR")
atlantic_stations <- c("143SUR" ,"144SUR", "145SUR","146SUR" , "147SUR", "148SUR","149SUR" ,"150SUR" ,
                       "151SUR", "152SUR", "155SUR")
ty <- 'T'
ty0 <- 'G'
tax_id <- 'taxo_groups3'
frac <-'GGZZ'
datat <- readRDS(paste('Meta', ty,'_',tax_id,'_', frac,'.rds', sep=''))
datag <- readRDS(paste('Meta', ty0,'_',tax_id,'_', frac,'.rds', sep=''))
datag <- datag[match(rownames(datat), rownames(datag)),]
clr_func <- function(x){
  log(x/exp(mean(log(x[x>0]))))
}
datat_clr <- apply(datat, 2, clr_func)
datag_clr <- apply(datag, 2, clr_func)

wilc_test <- function(x){
  na <- rownames(x)
  x <- as.numeric(x)
  t <-wilcox.test(x[colnames(datat) %in% arctic_stations], x[colnames(datat) %in% atlantic_stations])$p.value
  loc <- ifelse(median(x[colnames(datat) %in% arctic_stations])>median(x[colnames(datat) %in% atlantic_stations]), yes = 'arc', no='atl')
  return(c(na,t, loc))
}

cols <- readRDS('color_table_groups3.rds')
cols$col <- as.character(levels(cols$col))[cols$col]

cond <- colnames(datat_clr) %in% c(arctic_stations, atlantic_stations, '155SUR')
pca_clr <- PCA(t(datat_clr[,cond]), graph = F)
pca_clr_g <- PCA(t(datag_clr[,cond]), graph=F)
pdf(paste('PCA_clr_meta', ty,'_' ,tax_id,'_',frac,'.pdf', sep=''), width=14)
par(mfrow=c(1,2))
plot.PCA(pca_clr, choix='var', col.var = cols$col[match(rownames(pca_clr$var$coord), cols$taxon)])
plot.PCA(pca_clr, choix='ind')
plot.PCA(pca_clr_g, choix='var', col.var = cols$col[match(rownames(pca_clr_g$var$coord), cols$taxon)])
plot.PCA(pca_clr_g, choix='ind')
dev.off()
test_t <- data.frame(t(apply(datat_clr, 1, wilc_test)), stringsAsFactors = F)
test_t$X3 <- p.adjust(test_t$X1, 'hommel')
test_t$loc <- sapply(1:length(test_t$X3),FUN = function(x){ifelse(test_t$X3[x]>0.05, yes='0', no=test_t$X2[x])})
test_g <- data.frame(t(apply(datag_clr, 1, wilc_test)), stringsAsFactors = F)
test_g$X3 <- p.adjust(test_g$X1, 'hommel')
test_g$loc <- sapply(1:length(test_g$X3),FUN = function(x){ifelse(test_g$X3[x]>0.05, yes='0', no=test_g$X2[x])})

test_t0 <- data.frame(t(apply(log(datat), 1, wilc_test)), stringsAsFactors = F)
test_t0$X3 <- p.adjust(test_t0$X1, 'hommel')
test_t0$loc <- sapply(1:length(test_t0$X3),FUN = function(x){ifelse(test_t0$X3[x]>0.05, yes='0', no=test_t0$X2[x])})
test_g0 <- data.frame(t(apply(log(datag), 1, wilc_test)), stringsAsFactors = F)
test_g0$X3 <- p.adjust(test_g0$X1, 'hommel')
test_g0$loc <- sapply(1:length(test_g0$X3),FUN = function(x){ifelse(test_g0$X3[x]>0.05, yes='0', no=test_g0$X2[x])})

logratio_func <- function(x){
  na <- rownames(x)
  x <- as.numeric(x)
  lr <-median(x[colnames(datat) %in% arctic_stations])-median(x[colnames(datat) %in% atlantic_stations])
  return(c(na, lr))
}
lr_t <- t(apply(log(datat), 1, logratio_func))
lr_t <- rbind(lr_t, t(apply(log(datag), 1, logratio_func)))

rownames(lr_t)<-c('metatransciptomes', 'metagenomes')
pvals <- rbind(test_t0$X3, test_g0$X3)
#set <-magma(50)
#set <- colorRampPalette(colors = c("red", "yellow", 'blue'))(50)
set <- colorRampPalette(colors = c("orange4", 'gray', 'darkorchid'))(50)
pdf(paste('heatmap_metaTG_ratios_', frac,'_' , tax_id,'.pdf',sep=''), height=3)
heatmap.2(lr_t, trace="none",symm=TRUE, Rowv = NA, Colv = NA,cexRow=0.5, cexCol=0.5,
          dendrogram = "none", keysize=4, col = set,margins = c(7,5),
          symkey = F, br=seq(-max(abs(c(lr_t)), na.rm=T), max(abs(c(lr_t)), na.rm=T),length.out =  51),denscol = 'black',
          cellnote = ifelse(pvals>0.05, NA, '*'),notecol='red', lhei=rep(1, nrow(lr_t)), 
          notecex = 1 , sepwidth=c(2,0.5)
)
dev.off()




clr_t<- t(apply(datat_clr, 1, logratio_func))
clr_t <- rbind(clr_t, t(apply(datag_clr, 1, logratio_func)))

rownames(clr_t)<-c('metatransciptomes', 'metagenomes')
pvals <- rbind(test_t$X3, test_g$X3)
#set <-magma(50)
pdf(paste('heatmap_metaTG_clr_ratios_', frac,'_' , tax_id,'.pdf',sep=''), height=3)
heatmap.2(clr_t, trace="none",symm=TRUE, Rowv = NA, Colv = NA,cexRow=0.5, cexCol=0.5,
          dendrogram = "none", keysize=4, col = set,margins = c(7,5),
          symkey = F, br=seq(-max(abs(c(clr_t)), na.rm=T), max(abs(c(clr_t)), na.rm=T),length.out =  51),denscol = 'black',
          cellnote = ifelse(pvals>0.05, NA, '*'),notecol='red', lhei=rep(1, nrow(lr_t)),
          notecex = 1 , sepwidth=c(2,0.5)
)
dev.off()
