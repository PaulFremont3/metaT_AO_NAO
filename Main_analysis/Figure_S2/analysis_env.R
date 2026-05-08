#!/bin/env/usr/env Rscript
library('treemap')
#library("readxl")
library('gplots')
library('FactoMineR')
library('stringr')
library('RColorBrewer')
# setwd("~/Groups_metaT")
#setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4/Groups_metaT/")

n_var <- 9
env_data <- read.table('env_arctic_3.txt', header = T)
colnames(env_data)[11:13] <- c('Si*', 'N*', 'Iron_PISCESv2')
env_data <- env_data[,c(-6,-11, -12,-14,-15 )]

env_data <- env_data[!(env_data$Station %in% c(142,143,146,149,191,201,205, 206, 208, 209, 210)),]
env_data$season=c(rep(1, 7), rep(2,3), rep(3, 10))
cor_map <- cor(env_data[,4:(4+n_var-1)])
pdf('variables_correlation_map.pdf')
heatmap.2(cor_map, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1,margins=c(10,10),
          labCol=colnames(env_data)[4:(4+n_var-1)],
          labRow=colnames(env_data)[4:(4+n_var-1)], col= redgreen(200), symkey = F, density.info = 'none')
v <- env_data[,c(4:(4+n_var-1))]
rownames(v)<-env_data$Station
set.seed(1)
u <- PCA(v, graph = F)
res.km <- kmeans(u$ind$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
u$ind$coord[,1] <- (u$ind$coord[,1]-min(u$ind$coord[,1]))*2/(max(u$ind$coord[,1])-min(u$ind$coord[,1]))-1
u$ind$coord[,2] <- (u$ind$coord[,2]-min(u$ind$coord[,2]))*2/(max(u$ind$coord[,2])-min(u$ind$coord[,2]))-1
cols <- c('red', 'violet', 'lightblue', 'darkblue')
par(mar=c(5.1, 6.1, 4.1, 2.1))
#plot.PCA(u, axes=c(1, 2), choix="ind", habillage="none", col.ind = cols[grp], cex.lab=2, cex.axis=2, cex=1.5)
# dat <- NULL
# for (i in 1:2){
#   dat <- cbind(dat, spread.labs(u$ind$coord[,i], mindiff = 0.08))
# }
# text(dat[,1], dat[,2],labels=rownames(u$ind$coord), cex=1.5)
#arrows(x0=rep(0, 10),y0= rep(0, 10),x1=u$var$coord[,1], y1=u$var$coord[,2])
#text(u$var$coord[,1]+0.1*sign(u$var$coord[,1]) , u$var$coord[,2]+0.1*sign(u$var$coord[,2]),
#     labels=colnames(env_data)[4:(4+n_var-1)], cex=1.5)

grp1 <- rownames(u$ind$coord)
grp0 <- rownames(u$ind$coord)
grp0[grp1>'155'] <- 3
grp0[grp1 <='155'] <- 1
#grp0[grp1>='155' & grp1<='163']<- 2
cols0 <- c('darkorange', 'red', 'darkblue')
e1=round(100*u$eig[1,1]/sum(u$eig[,1]), 1)
e2=round(100*u$eig[2,1]/sum(u$eig[,1]), 1)
par(mar=c(5.1, 7.1, 4.1, 2.1))
plot(u$ind$coord[,1], u$ind$coord[,2], col=cols0[as.numeric(grp0)], cex.lab=2, cex.axis=2, cex=1.5, pch=16,
     xlab = paste('PCA 1 (',e1,'%)', sep=''), ylab=paste('PCA 2 (',e2,'%)', sep=''))
#plot.PCA(u, axes=c(1, 2), choix="ind", habillage="none", 
#         col.ind = cols0[as.numeric(grp0)], cex.lab=2, cex.axis=2, cex=1.5)
# dat <- NULL
# for (i in 1:2){
#   dat <- cbind(dat, spread.labs(u$ind$coord[,i], mindiff = 0.08))
# }
arrows(x0=rep(0, 10),y0= rep(0, 10),x1=u$var$coord[,1], y1=u$var$coord[,2])
text(u$var$coord[,1]+0.1*sign(u$var$coord[,1]) , u$var$coord[,2]+0.1*sign(u$var$coord[,2]),
     labels=colnames(env_data)[4:(4+n_var-1)], cex=1.5)
text(u$ind$coord[,1]+0.05*sign(u$ind$coord[,1]), u$ind$coord[,2]+0.05*sign(u$ind$coord[,2]),labels=env_data$Station,
     col=cols0[as.numeric(grp0)])

cols1=c('blue', 'green','red')
plot(u$ind$coord[,1], u$ind$coord[,2], col=cols1[env_data$season], cex.lab=2, cex.axis=2, cex=1.5, pch=16,
     xlab = paste('PCA 1 (',e1,'%)', sep=''), ylab=paste('PCA 2 (',e2,'%)', sep=''))
arrows(x0=rep(0, 10),y0= rep(0, 10),x1=u$var$coord[,1], y1=u$var$coord[,2])
text(u$var$coord[,1]+0.1*sign(u$var$coord[,1]) , u$var$coord[,2]+0.1*sign(u$var$coord[,2]),
     labels=colnames(env_data)[4:(4+n_var-1)], cex=1.5)
text(u$ind$coord[,1]+0.05*sign(u$ind$coord[,1]), u$ind$coord[,2]+0.05*sign(u$ind$coord[,2]),labels=env_data$Station,
     col=cols0[as.numeric(grp0)])

cols2=hcl.colors(12, "Set3")
plot(u$ind$coord[,1], u$ind$coord[,2], col=cols2[env_data$Month], cex.lab=2, cex.axis=2, cex=1.5, pch=16,
     xlab = paste('PCA 1 (',e1,'%)', sep=''), ylab=paste('PCA 2 (',e2,'%)', sep=''))
arrows(x0=rep(0, 10),y0= rep(0, 10),x1=u$var$coord[,1], y1=u$var$coord[,2])
text(u$var$coord[,1]+0.1*sign(u$var$coord[,1]) , u$var$coord[,2]+0.1*sign(u$var$coord[,2]),
     labels=colnames(env_data)[4:(4+n_var-1)], cex=1.5)
text(u$ind$coord[,1]+0.05*sign(u$ind$coord[,1]), u$ind$coord[,2]+0.05*sign(u$ind$coord[,2]),labels=env_data$Station,
     col=cols0[as.numeric(grp0)])
dev.off()


h <- v
sel_var=c(1:6,9:11)
hi=h[,sel_var]
c=1
for (i in sel_var){
  hi[,c] <- (h[,i]-mean(h[,i]))/sd(h[,i])
  c=c+1
}
env_dist <- as.matrix(dist(hi))
saveRDS(env_dist, 'environmental_distances_arctic_3.rds')

list_vars <- list(c(1,2,3), c(4:8), c(9), c(10:11))
set.seed(3)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- rep(as.vector(col_vector), 2)
col_vector <- sample(col_vector)
col_vector <- c('red', 'lightblue', 'darkgreen', 'darkorange', 'cyan', 'magenta4',
                'saddlebrown', 'peru', 'royalblue4', 'orangered', 'mediumorchid1',
                'forestgreen', 'firebrick4', 'darkturquoise')
#col_vector[1] <- 'black'
list_v <- c('Physical', 'Nutrients', 'Iron',  'Seasonality indexes')
pdf('z-scores_variables.pdf', width = 12)
cu=1
for (u in list_vars){
  c=1
  plot(as.factor(env_data$Station),rep(10,length(env_data$Station)), main=list_v[[cu]], 
       las=2, col=scales::alpha('white', 0), ylim = c(-max(abs(h)),max(abs(h))),
       xlim=c(1,length(env_data$Station)))
  for (t in u){
    points(as.factor(env_data$Station), h[,t], pch=t, col=col_vector[t], cex=2)
    
    c=c+1
  }
  cu=cu+1
  legend('topleft', colnames(h)[u], pch = u, col = col_vector[u], ncol=2)
}
dev.off()

