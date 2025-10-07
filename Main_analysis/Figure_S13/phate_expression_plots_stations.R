library('stringr')
library('ggplot2')
library('tidyverse')
library('ggpointdensity')
library('paletteer')
library('pals')
library('RColorBrewer')
frac <- commandArgs(trailingOnly = T)[1]
zeros <- commandArgs(trailingOnly = T)[2]
subs <- commandArgs(trailingOnly = T)[3]
norm <- commandArgs(trailingOnly = T)[4]
sub2 <- commandArgs(trailingOnly = T)[5]
knn <- commandArgs(trailingOnly = T)[6]
taxo <- commandArgs(trailingOnly = T)[7]


u=read.table(paste('../data/phate_training_expression_',frac,'_',zeros,'_',norm,'_',sub2,'.txt', sep=''))
v=read.table(paste('../data/phate_fit_expression_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'.txt', sep=''))
unigenes <- readLines(paste('../data/unigenes_phate_',frac,'_',zeros,'_',norm,'_',sub2,'.txt', sep=''))
unigenes <- strsplit(unigenes, ' ')[[1]]


stats <- c(144, 145,  147, 148 ,
             151, 152, 155,158,163, 150,
             168 ,173, 175, 178,180, 188, 189,
             193, 194, 196)

dat <- readRDS(paste('../data/metat_analysis_',frac,'.rds',sep= ''))
dat <- dat[dat[,111]>4,]

dat <- dat[,c(1,49, 110, 112)]
new_data <- matrix(NA, ncol=length(stats), nrow=dim(dat)[1])
unigenes_ids <- dat$V1
for (j in 1:dim(dat)[1]){
  stts <- as.numeric(strsplit(as.character(dat$V110[j]), '_')[[1]])
  exprs <- as.numeric(strsplit(as.character(dat$V112[j]), '_')[[1]])
  new_data[j, match(stts, stats)] <- exprs
}


mx=max(u, na.rm=T)
mi=min(u, na.rm=T)
mxo=max(abs(c(mx, mi)))
mio=-mxo


plot_f=function(data, cols, a,b,c,d, sta){
  plot(0, 0,xlim=c(a,b),ylim=c(c,d), xlab='PHATE 1',
     ylab='PHATE 2', col=scales::alpha('white', 0), cex=0.5,
     pch=19, cex.axis=1.3, cex.lab=1.3, main=sta)
  points(data$x, data$y, pch=20, cex=0.2,col=cols[data$co] )
}

legend_function = function(bks,colos, mino, maxo ){
  ft=data.frame(x=0, y=0)
  ggobj=ggplot(data = ft, mapping = aes(x = x, y = y,color=y)) +
  geom_point(size=0.3, alpha=1)+
  scale_color_gradientn(colours=colos, na.value = "transparent",
                       breaks=bks,labels=bks,
                       limits=c(mino,maxo)) + theme_bw(base_size = 23) +theme(panel.grid.minor = element_blank())
  print(ggobj)
}


quart1=(mxo-mio)/4 +mio
quart= abs(3*(mxo-mio)/4 +mio)
mid= (mxo-mio)/2 +mio

bks=c(round(mio, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(mxo, 1))

col_set3=colorRampPalette(c( "darkblue","blue", "gray", "red", "darkred"))(100)
col_set4=c(col_set3, 'black')

pdf('legend_phate_expression_stations.pdf')
legend_function(bks,col_set3,mio,mxo)
dev.off()


a=min(v[,1])
b=max(v[,1])
c=min(v[,2])
d=max(v[,2])



for (i in 1:dim(u)[2]){
  st=stats[i]
  vals=u[,i]
  ord=order(abs(vals) , decreasing=F)
  vals=vals[ord]
  co=round( (vals-mio)*99/(mxo-mio)+1 )
  
  v2=new_data[,i]
  ord2=order(abs(v2) , decreasing=F)
  
  vals2=v2[ord2]
  co2=round( (vals2-mio)*99/(mxo-mio)+1 )  
  co2[is.na(co2)] = 101
  id=which(is.na(vals2))
  reord2= c(id:length(vals2) , 1:(id-1))
  co2=co2[reord2]

  dataf=data.frame(x=v[ord,1], y=v[ord,2], co=co)
  pdf(paste('expression_phate_station_', st, ".pdf", sep=''))
  plot_f(dataf, col_set3,a,b,c,d, st)
  dev.off()

  lc=length(co2)
  print(co2[(lc-100):lc])
  phate_coords=v
  phate_coords=phate_coords[ord2,]
  phate_coords=phate_coords[reord2,]
  dataf=data.frame(x=phate_coords[,1], y=phate_coords[,2], co=co2)
  pdf(paste('expression_phate_station_zeros_', st, ".pdf", sep=''))
  plot_f(dataf, col_set4,a,b,c,d, st)
  dev.off()
}


