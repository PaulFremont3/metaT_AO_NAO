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
bis <- commandArgs(trailingOnly = T)[8]

if (bis=='0'){
 bis=''
}
u=read.table(paste('phate_training_expression_',frac,'_',zeros,'_',norm,'_',sub2,'.txt', sep=''))
v=read.table(paste('phate_fit_expression_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'.txt', sep=''))
unigenes <- readLines(paste('unigenes_phate_',frac,'_',zeros,'_',norm,'_',sub2,'.txt', sep=''))
unigenes <- strsplit(unigenes, ' ')[[1]]


u1=readRDS(paste('sub_phate_fit_expression_py_subtax_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.rds', sep=''))
tax_pfam_file <- readRDS(paste('sub_phate_pfam_taxo_subtax_', knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,bis,'.rds', sep=''))
uni_sel= unigenes %in% tax_pfam_file$uid
ub=u[uni_sel,]
vb=v[uni_sel,]

dat_t=list(u,u, ub, ub )
dat_f=list(v, v, vb, vb)

funco=function(x){
  sum(x!=0)
}

func_perc=function(x){
  round(100*sum(x!=0)/length(x)+1, 0)
}

data=NULL
to_plot=list(NULL, NULL,NULL, NULL, NULL, NULL, NULL)
to_plot1=list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
for (i in 1:4){
  ui=dat_t[[i]]
  vi=dat_f[[i]]
  if (i %in% c(2, 4)){
    sel1=1:8
    ui1=ui[,sel1]
    sel2=9:20
    ui2=ui[,sel2]
    means1 = apply(ui1, 1, mean)
    means2 = apply(ui2, 1, mean)
    means=abs(means1-means2)
    means3=means1-means2
 
    means4 = apply(ui1, 1, function(x){mean(x[x!=0])})
    means5 = apply(ui2, 1, function(x){mean(x[x!=0])})
    means_ter=means4-means5   

    
    perc4=apply(ui1, 1, func_perc)
    perc5=apply(ui2, 1, func_perc)
    percs= log10(perc4/perc5)
    #means1=abs(means1-0)
    #means2=abs(means2-0)
  } else {
    means = apply(ui, 1, mean)
    count_stats=apply(ui, 1, funco)
    #means3= apply(ui, 1, mean)
  }
  vars=apply(ui, 1, function(x){sd(x[x!=0])})
  vars1=apply(ui, 1,sd)
  #medians = apply(ui, 1, median)
  means_bis = apply(ui, 1, function(x){mean(x[x!=0])})
  #medians_bis = apply(ui, 1, function(x){median(x[x!=0])})
  if (i ==1){
    to_plot[[1]] = means
    to_plot[[2]] = vi[,1]
    to_plot[[3]] = vi[,2]
    to_plot[[4]] = vars
    to_plot[[5]] = count_stats
    to_plot[[6]] = means_bis
    to_plot[[7]] = vars1
  }
  if (i ==2){
    to_plot1[[1]] = means
    to_plot1[[2]] = vi[,1]
    to_plot1[[3]] = vi[,2]
    to_plot1[[4]] = means3
    to_plot1[[5]] = percs
    to_plot1[[6]] = means_bis
    to_plot1[[7]] = means_ter
  }  
  c1=cor.test(means, vi[,1])
  #c2=cor.test(medians, vi[,1])
  c2=cor.test(means, vi[,2])
  c3=cor.test(means, abs(vi[,1]))
  c4=cor.test(means, abs(vi[,2]))
  c5=cor.test(vars, abs(vi[,1]))
  c6=cor.test(vars, abs(vi[,2]))
  #c4=cor.test(medians, vi[,2])
  print(c(c1$estimate, c1$p.value, c2$estimate, c2$p.value, c3$estimate, c3$p.value, c4$estimate, c4$p.value, c5$estimate, c5$p.value, c6$estimate, c6$p.value))  

  #c5=cor.test(means_bis, vi[,1])
  #c6=cor.test(medians_bis, vi[,2])
  #c7=cor.test(medians_bis, vi[,1])
  #c8=cor.test(medians_bis, vi[,2])
   
  data=c(data, c(c1$estimate, c1$p.value, c2$estimate, c2$p.value, c3$estimate, c3$p.value, c4$estimate, c4$p.value, c5$estimate, c5$p.value, c6$estimate, c6$p.value))
}



data=matrix(data, ncol=2, byrow=T)
axis=rep(c(1,2), 6*2)
type=c(rep('all', 12), rep('sub', 12))
type2=rep( c( rep('mean', 4),rep('var', 2) ,rep('delta_mean', 4), rep('var', 2) ), 2)
axis_val=rep(c('signed', 'signed', 'abs', 'abs', 'abs', 'abs'), 4)
#zeros=rep(c('yes','yes','yes','yes', 'no', 'no','no','no'), 4)
#domain=rep(c(rep('all', 4), rep('atlantic', 4), rep('arctic', 4)) , 2)
data=as.data.frame(data)
data$axis=axis
data$axis_val=axis_val
data$type=type
data$type2=type2
#data$domain=domain
write.table(data, paste('phate_axis_correlations_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo, '.txt', sep=''), row.names=F)

#pdf('phate_plots_and_correlations.pdf')
#plot(to_plot[[2]], to_plot[[3]], xlab='PHATE 1', ylab='PHATE 2', cex=0.5, col=scales::alpha('black', 0.3))
#plot(to_plot[[2]], to_plot[[1]], xlab='PHATE 1', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
#plot(to_plot[[3]], to_plot[[1]], xlab='PHATE 2', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
#plot(abs(to_plot1[[2]]), to_plot1[[1]], xlab='|PHATE 1|', ylab='|delta| of mean expression (AO/NAO)', cex=0.5, col=scales::alpha('black', 0.3))
#plot(abs(to_plot1[[3]]), to_plot1[[1]], xlab='|PHATE 2|', ylab='|delta| of mean expression (AO/NAO)', cex=0.5, col=scales::alpha('black', 0.3))
#dev.off()
saveRDS(to_plot, 'to_plot_phate.rds')
saveRDS(to_plot1, 'to_plot1_phate.rds')

col_set=colorRampPalette(c("darkblue","blue", "gray", "red", "darkred"))(99)

#col_set1=magma(99)
col_set1=colorRampPalette(c("gray", 'lightblue', 'blue', 'darkblue'))(99)
col_set0=colorRampPalette(c("darkturquoise","turquoise", "gray", "orange", "darkorange"))(99)
col_set2=colorRampPalette(c( 'deepskyblue', 'blue', 'orange', 'red'))(15)
col_set3=colorRampPalette(c( 'cyan', 'deepskyblue','gray' ,'orange', 'firebrick'))(99)

plot_f=function(data, cols, a,b,c,d){
  plot(0, 0,xlim=c(a,b),ylim=c(c,d), xlab='PHATE 1',
     ylab='PHATE 2', col=scales::alpha('white', 0), cex=0.5,
     pch=19, cex.axis=1.3, cex.lab=1.3)
  points(data$x, data$y, pch=20, cex=0.2,col=scales::alpha(cols[data$co]) )
}

pdf('phate_plots_and_correlations_density_bis.pdf')
a=min(c(to_plot[[2]]))
b=max(c(to_plot[[2]]))
c=min(c(to_plot[[3]]))
d=max(c(to_plot[[3]]))
data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot[[1]])
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
mino= -maxo
quart1=(maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
data$co = round( 99*(data$z-mino)/(maxo- mino)+1 )
ordi=order(abs(data$z), decreasing=F)
data=data[ordi,]
plot_f(data, col_set,a,b,c,d)

a=min(c(to_plot[[2]]))
b=max(c(to_plot[[2]]))
c=min(c(to_plot[[3]]))
d=max(c(to_plot[[3]]))
data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot1[[6]])
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
mino= -maxo
quart1=(maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
data$co = round( 99*(data$z-mino)/(maxo- mino)+1 )
ordi=order(abs(data$z), decreasing=F)
data=data[ordi,]
plot_f(data, col_set,a,b,c,d)

data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot1[[4]])
#maxo= max(data$z, na.rm=T)
#mino= min(data$z, na.rm=T)
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
data$co = round( 99*(data$z-mino)/(maxo- mino)+1 )
ordi=order(abs(data$z), decreasing=F)
data=data[ordi,]
plot_f(data, col_set0,a,b,c,d)

data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot1[[7]])
#maxo= max(data$z, na.rm=T)
#mino= min(data$z, na.rm=T)
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
data$co = round( 99*(data$z-mino)/(maxo- mino)+1 )
ordi=order(abs(data$z), decreasing=F)
data=data[ordi,]
plot_f(data, col_set0,a,b,c,d)



#data=data.frame(x=to_plot[[2]], y=to_plot[[3]], co=to_plot[[5]]-4 )
#print(max(to_plot[[5]]-4))
#print(min(to_plot[[5]]-4))
#ordi=order(data$co, decreasing=F)
#data=data[ordi,]
#plot_f(data,col_set2,a,b,c,d)


zo=to_plot[[7]]
zo[zo>5.2]=5.2
data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=zo)
maxo= max(data$z, na.rm=T)
mino= 0
#maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
#mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
data$co = round( 99*(data$z-mino)/(max(data$z)- mino)+1 )
ordi=order(abs(data$z), decreasing=F)
data=data[ordi,]
plot_f(data, col_set1,a,b,c,d)

#data=data.frame(x=to_plot[[2]], y=to_plot[[3]], co=to_plot[[5]]-4 )
#print(max(to_plot[[5]]-4))
#print(min(to_plot[[5]]-4))
#ordi=order(data$co, decreasing=F)
#data=data[ordi,]
#plot_f(data,col_set2,a,b,c,d)

zo=to_plot[[4]]
zo[zo>5.2]=5.2
data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=zo)
maxo= max(data$z, na.rm=T)
mino= 0
#maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
#mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
data$co = round( 99*(data$z-mino)/(max(data$z)- mino)+1 )
ordi=order(abs(data$z), decreasing=F)
data=data[ordi,]
plot_f(data, col_set1,a,b,c,d)

data=data.frame(x=to_plot[[2]], y=to_plot[[3]], co=to_plot[[5]]-4 )
#print(max(to_plot[[5]]-4))
#print(min(to_plot[[5]]-4))
ordi=order(data$co, decreasing=F)
data=data[ordi,]
plot_f(data,col_set2,a,b,c,d) 

data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot1[[5]])
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))) )
mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
data$co = round( 99*(data$z-mino)/(maxo- mino)+1 )
ordi=order(abs(data$z), decreasing=F)
data=data[ordi,]
plot_f(data, col_set3,a,b,c,d)
dev.off()


legend_function = function(bks,colos, mino, maxo ){
  ft=data.frame(x=0, y=0)
  ggobj=ggplot(data = ft, mapping = aes(x = x, y = y,color=y)) +
  geom_point(size=0.3, alpha=1)+
  scale_color_gradientn(colours=colos, na.value = "transparent",
                       breaks=bks,labels=bks,
                       limits=c(mino,maxo)) + theme_bw(base_size = 23) +theme(panel.grid.minor = element_blank())
  print(ggobj)
}

legs=rep(list(NULL), 8)
pdf('phate_plots_and_expression_values.pdf')
data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot[[1]])
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
mino= -maxo
quart1=(maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino

bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))
a=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
  geom_point(size=0.3, alpha=0.3) +
  scale_color_gradientn(colours=col_set, na.value = "transparent",
                       breaks=bks,labels=bks,
                       limits=c(mino,maxo)) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23) +theme(panel.grid.minor = element_blank(), legend.position="none")
print(a)

legs[[1]]=list(bks, col_set, mino, maxo)

data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot1[[6]])
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
mino= -maxo
quart1=(maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino

bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))
a=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
  geom_point(size=0.3, alpha=0.3) +
  scale_color_gradientn(colours=col_set, na.value = "transparent",
                       breaks=bks,labels=bks,
                       limits=c(mino,maxo)) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23) +theme(panel.grid.minor = element_blank(), legend.position="none")
print(a)

legs[[2]]=list(bks, col_set, mino, maxo)

print('cors with zeroes:')
print('delta ao-nao versus sd')
print(cor.test( abs(to_plot1[[4]]), to_plot[[7]] ))
print('delta ao-nao versus mean expr')
print(cor.test( abs(to_plot1[[4]]), abs(to_plot[[1]]) ))
print('sd versus mean expr')
print(cor.test( abs(to_plot[[1]]) , to_plot[[7]] ))
print('mean expr vs phate 1')
print(cor.test(to_plot[[2]],to_plot[[1]]))
print('cors with no zeroes:')
print('delta ao-nao versus sd')
print(cor.test( abs(to_plot1[[7]]), to_plot[[4]] ))
print('delta ao-nao versus mean expr')
print(cor.test( abs(to_plot1[[7]]), abs(to_plot1[[6]]) ))
print('sd versus mean expr')
print(cor.test( abs(to_plot1[[6]]) , to_plot[[4]] ))
print('mean expr vs phate 1')
print(cor.test(to_plot[[2]],to_plot1[[6]]))

print('mean expr vs phate 1 (zeros):')
print(cor.test( to_plot[[2]] ,to_plot[[1]] ))
print('mean expr vs phate 1 (no zeros):')
print(cor.test( to_plot[[2]] ,to_plot1[[6]] ))
print('abs delta mean expr vs phate 1 (zeros):')
print(cor.test( to_plot[[2]] ,abs(to_plot1[[4]]) ))
print('abs delta mean expr vs phate 1 (no zeros):')
print(cor.test( to_plot[[2]] ,abs(to_plot1[[7]]) ))
print('biogeography vs phate 2:')
print(cor.test(to_plot[[3]],to_plot1[[5]]))

data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot1[[4]])
#maxo= max(data$z, na.rm=T)
#mino= min(data$z, na.rm=T)
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))
b=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
  geom_point(size=0.3, alpha=0.3) +
  scale_color_gradientn(colours=col_set0,na.value = "transparent",
                       breaks=bks,labels=bks,
                       limits=c(mino,maxo) ) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23) +theme(panel.grid.minor = element_blank(), legend.position="none")
print(b)
legs[[3]]=list(bks, col_set0, mino, maxo)

data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot1[[7]])
#maxo= max(data$z, na.rm=T)
#mino= min(data$z, na.rm=T)
maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))
b=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
  geom_point(size=0.3, alpha=0.3) +
  scale_color_gradientn(colours=col_set0,na.value = "transparent",
                       breaks=bks,labels=bks,
                       limits=c(mino,maxo) ) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23) +theme(panel.grid.minor = element_blank(), legend.position="none")
print(b)
legs[[4]]=list(bks, col_set0, mino, maxo)

zo=to_plot[[7]]
zo[zo>5.2]=5.2
data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=zo)
maxo= max(data$z, na.rm=T)
mino= min(data$z, na.rm=T)
#maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
#mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))
b1=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
  geom_point(size=0.3, alpha=0.3) +
  scale_color_gradientn(colours=col_set1,na.value = "transparent",
                       breaks=bks,labels=bks,
                       limits=c(mino,maxo) ) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23 ) +theme(panel.grid.minor = element_blank(), legend.position="none")
print(b1)
legs[[5]]=list(bks, col_set1, mino, maxo)

zo=to_plot[[4]]
zo[zo>5.2]=5.2
data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=zo)
maxo= max(data$z, na.rm=T)
mino= min(data$z, na.rm=T)
#maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
#mino= -maxo
quart1= (maxo-mino)/4 +mino
quart= abs(3*(maxo-mino)/4 +mino)
mid= (maxo-mino)/2 +mino
bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))
b1=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
  geom_point(size=0.3, alpha=0.3) +
  scale_color_gradientn(colours=col_set1,na.value = "transparent",
                       breaks=bks,labels=bks,
                       limits=c(mino,maxo) ) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23 ) +theme(panel.grid.minor = element_blank(), legend.position="none")
print(b1)
legs[[6]]=list(bks, col_set1, mino, maxo)
dev.off()

legs[[7]] =list(c(5:20), col_set2, 5, 20)
legs[[8]] =list(c(-2, -1, 0,1,2), col_set3, -2, 2)
pdf('legend_phate.pdf')
for (k in 1:8){
  legend_function(legs[[k]][[1]], legs[[k]][[2]],legs[[k]][[3]],legs[[k]][[4]])
}
dev.off()


My_Theme = theme(
  axis.title.x = element_text(size = 23),
  axis.text.x = element_text(size = 23),
  axis.title.y = element_text(size = 23), 
  axis.text.y = element_text(size = 23),
  legend.text = element_text(size = 23),
  legend.title =element_text(size = 23),
  panel.background = element_rect(fill = "white", colour = "grey50"))

pdf('phate_plots_and_correlations_density.pdf', width=10)
data=data.frame(x=to_plot[[2]], y=to_plot[[3]])# z=to_plot[[1]], w=to_plot1[[1]])
#plot(to_plot[[2]], to_plot[[3]], xlab='PHATE 1', ylab='PHATE 2', cex=0.5, col=scales::alpha('black', 0.3))
ab=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3)  +
  scale_color_viridis_c(option='magma') +labs(x ='PHATE 1', y = 'PHATE 2') + My_Theme
print(ab)



#data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot[[1]])
#maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
#mino= -maxo
#quart1=(maxo-mino)/4 +mino
#quart= abs(3*(maxo-mino)/4 +mino)
#mid= (maxo-mino)/2 +mino

#bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))

#a=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
#  geom_point(size=0.3, alpha=0.3) +
#  scale_color_gradientn(colours=col_set, na.value = "transparent",
#                       breaks=bks,labels=bks,
#                       limits=c(mino,maxo)) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23) +theme(panel.grid.minor = element_blank(), legend.position="none")
#print(a)

#print(cor.test( abs(to_plot1[[4]]), to_plot[[4]] ))
#print(cor.test( abs(to_plot1[[4]]), abs(to_plot[[1]]) ))
#print(cor.test( abs(to_plot[[1]]) , to_plot[[4]] ))
#data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=to_plot1[[4]])
#maxo= max(data$z, na.rm=T)
#mino= min(data$z, na.rm=T)
#maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
#mino= -maxo
#quart1= (maxo-mino)/4 +mino
#quart= abs(3*(maxo-mino)/4 +mino)
#mid= (maxo-mino)/2 +mino
#bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))
#b=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
#  geom_point(size=0.3, alpha=0.3) +
#  scale_color_gradientn(colours=col_set0,na.value = "transparent",
#                       breaks=bks,labels=bks,
#                       limits=c(mino,maxo) ) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23) +theme(panel.grid.minor = element_blank(), legend.position="none")
#print(b)

#zo=to_plot[[4]]
#zo[zo>5.2]=5.2
#data=data.frame(x=to_plot[[2]], y=to_plot[[3]], z=zo)
#maxo= max(data$z, na.rm=T)
#mino= min(data$z, na.rm=T)
#maxo= max( c( abs(min(data$z, na.rm=T)) , abs(max(data$z, na.rm=T))))
#mino= -maxo
#quart1= (maxo-mino)/4 +mino
#quart= abs(3*(maxo-mino)/4 +mino)
#mid= (maxo-mino)/2 +mino
#bks=c(round(mino, 1), round(quart1, 1),round(mid, 1),round(quart, 1),round(maxo, 1))
#b1=ggplot(data = data, mapping = aes(x = x, y = y, color=z)) +
#  geom_point(size=0.3, alpha=0.3) +
#  scale_color_gradientn(colours=col_set1,na.value = "transparent",
#                       breaks=bks,labels=bks,
#                       limits=c(mino,maxo) ) +labs(x ='PHATE 1', y = 'PHATE 2') + theme_bw(base_size = 23 ) +theme(panel.grid.minor = element_blank(), legend.position="none")
#print(b1)

data=data.frame(x=to_plot[[7]], y=abs(to_plot[[1]] ) )
#plot(to_plot[[2]], to_plot[[1]], xlab='PHATE 1', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
b2=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='plasma') +labs(x ='sd', y = '|mean expression|') + My_Theme
print(b2)

data=data.frame(x=to_plot[[7]], y=to_plot[[1]]  )
#plot(to_plot[[2]], to_plot[[1]], xlab='PHATE 1', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
b3=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='plasma') +labs(x ='sd', y = 'mean expression') + My_Theme
print(b3)

data=data.frame(x=to_plot[[4]], y=abs( to_plot1[[6]] ) )
#plot(to_plot[[2]], to_plot[[1]], xlab='PHATE 1', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
b2=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='plasma') +labs(x ='sd', y = '|mean expression|') + My_Theme
print(b2)

data=data.frame(x=to_plot[[4]], y=to_plot1[[6]]  )
#plot(to_plot[[2]], to_plot[[1]], xlab='PHATE 1', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
b3=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='plasma') +labs(x ='sd', y = 'mean expression') + My_Theme
print(b3)


data=data.frame(x=to_plot[[2]], y=to_plot[[1]])
#plot(to_plot[[2]], to_plot[[1]], xlab='PHATE 1', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
c=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='inferno') +labs(x ='PHATE 1', y = 'mean expression') + My_Theme
print(c)

data=data.frame(x=to_plot[[2]], y=to_plot1[[6]])
#plot(to_plot[[2]], to_plot[[1]], xlab='PHATE 1', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
c=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='inferno') +labs(x ='PHATE 1', y = 'mean expression') + My_Theme
print(c)

data=data.frame(x=to_plot[[3]], y=to_plot1[[5]])
#plot(to_plot[[2]], to_plot[[1]], xlab='PHATE 1', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
c=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='viridis') +labs(x ='PHATE 2', y = '') + My_Theme
print(c)


data=data.frame(x=to_plot[[3]], y=to_plot[[1]])
#plot(to_plot[[3]], to_plot[[1]], xlab='PHATE 2', ylab='mean expression', cex=0.5, col=scales::alpha('black', 0.3))
d=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='plasma') +labs(x ='PHATE 2', y = 'mean expression') + My_Theme
print(d)

data=data.frame(x=abs(to_plot1[[2]]), y=to_plot1[[1]])
#plot(abs(to_plot1[[2]]), to_plot1[[1]], xlab='|PHATE 1|', ylab='|delta| of mean expression (AO/NAO)', cex=0.5, col=scales::alpha('black', 0.3))
e=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='viridis') +labs(x ='|PHATE 1|', y = '|delta| of mean expression (AO/NAO)') + My_Theme
print(e)

data=data.frame(x=abs(to_plot1[[3]]), y=to_plot1[[1]])
#plot(abs(to_plot1[[3]]), to_plot1[[1]], xlab='|PHATE 2|', ylab='|delta| of mean expression (AO/NAO)', cex=0.5, col=scales::alpha('black', 0.3))
f=ggplot(data = data, mapping = aes(x = x, y = y)) +
  geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option='cividis') +labs(x ='|PHATE 2|', y = '|delta| of mean expression (AO/NAO)') + My_Theme
print(f)
dev.off()

