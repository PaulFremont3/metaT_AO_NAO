setwd("~/phate_contour_plot")

dat=readRDS('data-phate-contour_plot.rds')
a=min(dat[[1]][,1])
b=max(dat[[1]][,1])
c=min(dat[[1]][,2])
d=max(dat[[1]][,2])

u=readRDS('color_table_groups3.rds')

taxos=c('unknown', "Eukaryota (unclassified)", 
        "other Opisthokonta","Hexanauplia" , "Ciliophora", "other Haptophyta",
        "Phaeocystales", "Pelagophyceae", "Mamiellales", "other Stramenopiles", 
        "Tunicata","Bacillariophyta" , 'Dinophyceae')

no_co=scales::alpha('white', 0)
pdf('contour_phate.pdf')
plot(0, 0,xlim=c(a,b),ylim=c(c,d), xlab='PHATE 1',
     ylab='PHATE 2', col=scales::alpha('white', 0), cex=0.5,
     pch=19, cex.axis=1.3, cex.lab=1.3)
k=1
for (t in taxos){
  #x=data[colos==co,1]
  #y=data[colos==co,2]
  #print(length(x))
  co=u$col[u$taxon==t]
  z <- dat[[2]][[co]]
  #data_cont[[co]] = z
  contour(z, drawlabels=FALSE, nlevels=9,lwd=3, lty = 2,
          col= scales::alpha(rep(co, 8), 0.6), add=TRUE)
  
  k=k+1
}
dev.off()

#scales::alpha(rep(co, 8), 0.5)
pdf('contour_phate_bis.pdf')

k=1
for (t in taxos){
  plot(0, 0,xlim=c(a,b),ylim=c(c,d), xlab='PHATE 1', ylab='PHATE 2', 
       col=scales::alpha('white', 0), cex=0.5, pch=19, cex.axis=1.3, cex.lab=1.3)
  #x=data[colos==co,1]
  #y=data[colos==co,2]
  #print(length(x))
  co=u$col[u$taxon==t]
  z <- dat[[2]][[co]]
  #data_cont[[co]] = z
  contour(z, drawlabels=FALSE, nlevels=9,lwd=3, col= rep(co, 8), add=TRUE)
  
  k=k+1
}
dev.off()


pdf('contour_phate_ter.pdf', width=5, height=5)

k=1
for (t in taxos){
  plot(0, 0,xlim=c(a,b),ylim=c(c,d), xlab='PHATE 1', ylab='PHATE 2', 
       col=scales::alpha('white', 0), cex=0.5, pch=19, cex.axis=1.3, cex.lab=1.3, )
  #x=data[colos==co,1]
  #y=data[colos==co,2]
  #print(length(x))
  co=u$col[u$taxon==t]
  z <- dat[[2]][[co]]
  #data_cont[[co]] = z
  contour(z, drawlabels=FALSE, nlevels=9,lwd=3, col= rep(co, 8), add=TRUE)
  
  k=k+1
}
dev.off()

