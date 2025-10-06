library('gplots')
library('matlab')

u=read.table('minAijji_tarrive_min_surface_1000.csv')
v=read.table('Lagr_distances.txt', sep='\t', header = T, row.names = 1)

stats=c('144', '145','147', '148', '150', '151', '152', '155', '158', '163', '168', '173', '175', '178', '180', '188', '189', '193', '194', '196' )
stats_n=paste(stats, '_SUR', sep='')
sel=rownames(v) %in% stats

v=v[sel,sel]
rownames(v)=stats
colnames(v)=stats
col=colorRampPalette(c("red", "orange", "yellow", "green", "blue"))(100)
pdf('lagr_time_matrix.pdf')
heatmap.2(x = as.matrix(v),trace="none",symm=F, Rowv = NA, Colv = NA,cexRow=1, cexCol=1,
          dendrogram = "none", keysize=1,margins = c(7,5),col=col,
          symkey = F, br=seq(0, max(u, na.rm=T),length.out =  101))

dev.off()
