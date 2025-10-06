library('MASS')
library('RColorBrewer')
library('ggplot2')
library('plyr')
library('shipunov')
frac <- commandArgs(trailingOnly = T)[1]
zeros <- commandArgs(trailingOnly = T)[2]
subs <- commandArgs(trailingOnly = T)[3]
norm <- commandArgs(trailingOnly = T)[4]
sub2 <- commandArgs(trailingOnly = T)[5]
knn <- commandArgs(trailingOnly = T)[6]
taxo <- commandArgs(trailingOnly = T)[7]

photosynthesis_pfams <- c('PF00504' ,'PF00127','PF11623','PF02276' , 'PF01716','PF00556' ,'PF00124', 'PF03967', 'PF00223', 'PF02605', 'PF07465', 'PF05479', 'PF00737', 'PF02532', 'PF02533', 'PF05151'
                          , 'PF02468', 'PF04725', 'PF01405', 'PF06514', 'PF07123', 'PF06596', 'PF06298' , 'PF00421', 'PF05969', 'PF13326', 'PF00796',
                          'PF02427', 'PF02507', 'PF03244', 'PF01701', 'PF01241', 'PF10657', 'PF11947')
#nitrogen_pfams <-c('PF04891', 'PF03206', 'PF04319', 'PF06988', 'PF07732', 'PF13473', 'PF00115',
#                   'PF04879', 'PF02335', 'PF01077','PF03460', 'PF00384',
#                   'PF01568', 'PF04324', 'PF00174', 'PF13435', 'PF07992', 'PF00355', 'PF02665')

nitrogen_pfams <- c('PF02665', 'PF01077','PF03460','PF08376'
                     ,'PF07732', 'PF00384',
                   'PF01568', 'PF04324', 'PF00174', 'PF07992', 'PF00355') 
#carbon_fixation_pfams <- c('PF13452', 'PF00485', 'PF00016', 'PF00101', 'PF06240', 'PF07690' , 'PF00194', 'PF00126')
carbon_fixation_pfams <- c('PF00485', 'PF00016', 'PF00101', 'PF00194')
iron_functions <- c( 'PF00111',  'PF13085', 'PF00258')
#other_functions <- c('PF03842', 'PF16867', 'PF00127', 'PF16983','PF11124',  'PF01384', 'PF16974')
other_functions <- c('PF01061','PF00005', 'PF00854','PF07690', 'PF02421', 'PF04145', 'PF16983', 'PF02535', 'PF03842', 'PF01384','PF16974', 'PF00909')
#temperature_pfams <- c('PF05971', 'PF01180', 'PF05351', 'PF05035', 'PF01346', 'PF03104', 'PF11999',
#                       'PF00249', 'PF03259', 'PF02209', 'PF00241', 'PF00626', 'PF01119','PF10551', 'PF01979',
#                       'PF01964', 'PF02358','PF14249', 'PF02649',
#                       'PF00313', 'PF05562', 'PF14169')
temperature_pfams <- c('PF00360','PF03952','PF01786','PF06415',
                       'PF05971', 'PF01180', 'PF05351', 'PF05035', 'PF01346', 'PF03104', 'PF11999',
                       'PF00249', 'PF03259', 'PF02209', 'PF00241', 'PF00626', 'PF01119', 'PF01979',
                       'PF01964', 'PF02358','PF14249', 'PF02649',
                       'PF00313', 'PF05562', 'PF14169')

list_all_pfams <- list(photosynthesis_pfams, nitrogen_pfams, carbon_fixation_pfams, iron_functions, other_functions, temperature_pfams)
names_pfams_list <- c('Photosynthesis', 'Nitrogen metabolism', 'Carbon fixation', 'Iron metabolism', 'Nutrient transporters',                        'Temperature sensitive functions')
cols_pfams <- c('green', 'purple', 'firebrick1', 'goldenrod1','darkorange1','deepskyblue')
if (frac =='GGZZ'){
  stats <- c(144, 145,  147, 148 ,
             151, 152, 155,158,163, 150,
             168 ,173, 175, 178,180, 188, 189,
             193, 194, 196)
} else if (frac=='QQSS'){
  stats <- c(143, 144,145,146, 147,148, 149, 150 ,
            151 ,152 ,155, 158   ,168 ,173, 175, 178,180, 188, 189, 191,
            193, 194, 196)
}

dat <- readRDS(paste('metat_analysis_',frac,'.rds',sep= ''))
dat <- dat[dat[,111]>4,]
unigenes <- readLines(paste('unigenes_phate_',frac,'_',zeros,'_',norm,'_',sub2,'.txt', sep=''))
unigenes <- strsplit(unigenes, ' ')[[1]]
dat <- dat[dat[,1] %in% unigenes, ]
if (subs=='_subset'){
  new_data <- matrix(NA, ncol=length(stats), nrow=dim(dat)[1])
  for (j in 1:dim(dat)[1]){
    stts <- as.numeric(strsplit(as.character(dat$V110[j]), '_')[[1]])
    exprs <- as.numeric(strsplit(as.character(dat$V112[j]), '_')[[1]])
    new_data[j, match(stts, stats)] <- exprs
  }
 
  sum_func <- function(x){
    sum(!is.na(x))
  }
  
  sums <- apply(new_data, 1, sum_func)
  #saveRDS(sums, paste('n_stats_unigenes_',frac,'.rds'))
  new_data <- new_data[sums>12,]
  dat <- dat[sums>12,]
}

data <- read.table(paste('phate_fit_expression_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2 ,'.txt', sep=''))


data_uni_tax <- readRDS(paste('data_uni_',taxo,'.rds', sep=''))
colnames(data_uni_tax)<-c('taxName', 'geneID')
data_uni_tax<- data_uni_tax[!duplicated(data_uni_tax$geneID),]
data_uni_tax <- data_uni_tax[data_uni_tax$geneID %in% unigenes,]
if (is.factor(data_uni_tax$taxName)){
  data_uni_tax$taxName <- as.character(levels(data_uni_tax$taxName))[data_uni_tax$taxName]
}

cols <- readRDS(paste('color_table_',taxo,'.rds', sep=''))
if (is.factor(cols$col)){
 cols$col <- as.character(levels(cols$col))[cols$col]
}
if (is.factor(cols$taxon)){
 cols$taxon <- as.character(levels(cols$taxon))[cols$taxon]
}

a= min(data[,1])
b =max(data[,1])
c= min(data[,2])
d =max(data[,2])
pdf(paste('phate_fit_expression_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''))

taxos_vec  <- cols$taxon[match(data_uni_tax$taxName[match(dat[,1], data_uni_tax$geneID)],cols$taxon)]
taxos_vec[is.na(taxos_vec)] <- 'unknown'
#taxos_vec <- as.character(levels(taxos_vec))[taxos_vec]
colos <- cols$col[match(data_uni_tax$taxName[match(dat[,1], data_uni_tax$geneID)],cols$taxon)]
colos[is.na(colos)] <- cols$col[cols$taxon=='unknown']
if (is.factor(colos)){
  colos <- as.character(levels(colos))[colos]
}

c_no_pfam <- scales::alpha('black', 0.4)
col_not <- c_no_pfam
colos_bis <- rep(col_not, length(dat[,1]))
h=1
for (pfam_types in list_all_pfams){
  for (pf in pfam_types){
    colos_bis[grep(pf, dat[,49])] <-cols_pfams[h]
  }
  h=h+1
}

plot(data[,1], data[,2], xlab='PHATE 1', ylab='PHATE 2', col=scales::alpha(colos, 0.5), cex=0.5, pch=19, cex.axis=1.3, cex.lab=1.3)
dev.off()




if (taxo=='groups3'){
  sub_taxs <- c('Hexanauplia', 'Mamiellales', 'Phaeocystales', 'Pelagophyceae', 'Bacillariophyta')
  if (sub2=='BMPP'){
    sub_taxs <- c('Mamiellales', 'Phaeocystales', 'Pelagophyceae', 'Bacillariophyta')
  }
} else if (taxo=='MGT-v2'){
  sub_taxs <- c('3_Phaeocystales','439_Bathycoccaceae', '79_Bathycoccaceae', '195_Cymatosiraceae', '76_root', '4_Haptista',
                '42_Pelagomonadales')
}

sub_taxs_cols <-cols$col[match(sub_taxs,cols$taxon )] 
taxos_vec_bis <- taxos_vec
#taxos_vec_bis[!(taxos_vec_bis %in% sub_taxs)]
colos_ter <- colos
colos_ter[!(taxos_vec_bis %in% sub_taxs)] <- 'black'

if (taxo=='groups3'){
  selection <- colos_ter != 'black' & colos_ter != cols$col[cols$taxon=='Hexanauplia'] & colos_bis != c_no_pfam
} else if (taxo=='MGT-v2'){
  selection <- colos_ter != 'black' & colos_bis != c_no_pfam
}
data_sub <- data[selection,]
uni_sub <- unigenes[selection]
pfams <- dat[selection,37]
pfams_all <- dat[selection,49]
taxos_sub <- taxos_vec_bis[selection]
annot <- colos_bis[selection]
#print(annot)
annot <- sapply(annot, function(x){names_pfams_list[which(cols_pfams==x)]})
sig_cor <- readRDS('Significant_uids_all_GGZZ_0_05_clr.rds')
sig_cor$vari <- as.character(levels(sig_cor$vari))[sig_cor$vari]
sig_cor <- sig_cor[!(sig_cor$vari %in% c('Metagenome', 'Environment', 'Travel_time')),]
sig_cor$cor_s <- sign(sig_cor$cor)
sig_cor$cor_s[sig_cor$cor_s==-1] = '-'
sig_cor$cor_s[sig_cor$cor==1] = '+'
sig_cor <- sig_cor[order(abs(sig_cor$cor), decreasing=T),]
sig_cor <- sig_cor[!(duplicated(sig_cor$uni)),]
sig_cor$vari_sign<- paste(sig_cor$vari, sign(sig_cor$cor), sep='_')

to_save<- data.frame('uid'=uni_sub, 'pfam'=pfams,'pfams'=pfams_all, 'taxo'=taxos_sub, 'annot'=annot)
to_save$vari_sign <- NA
to_save$vari_sign <- sig_cor$vari_sign[match(to_save$uid, sig_cor$uni)]
saveRDS(to_save, paste('sub_phate_pfam_taxo_subtax_', knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.rds', sep='')) 
saveRDS(data_sub, paste('sub_phate_fit_expression_py_subtax_', knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.rds', sep=''))

pdf(paste('phate_fit_expression_py_subtax_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''))
plot(0, 0, xlab='PHATE 1', ylab='PHATE 2',xlim=c(a,b), ylim=c(c,d), col=scales::alpha('white', 0), cex=0.5, pch=19, cex.axis=1.3, cex.lab=1.3)
if (length(data[colos_ter=='black',1]) != 0){
  points(data[colos_ter=='black',1], data[colos_ter=='black',2], col=scales::alpha('black', 0.05), cex=0.5, pch=19)
}
xa=data[colos_ter!='black',1]
ya=data[colos_ter!='black',2]
cls=scales::alpha(colos_ter[colos_ter!='black'])
#print(length(xa))
set.seed(1)
shuf=sample.int(length(xa))
points(xa[shuf], ya[shuf], col=scales::alpha(cls[shuf], 0.5)) 
#for (co in sub_taxs_cols){
#  points(data[colos_ter==co,1], data[colos_ter==co,2], col=scales::alpha(colos_ter[colos_ter==co], 0.5) )
#}
dev.off()

if (taxo=='groups3'){
  cols_low_d <- c('yellow2', 'moccasin','plum1','steelblue1','lightcoral')
} else if (taxo=='MGT-v2'){
  cols_low_d <- rep("grey", length(sub_taxs))
}
my.cols_tax <- rep(list(NULL) , length(sub_taxs))
for (i in 1:length(sub_taxs_cols)){
  my.cols_tax[[i]] <- colorRampPalette(c(cols_low_d[i], sub_taxs_cols[i]))(8)
}
pdf(paste('phate_fit_expression_py_subtax_contours_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''))
plot(0, 0,xlim=c(a,b),ylim=c(a,b), xlab='PHATE 1', ylab='PHATE 2', col=scales::alpha('white', 0), cex=0.5, pch=19, cex.axis=1.3, cex.lab=1.3)
k=1
for (co in sub_taxs_cols){
  x=data[colos_ter==co,1]
  y=data[colos_ter==co,2]
  print(length(x))
  z <- kde2d(x, y, n=100)
  contour(z, drawlabels=FALSE, nlevels=5,lwd=3, col=my.cols_tax[[k]], add=TRUE)
  k=k+1
}
dev.off()

data_cont=list()
no_co=scales::alpha('white', 0)
pdf(paste('phate_fit_expression_py_taxos_contours_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''))
plot(0, 0,xlim=c(a,b),ylim=c(a,b), xlab='PHATE 1', ylab='PHATE 2', col=scales::alpha('white', 0), cex=0.5, pch=19, cex.axis=1.3, cex.lab=1.3)
k=1
for (co in unique(colos)){
  x=data[colos==co,1]
  y=data[colos==co,2]
  #print(length(x))
  z <- kde2d(x, y, n=1000)
  data_cont[[co]] = z
  contour(z, drawlabels=FALSE, nlevels=8,lwd=3, col= c(rep(no_co, 6), co), add=TRUE)
  
  k=k+1
}
data0 =cbind(data,colos) 
to_save_cont=list(a,b,unique(colos), data_cont)
saveRDS(to_save_cont, 'data-phate-contour_plot.rds')
#df=data.frame(x=data[,1], y=data[,2], z=colos)
#find_hull <- function(df) df[chull(df$x, df$y), ]
#hulls <- ddply(df, "z", find_hull)
#ploto <- ggplot(data = dot, aes(x = x, y = y)) + 
#geom_polygon(data = hulls, alpha = 0.5) +
#labs(x = "PHATE 1", y = "PHATE 2") +theme_bw(base_size = 23)
#plot(0, 0,xlim=c(a,b),ylim=c(a,b), xlab='PHATE 1', ylab='PHATE 2', col=scales::alpha('white', 0), cex=0.5, pch=19, cex.axis=1.3, cex.lab=1.3)
#Hulls(pts = df[,1:2], df[, 3], centers=T, outliers = F, c.pch=19, usecolors=unique(colos))
dev.off()


pdf(paste('legend_contours_subtax_',taxo,'.pdf', sep='')) 
pnt <- cbind(x =c(0,5,5,0), y =c(0,50,50,0))
plot(0,0, col='white', xlim=c(0,30.5*length(sub_taxs_cols)), ylim=c(0, 60), axes=FALSE, frame.plot=F, xlab = '', ylab='')
for (k in 1:length(my.cols_tax)){
  pnt0 <- pnt
  pnt0[,1] <- pnt0[,1]+(k-1)*30 
  SDMTools::legend.gradient(pnt0, my.cols_tax[[k]], title =sub_taxs[k],limits=c('',''), cex=1)
}
dev.off()


txs=NULL
sums=NULL
pdf(paste('phate_fit_expression_decomposed_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''), width=10, height=10)
par(mfrow=c(2,2))
for (tx in unique(cols$taxon)){
  selo <- taxos_vec==tx
  if (sum(selo)>200){
    plot(data[selo,1], data[selo,2], xlab='PHATE 1',main=tx, ylab='PHATE 2', col=scales::alpha(colos[selo], 1), cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))
    txs=c(txs, tx)
    sums=c(sums, sum(selo))
  }
}
dev.off()

ordinat = order(sums,decreasing=T )
txs=txs[ordinat]

lists=rep(list(NULL), length(txs))
for (j in 1:length(txs)){
  lists[[j]] = txs[1:j]
}
print('ok')
print(lists)
for (j in 1:length(txs)){
  li_tx=lists[[j]]
  if (j<10){
    png(paste('phate_fit_expression_decomposed_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'_0',j,'.png', sep=''), width=4.5,height=4.5,units="in",res=2000)
  } else {
    png(paste('phate_fit_expression_decomposed_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'_',j,'.png', sep=''), width=4.5,height=4.5,units="in",res=2000)
  }
  cou=1
  for (tx in li_tx){
    selo <- taxos_vec==tx
    if (cou==1){
      m_tx=li_tx[length(li_tx)]
      plot(data[selo,1], data[selo,2], xlab='PHATE 1',main=m_tx, ylab='PHATE 2', col=scales::alpha(colos[selo], 1), cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))
    } else{
      points(data[selo,1], data[selo,2] ,  col=scales::alpha(colos[selo], 1), cex=0.2, pch=19)
    }
    cou=cou+1
  }
  print(j)
  dev.off()
}


for (tx in unique(cols$taxon)){
  selo <- taxos_vec==tx
  if (sum(selo)>200){
    png(paste('phate_fit_expression_decomposed_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'_',tx,'.png', sep='') ,width=4.5,height=4.5,units="in",res=2000)
    plot(data[selo,1], data[selo,2], xlab='PHATE 1',main=tx, ylab='PHATE 2', col=scales::alpha(colos[selo], 1), cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))
    dev.off()
  }
}


sel_bis <- colos_bis==col_not
pdf(paste('phate_fit_expression_pfams_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''), width=5, height=5)
plot(data[sel_bis,1], data[sel_bis,2], xlab='PHATE 1', ylab='PHATE 2', col=colos_bis[sel_bis], cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))
points(data[!sel_bis ,1], data[!sel_bis,2], col=colos_bis[!sel_bis], cex=0.7, pch=16)
dev.off()

pdf(paste('phate_fit_expression_pfams_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'_T.pdf', sep=''), width=5, height=5)
par(mfrow=c(1,1))
comp <- 'T'
signif_uids <- readRDS(paste('Significant_uids_all_',frac,'_0_05_',comp,'_clr.rds', sep=''))
signif_uids$cor_sign <- as.numeric(signif_uids$cor>0)
signif_uids$cor_sign[signif_uids$cor_sign==0] <-'-'
signif_uids$cor_sign[signif_uids$cor_sign==1] <-'+'
signif_uids_moins <- signif_uids$uni[signif_uids$cor_sign=='-']
signif_uids_plus <- signif_uids$uni[signif_uids$cor_sign=='+']
sel_m <- dat[,1] %in%  signif_uids_moins
sel_p <- dat[,1] %in%  signif_uids_plus

plot(data[sel_bis,1], data[sel_bis,2], xlab='PHATE 1', ylab='PHATE 2', col=colos_bis[sel_bis], cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3 ,xlim=c(a, b), ylim=c(c, d))

points(data[sel_m &  !sel_bis ,1], data[sel_m & !sel_bis,2], col=colos_bis[sel_m & !sel_bis], cex=0.6)
points(data[sel_p &  !sel_bis ,1], data[sel_p & !sel_bis,2], col=colos_bis[sel_p & !sel_bis], cex=0.6, pch=6)
dev.off()

pdf(paste('phate_fit_expression_decomposed_pfams_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''), width=5, height=5)
par(mfrow=c(1,1))
for (tx in unique(cols$taxon)){
  selo <- taxos_vec==tx
  sel_bis <- colos_bis==col_not
  if (sum(selo)>200){
    plot(data[selo & sel_bis,1], data[selo & sel_bis,2], xlab='PHATE 1',main=tx, ylab='PHATE 2', col= scales::alpha(colos_bis[selo & sel_bis], 0.2), cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))
    points(data[selo & !sel_bis,1], data[selo & !sel_bis,2],pch=16, col=colos_bis[selo & !sel_bis], cex=0.7)
  }
}
dev.off()

# c('green', 'purple', 'firebrick1', 'goldenrod1','darkorange1','deepskyblue')
my.cols=list(colorRampPalette(brewer.pal(9,"Greens"))(8), colorRampPalette(brewer.pal(9,"Purples"))(8),
             colorRampPalette(brewer.pal(9,"Reds"))(8), colorRampPalette(brewer.pal(9,"YlOrBr"))(8),
             colorRampPalette(brewer.pal(9,"Oranges"))(8), colorRampPalette(brewer.pal(9,"PuBu"))(8) )

pdf(paste('phate_fit_density_metabolisms_decomposed_pfams_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''), width=5, height=5)
par(mfrow=c(1,1))
if (taxo =='groups3'){
  txs <- unique(cols$taxon)
} else if (taxo =='MGT-v2'){
  txs <- sub_taxs
}
for (tx in txs){
  selo <- taxos_vec==tx
  sel_bis <- colos_bis==col_not
  contours_list <- rep(list(NULL), 8)
  centr_list <- rep(list(NULL), 8)
  for (k in 1:length(cols_pfams)){
    func <- names_pfams_list[k]
    sel_pf <- colos_bis==cols_pfams[k]
    if (sum(selo & sel_pf)>4){
      xu = data[selo & sel_pf & !sel_bis,1]
      yu = data[selo & sel_pf & !sel_bis,2]
      z <- kde2d(xu, yu, n=100)
      contours_list[[k]] <- z
      centr_list[[k]] <- append(centr_list[[k]], c(median(xu), median(yu)))
      plot(xu, yu, xlab='PHATE 1',main= paste(tx,func , sep=' '), ylab='PHATE 2', col=cols_pfams[k], cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))
      contour(z, drawlabels=FALSE, nlevels=8, col=my.cols[[k]], add=TRUE)
    }
  }
  plot(0, 0, xlab='PHATE 1',main= tx, ylab='PHATE 2', col=scales::alpha('white', 0), cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3 ,xlim=c(a, b), ylim=c(c, d))
  for (k in 1:length(cols_pfams)){
    if (!is.null(contours_list[[k]])){
      #points(centr_list[[k]][1] , centr_list[[k]][2], col=cols_pfams[k], pch=19, cex=0.8)
      contour(contours_list[[k]], drawlabels=FALSE,lwd=1.8, nlevels=4, col=scales::alpha(my.cols[[k]], 0.6), add=TRUE)
    }
  }
  for (k in 1:length(cols_pfams)){
    if (!is.null(contours_list[[k]])){
      points(centr_list[[k]][1] , centr_list[[k]][2], col=cols_pfams[k], pch=19, cex=1)
      #contour(contours_list[[k]], drawlabels=FALSE,lwd=2, nlevels=2, col=my.cols[[k]], add=TRUE)
    }
  }
}
dev.off()


pdf(paste('phate_fit_density_metabolisms_pfams_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'.pdf', sep=''), width=5, height=5)
par(mfrow=c(1,1))
contours_list <- rep(list(NULL), 8)
centr_list <- rep(list(NULL), 8)
for (k in 1:length(cols_pfams)){
  func <- names_pfams_list[k]
  sel_pf <- colos_bis==cols_pfams[k] 
  if (sum(sel_pf)>4){
    sel_pf <- colos_bis==cols_pfams[k]
    xu = data[sel_pf & !sel_bis,1]
    yu = data[sel_pf & !sel_bis,2]
    z <- kde2d(xu, yu, n=100)
    contours_list[[k]] <- z
    centr_list[[k]] <- append(centr_list[[k]], c(median(xu), median(yu)))
    plot(xu, yu, xlab='PHATE 1',main= func, ylab='PHATE 2', col=cols_pfams[k], cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))
    contour(z, drawlabels=FALSE,lwd=2, nlevels=8, col=my.cols[[k]], add=TRUE)
  }
}
plot(0, 0, xlab='PHATE 1',main= 'all', ylab='PHATE 2', col=scales::alpha('white', 0), cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))
for (k in 1:length(cols_pfams)){
  if (!is.null(contours_list[[k]])){
    #points(centr_list[[k]][1] , centr_list[[k]][2], col=cols_pfams[k], pch=19, cex=0.8)    
    contour(contours_list[[k]], drawlabels=FALSE,lwd=1.4, nlevels=4, col=scales::alpha(my.cols[[k]], 0.75), add=TRUE)
  }
}
for (k in 1:length(cols_pfams)){
  if (!is.null(contours_list[[k]])){
    points(centr_list[[k]][1] , centr_list[[k]][2], col=cols_pfams[k], pch=19, cex=1)
    #contour(contours_list[[k]], drawlabels=FALSE,lwd=2, nlevels=2, col=my.cols[[k]], add=TRUE)
  }
}
dev.off()

pdf(paste('legend_contours_metabolisms.pdf', sep=''))
pnt <- cbind(x =c(0,5,5,0), y =c(0,50,50,0))
plot(0,0, col='white', xlim=c(0,30.5*length(my.cols)), ylim=c(0, 60), axes=FALSE, frame.plot=F, xlab = '', ylab='')
for (k in 1:length(my.cols)){
  pnt0 <- pnt
  pnt0[,1] <- pnt0[,1]+(k-1)*30
  SDMTools::legend.gradient(pnt0, my.cols[[k]], title =names_pfams_list[k],limits=c('','') , cex=1) 
}
dev.off()

pdf(paste('phate_fit_expression_decomposed_pfams_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'_T.pdf', sep=''), width=5, height=5)
par(mfrow=c(1,1))
for (tx in unique(cols$taxon)){
  selo <- taxos_vec==tx
  sel_bis <- colos_bis==col_not
  if (sum(selo)>200){
    plot(data[selo & sel_bis,1], data[selo & sel_bis,2], xlab='PHATE 1',main=tx, ylab='PHATE 2', col=scales::alpha(colos_bis[selo & sel_bis], 0.2), cex=0.2, pch=19, cex.axis=1.3,cex.lab=1.3, xlim=c(a, b), ylim=c(c, d))

    points(data[selo & sel_m &  !sel_bis ,1], data[selo & sel_m & !sel_bis,2], col=colos_bis[selo & sel_m & !sel_bis], cex=0.7)
    points(data[selo & sel_p &  !sel_bis ,1], data[selo & sel_p & !sel_bis,2], col=colos_bis[selo & sel_p & !sel_bis], cex=0.7, pch=6)
    #plot(data[selo & sel_bis,1], data[selo & sel_bis,2], xlab='PHATE1',main=tx, ylab='PHATE2', col=colos_bis[selo & sel_bis], cex=0.2, pch=19, cex.axis=2, xlim=c(a, b), ylim=c(c, d))
    #points(data[selo & !sel_bis,1], data[selo & !sel_bis,2], col=colos_bis[selo & !sel_bis], cex=0.6)
  }
}
dev.off()

#comps <- c('T', 'Sal', 'SSD','O2' ,'NO3', 'Phos','Si', 'NO2', 'NH4', 'N.', 'Si.', 'Fe')
comps <- c('T')
for (comp in comps){  
  signif_uids <- readRDS(paste('Significant_uids_all_',frac,'_0_05_',comp,'_clr.rds', sep=''))
  signif_uids$cor_sign <- as.numeric(signif_uids$cor>0)
  signif_uids$cor_sign[signif_uids$cor_sign==0] <-'-'
  signif_uids$cor_sign[signif_uids$cor_sign==1] <-'+'
  signif_uids_moins <- signif_uids$uni[signif_uids$cor_sign=='-']
  signif_uids_plus <- signif_uids$uni[signif_uids$cor_sign=='+']
  sel_m <- dat[,1] %in%  signif_uids_moins
  sel_p <- dat[,1] %in%  signif_uids_plus
  x <- data[!sel_m & !sel_p,1]
  y <- data[!sel_m & !sel_p,2]
  if (length(x) ==0){
    x <-0
    y <- 0
    alph <- 0
  } else{
    alph <- 0.1
  } 
  pdf(paste('phate_fit_expression_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'_',comp,'.pdf', sep=''))
  plot(x, y, xlab='PHATE 1',xlim=c(min(data[,1]), max(data[,1])),ylim=c(min(data[,2]), max(data[,2])) , ylab='PHATE 2', col=scales::alpha('black', alph), cex=0.5, pch=19, cex.axis=1.3,cex.lab=1.3, main=comp,
       cex.main=2, cex.lab=2)
  points(data[sel_m,1], data[sel_m,2], col=scales::alpha('deepskyblue', 0.7), pch=15, cex=0.6)
  points(data[sel_p,1], data[sel_p,2], col=scales::alpha('red', 0.7), pch=17, cex=0.6)
  #plot(data[!sel_p,1], data[!sel_p,2], xlab='PHATE1', ylab='PHATE2', col=scales::alpha('black', 0.1), cex=0.5, pch=19, cex.axis=2, main=paste(comp,'_+', sep=''),
  #     cex.main=2, cex.lab=2)
  #points(data[sel_p,1], data[sel_p,2], col=colos[sel_m], pch=19)
  dev.off()
  pdf(paste('phate_fit_expression_decomposed_py_',knn,'_',frac,'_',zeros,'_',norm,'_',sub2,'_',taxo,'_',comp,'.pdf', sep=''), width=10, height=10)
  par(mfrow=c(2,2))
  for (tx in unique(cols$taxon)){
    selo <- taxos_vec==tx
    if (sum(selo & (sel_m | sel_p))>200){
      plot(data[selo & sel_m,1], data[selo & sel_m,2], xlab='PHATE 1',main=paste(tx, ' ', comp, sep=''), ylab='PHATE 2', col='deepskyblue', cex=0.6, pch=15, cex.axis=1.3,cex.lab=1.3,xlim=c(a, b), ylim=c(c, d))
      points(data[selo & sel_p,1], data[selo & sel_p,2], pch=17, col= 'red')
      #plot(data[selo & sel_p,1], data[selo & sel_p,2], xlab='PHATE1',main=paste(tx, ' ', comp, '_p', sep=''), ylab='PHATE2', col=scales::alpha(colos[selo & sel_p], 1), cex=0.6, pch=17, cex.axis=2, xlim=c(a, b), ylim=c(c, d))
    }
    #points(data[selo & sel_p,1], data[ selo & sel_p,2], col=scales::alpha(colos[selo & sel_p], 1), pch=17, cex=0.2)
  }
  dev.off()
  #pdf(paste('phate_fit_expression_decomposed_',frac,'_',nd,'_',knn,'_',di,'_',taxo,'_',comp,'_-.pdf', sep=''), width=15, height=15)
  #par(mfrow=c(3,3))
  #for (tx in unique(cols$taxon)){
  #  selo <- colos==cols$col[cols$taxon==tx]
  #  plot(data[selo & sel_p,1], data[selo & sel_p,2], xlab='PHATE1',main=tx, ylab='PHATE2', col=scales::alpha(colos[selo & sel_p], 1), cex=0.4, pch=17, cex.axis=2,
  #       xlim=c(a, b), ylim=c(c, d))
    #points(data[selo & sel_p,1], data[ selo & sel_p,2], col=scales::alpha(colos[selo & sel_p], 1), pch=17, cex=0.2)
  #}
  #dev.off()
}



