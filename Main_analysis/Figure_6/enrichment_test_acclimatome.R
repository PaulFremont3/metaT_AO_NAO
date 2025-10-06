library(UpSetR)
library(tidyverse)
library(grid)
library(gplots)
source('vioplot1.R')
GO_table=readRDS('../data/GO_table.rds')
table=read.table('../data/training_set_phate_GO_all_GGZZ_physical_clr_T.txt')
GOs=readLines('../data/training_set_phate_GO_rows_GGZZ_physical_clr_T.txt')
GOs=strsplit(GOs, split = "\t")[[1]]
table=table[!(GOs %in% c('unknown', 'NO_GO')),]
GOs=GOs[!(GOs %in% c('unknown', 'NO_GO'))]
sub_GO_table=GO_table[!is.na(match(GO_table[,2], GOs)),]
sub_GO_table=sub_GO_table[match(GOs, sub_GO_table[,2]),]



baseline=readRDS('../data/baseline_GO_matou_groups3s_GGZZ_subset.rds')
baseline=baseline[!(baseline$Group.2 %in% c('unknown', 'NO_GO')),]



taxos=c('Bacillariophyta', 'Mamiellales', 'Pelagophyceae', 'Phaeocystales')
taxos_bis =c(taxos, 'Hexanauplia')
cols= c('Bacillariophyta T_-',	'Bacillariophyta T_+',
        'Hexanauplia T_-'	,'Hexanauplia T_+',	'Mamiellales T_-',	'Mamiellales T_+'	,
        'Pelagophyceae T_-',	'Pelagophyceae T_+'	,'Phaeocystales T_-'	,'Phaeocystales T_+')

ord=c(1:2, 5:10, 3:4)
colnames(table)=cols
table=table[,ord]


sels_ter=list()
c=1
for (i in taxos_bis){
  for (j in taxos_bis){
    ind_i=which(taxos_bis==i)
    ind_j=which(taxos_bis==j)
    if (ind_i>ind_j){
      sel1=which( (grepl(i, colnames(table)) | grepl(j, colnames(table)) == T ) & !grepl('T_-', colnames(table)))
    } else if (ind_i<ind_j){
      sel1=which( (grepl(i, colnames(table)) | grepl(j, colnames(table)) == T ) & grepl('T_-', colnames(table)))
    } else if (ind_i==ind_j){
      sel1=which(grepl(i, colnames(table)) | grepl(j, colnames(table)) == T)
    }
    sels_ter[[c]]=sel1
    c=c+1
  }
}

pvals=NULL
par(mfrow = c(5, 5))
c=1
alg_m_sel = c(2,3,4,8,9, 14)
alg_m = NULL
alg_p_sel = c(6, 11,12, 16,17,18)
alg_p= NULL
hex_m_sel = c(5, 10, 15, 20)
hex_m = NULL
hex_p_sel= c(21, 22,23,24)
hex_p = NULL
alg_mp =NULL
alg_mp_sel =c(1,7,13,19, 25)
all_data =NULL
all_data_lab=NULL
for (sel in sels_ter){
  t1=str_split(colnames(table)[sel[1]], pattern = ' T')[[1]][1]
  t2=str_split(colnames(table)[sel[2]], pattern = ' T')[[1]][1]
  
  base1=baseline[baseline$Group.1==t1, ]
  base2=baseline[baseline$Group.1==t2, ]
  baseall=rbind(base1, base2)
  m1=match(base1$Group.2, base2$Group.2)
  base1_s=base1[!is.na(m1),]
  m2=match(base2$Group.2, base1$Group.2)
  base2_s=base2[!is.na(m2),]
  base_s_table=cbind(base2_s$x, base1_s$x)
  min_s=apply(base_s_table, 1, min)
  NS=sum(min_s)
  #print(NS/(N1+N2))
  
  N1=sum(base1$x)
  N2=sum(base2$x)
  
  table1=table[, sel]
  min_s_exp=apply(table1, 1, min)
  Es=sum(min_s_exp)
  
  E1=sum(table[,sel[1]])
  E2=sum(table[,sel[2]])
  
  print(colnames(table)[sel[1]])
  print(colnames(table)[sel[2]])
  print(NS/(N1+N2-NS))
  print(Es/(E1+E2-Es))
  print(' ')
  
  q=Es
  m=NS
  n=N1+N2-2*NS
  k=E1+E2-2*Es
  
  if (c %in% alg_m_sel){
    alg_m = c(alg_m, 100*Es/(E1+E2-Es))
    all_data_lab =c(all_data_lab,1)
  } else if (c %in% alg_p_sel){
    alg_p = c(alg_p, 100*Es/(E1+E2-Es))
    all_data_lab =c(all_data_lab,2)
  } else if (c %in% hex_m_sel){
    hex_m = c(hex_m, 100*Es/(E1+E2-Es))
    all_data_lab =c(all_data_lab,3)
  } else if (c %in% hex_p_sel){
    hex_p = c(hex_p, 100*Es/(E1+E2-Es))
    all_data_lab =c(all_data_lab,4)
  } else if (c %in% alg_mp_sel){
    alg_mp = c(alg_mp, 100*Es/(E1+E2-Es))
    all_data_lab =c(all_data_lab,5)
  }
  all_data =c(all_data, 100*Es/(E1+E2-Es))
  pval=phyper(q,m,n,k, lower.tail = F)
  
  par(cex = 5)
  r1=colnames(table)[sel[1]]
  r2=colnames(table)[sel[2]]
  pvals =c(pvals, pval)
  c=c+1
}
dev.off()
pvals_ad=p.adjust(pvals, 'hommel')
pvals_ad[pvals_ad>0.05]=1
sig=pvals_ad

pvalsm=matrix(sig, nrow=5, byrow = T)
#heatmap.2(pvalsm, dendrogram = 'none', Rowv = NULL, Colv = NULL, trace = 'none', labRow = taxos_bis,
#          labCol = taxos_bis)


pdf('violin_plot_acclimatome.pdf', width=9)
par(mar=c(10.1, 4.1, 4.1, 2.1))
vioplot1(list(alg_m, alg_p, alg_mp_sel), col = rep('darkgray',5),
         names=c('Algae vs Algae T-','Algae vs Algae T+',
                 'T- vs T+'), ylab = '% common genes')
stripchart(list(alg_m, alg_p,  alg_mp_sel),cex=2, vertical = TRUE, method = "jitter",
           pch = 19, add = TRUE, col = rep('orange', 5))

dev.off()
