library('phateR')
frac <- commandArgs(trailingOnly = T)[1]
subs <- commandArgs(trailingOnly = T)[2]
zeros <- commandArgs(trailingOnly = T)[3]
norm <- commandArgs(trailingOnly = T)[4]
sub2 <- commandArgs(trailingOnly = T)[5]
dat <- readRDS(paste('metat_analysis_',frac,'.rds',sep= ''))
dat <- dat[dat[,111]>4,]
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
dat <- dat[,c(1,49, 110, 112)]
new_data <- matrix(NA, ncol=length(stats), nrow=dim(dat)[1])
unigenes_ids <- dat$V1
for (j in 1:dim(dat)[1]){
  stts <- as.numeric(strsplit(as.character(dat$V110[j]), '_')[[1]])
  exprs <- as.numeric(strsplit(as.character(dat$V112[j]), '_')[[1]])
  new_data[j, match(stts, stats)] <- exprs
}
sum_func <- function(x){
  sum(!is.na(x))
}
if (subs=='_subset'){
  sums <- apply(new_data, 1, sum_func)
  saveRDS(sums, paste('n_stats_unigenes_',frac,'.rds'))
  new_data <- new_data[sums>12,]
  unigenes_ids <- unigenes_ids[sums>12]
}
if (norm=='0') {
  new_data[is.na(new_data)] <- as.numeric(zeros)
  suffix <- ''
} else if (norm=='1'){
  new_max <- 1
  new_min <- 0.1
  for (j in 1:dim(new_data)[1]){
    dt <- new_data[j,]
    mi_dt <- min(dt, na.rm=T)
    mx_dt <- max(dt, na.rm=T)
    new_data[j,] <- (dt - mi_dt)*(new_max-new_min)/(mx_dt-mi_dt)+new_min
    #print(new_data[j,])
  }
  new_data[is.na(new_data)] <- as.numeric(zeros)
  suffix <- '_norm'
} else if (norm=='2'){ 
  mi_dt <- min(new_data, na.rm=T)
  mx_dt <- max(new_data, na.rm=T)
  new_min <- 0.1
  new_max <- (mx_dt-mi_dt)+new_min
  new_data <- (new_data - mi_dt)*(new_max-new_min)/(mx_dt-mi_dt)+new_min
  new_data[is.na(new_data)] <- as.numeric(zeros)
}


if (sub2 != '0'){
  taxo <- readRDS('data_uni_groups3.rds')
  if (sub2=='HBMPP'){
    taxo <- taxo[taxo$taxName %in% c('Hexanauplia', 'Bacillariophyta','Mamiellales', 'Pelagophyceae', 'Phaeocystales'),]
    sel <- unigenes_ids %in% taxo$geneID
    unigenes_ids <- unigenes_ids[sel]
    new_data <- new_data[sel,]
  } else if (sub2=='BMPP'){
    taxo <- taxo[taxo$taxName %in% c('Bacillariophyta','Mamiellales', 'Pelagophyceae', 'Phaeocystales'),]
    sel <- unigenes_ids %in% taxo$geneID
    unigenes_ids <- unigenes_ids[sel]
    new_data <- new_data[sel,]
  } else if (sub2=='HBMPP+-'){
    taxo <- taxo[taxo$taxName %in% c('Hexanauplia', 'Bacillariophyta','Mamiellales', 'Pelagophyceae', 'Phaeocystales'),]
    sel <- unigenes_ids %in% taxo$geneID
    unigenes_ids <- unigenes_ids[sel]
    new_data <- new_data[sel,]
    sigs <- readRDS('Significant_uids_all_GGZZ_0_05_T_clr.rds')
    sel0 <- unigenes_ids %in% sigs$uni
    unigenes_ids <- unigenes_ids[sel0]
    new_data <- new_data[sel0,]
  } else if (sub2=='BMPP+-'){
    taxo <- taxo[taxo$taxName %in% c('Bacillariophyta','Mamiellales', 'Pelagophyceae', 'Phaeocystales'),]
    sel <- unigenes_ids %in% taxo$geneID
    unigenes_ids <- unigenes_ids[sel]
    new_data <- new_data[sel,]
    sigs <- readRDS('Significant_uids_all_GGZZ_0_05_T_clr.rds')
    sel0 <- unigenes_ids %in% sigs$uni
    unigenes_ids <- unigenes_ids[sel0]
    new_data <- new_data[sel0,]
  } else if (sub2=='all+-'){
    sigs <- readRDS('Significant_uids_all_GGZZ_0_05_T_clr.rds')
    sel0 <- unigenes_ids %in% sigs$uni
    unigenes_ids <- unigenes_ids[sel0]
    new_data <- new_data[sel0,]
  }
  
}
pdf(paste('expr_distirbution_phate_',frac,'_',zeros,subs,'_',norm,'_',sub2,'.pdf', sep=''))
hist(new_data, breaks=30)
dev.off()

writeLines(as.character(stats), paste('stations_',frac,'.txt', sep=''), sep=' ')

if (subs=='0'){
  subs = ''
}
writeLines(as.character(unigenes_ids), paste('unigenes_phate_',frac,'_',zeros,subs,'_',norm,'_',sub2,'.txt', sep=''), sep=' ')
#}
write.table(new_data,paste('phate_training_expression_',frac,'_',zeros,subs,'_',norm,'_',sub2,'.txt', sep=''),
	   row.names = F,
            col.names = F )
set.seed(1)
phate_fit <- phate(new_data, ndim=nd, knn=knn, knn.dist.method = di)
saveRDS(phate_fit, paste('phate_fit_expression_',frac,'_',nd,'_',knn,'_',di,'_',zeros,subs,norm,'.rds', sep=''))
