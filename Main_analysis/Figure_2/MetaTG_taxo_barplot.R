library('dplyr')
library('stringr')

frac = commandArgs(trailingOnly = T)[1] #'GGZZ'
type = commandArgs(trailingOnly = T)[2] # c('T', 'uni_T', 'G', 'uni_G')
tax_id = commandArgs(trailingOnly = T)[3] # 'txo_groups3'
subset = commandArgs(trailingOnly = T)[4] # 0 or 1
classi = commandArgs(trailingOnly = T)[5] # 2 to 6 (classification of the unigenes (unique atlantic/arctic etc)
bis = commandArgs(trailingOnly = T)[6] #'_bis' or 0

if (bis!='_bis'){
  bis=''
}

if (tax_id=='taxo_MGT-v2'){
  colors_taxon <- readRDS('../data/color_table_MGT-v2.rds')
  os=15
} else if (tax_id=='taxo_groups3'){
  colors_taxon <- readRDS('../data/color_table_groups3.rds')
  os=15
}

if (type %in% c('uni_T', 'uni_G')){
  type0 <- paste(type, '1', sep='')
} else{
  type0 <- type
} 
clusts <- c(type0,paste(type0,c('com', 'atl', 'arc'), sep='_'))
clusts0 <- c('all','com', 'atl', 'arc')

c=1
list_to_plot <- list()
for (c in 1:length(clusts)){
  if (classi %in% c(5,6)) {
    data0<-readRDS(paste('../data/Meta', clusts[c],'_',tax_id,'_', frac,'_',classi,bis,'.rds', sep=''))
    if (c==2){
      h <- apply(data0, 2, sum)
      #print(h)
    }
  } else{
    data0<-readRDS(paste('../data/Meta', clusts[c],'_',tax_id,'_', frac,bis,'.rds', sep=''))
  }
  list_to_plot[[c]]=data0
  c=c+1
}

if (subset=='1'){
  for (i in 1:4){
    ltp <- list_to_plot[[i]]
    if (i==1){
      d_ltp <- dim(ltp)[1]
      root_row=which(rownames(ltp)=='root')
      sel=unique(c((d_ltp-os):d_ltp, root_row))
      to_take <- rownames(ltp)[sel]
    }
    ltp <- ltp[rownames(ltp) %in% to_take,]
    list_to_plot[[i]] <- ltp
  }
}



list_to_plot_in <- list()
list_to_plot_in_norm <- list()
if (subset=='0'){
  if (classi %in%  c(5,6)){
    pdf(paste('Meta', type0,'_',tax_id,'_', frac,'_',classi,bis,'.pdf', sep=''), width = 15, height = 7)
  } else{
    pdf(paste('Meta', type0,'_',tax_id,'_', frac,bis,'.pdf', sep=''), width = 15, height = 7)
  }
} else{
  if (classi %in%  c(5,6)){
    pdf(paste('Meta', type0,'_',tax_id,'_', frac,'_',subset,'_',classi,bis,'.pdf', sep=''), width = 15, height = 7)
  } else{
    pdf(paste('Meta', type0,'_',tax_id,'_', frac,'_',subset,bis,'.pdf', sep=''), width = 15, height = 7)
  }
}
data0 <- list_to_plot[[1]]
sts <- colnames(data0)
print(sts)
sts <- sapply(sts, function(x){str_replace(x, 'SUR', 'SRF')})
sts <- as.character(sts)
colnames(data0) <- sts

data0 <- data0[,!(colnames(data0) %in% c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]

sts <- colnames(data0)
if (type=='G' & subset == '0'){
  sums <- apply(data0, 2, sum)
  data0 <- rbind(data0, 1-sums)
  list_to_plot_in[[1]] <- data0
  rownames(data0)[dim(data0)[1]] <- 'Not expressed'
  colors_taxon <- rbind(colors_taxon, c('black', 'Not expressed'))
}
if (type=='uni_G' & subset == '0'){
  uni_count <- readRDS(paste('../data/count_unis_metaG_all_',frac,'.rds', sep=''))
  uni_count <- uni_count[!(names(uni_count) %in% c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]
  data0 <- rbind(data0, uni_count)
  data0[data0<0] <- 0
  rownames(data0)[dim(data0)[1]] <- 'Not expressed'
  list_to_plot_in[[1]] <- data0
  colors_taxon <- rbind(colors_taxon, c('black', 'Not expressed'))
}
par(mfrow=c(1,2))
if (subset=='0'){
  d_ltp <- dim(data0)[1]
  if (type %in% c('G' , 'uni_G')){
    d_ltp <- d_ltp-1
    root_row=which(rownames(data0)=='root')
    sel=unique(c((d_ltp-os):d_ltp, root_row))
    dt <- apply(data0[sel,], 2, sum)
    anti_sel=1:(d_ltp+1)
    anti_sel=anti_sel[!(anti_sel %in% sel)]
    data0=data0[anti_sel,]
    data0 <- rbind(data0, dt)
  } else{
    root_row=which(rownames(data0)=='root')
    sel=unique(c((d_ltp-os):d_ltp, root_row))  
    dt <- apply(data0[sel,], 2, sum)
    anti_sel=1:d_ltp
    anti_sel=anti_sel[!(anti_sel %in% sel)]
    data0=data0[anti_sel,]
    data0 <- rbind(data0, dt)
  }
  rownames(data0)[dim(data0)[1]] <- 'Other'
  if (sum(rownames(data0)=='other')+sum(rownames(data0)=='Other')==2 ){
    inds <- c(which(rownames(data0)=='other'), which(rownames(data0)=='Other'))
    data0[inds[2],] <- data0[inds[1],] + data0[inds[2],]
    data0 <- data0[-inds[1], ]
  }
  rownames(data0)[which(rownames(data0)=='Other')] ='other'
  if (type %in% c('G' , 'uni_G')){
    di <- dim(data0)[1]
    data0 <- data0[c(1:(di-2), di, di-1), ]
  }
}
u <- barplot(as.matrix(data0), las=2, names.arg=rep('', dim(data0)[2]),cex.axis=1.3,
        col = as.character(colors_taxon$col[match(rownames(data0), colors_taxon$taxon)]))

transition <- c('')
arctic <- sts[sts >'155SRF']
atlantic <- sts[sts <='155SRF']
axis(1, at = u[sts %in% transition], col.axis =  'red', labels = sts[sts %in% transition], tick = F,cex.axis=1.3, las=2)
axis(1, at = u[sts %in% arctic], col.axis =  'darkblue', labels = sts[sts %in% arctic], tick = F, cex.axis=1.3,las=2)
axis(1, at = u[sts %in% atlantic], col.axis =  'darkorange', labels = sts[sts %in% atlantic], tick = F,cex.axis=1.3, las=2)
print(rownames(data0))
plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
legend('topleft',legend = rownames(data0)[1:50], 
       fill=as.character(colors_taxon$col[match(rownames(data0), colors_taxon$taxon)])[1:50], 
       box.lty=0 , ncol=2, cex=1)
if (type %in% c('uni_T', 'uni_G')){
  
  sums_u <- apply(data0, 2, sum)
  for (i in 1:dim(data0)[2]){
    data0[,i] <- data0[,i]/sums_u[i]
  }
  if (subset=='0'){
    d_ltp <- dim(data0)[1]

    if (type=='uni_T'){
      rownames(data0)[dim(data0)[1]] <- 'Other'
      if (sum(rownames(data0)=='other')+sum(rownames(data0)=='Other')==2 ){
        inds <- c(which(rownames(data0)=='other'), which(rownames(data0)=='Other'))
        data0[inds[2],] <- data0[inds[1],] + data0[inds[2],]
        data0 <- data0[-inds[1], ]
      }
      rownames(data0)[which(rownames(data0)=='Other')] ='other'
    }

  }
  list_to_plot_in_norm[[1]] <- data0
  u<-barplot(as.matrix(data0), las=2,cex.axis=1.3,names.arg = rep('', dim(data0)[2]),
          col = as.character(colors_taxon$col[match(rownames(data0), colors_taxon$taxon)]))
  axis(1, at = u[sts %in% transition], col.axis =  'red', labels = sts[sts %in% transition], tick = F,cex.axis=1.3, las=2)
  axis(1, at = u[sts %in% arctic], col.axis =  'darkblue', labels = sts[sts %in% arctic], tick = F,cex.axis=1.3, las=2)
  axis(1, at = u[sts %in% atlantic], col.axis =  'darkorange', labels = sts[sts %in% atlantic], tick = F,cex.axis=1.3, las=2)
  plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
  legend('topleft',legend = rownames(data0)[1:50], 
         fill=as.character(colors_taxon$col[match(rownames(data0), colors_taxon$taxon)])[1:50], 
         box.lty=0 , ncol=2, cex=1)
}


par(mfrow=c(2,2), mar=c(5.1, 5.1, 4.1, 2.1))
leg <- NULL
for (j in 2:4){
  data0 <- list_to_plot[[j]]
  sts <- colnames(data0)
  sts <- sapply(sts, function(x){str_replace(x, 'SUR', 'SRF')})
  sts <- as.character(sts)
  colnames(data0) <- sts
  data0 <- data0[,!(colnames(data0) %in% c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]
  sts<-colnames(data0)
  h <- apply(data0, 2, sum)
  #print(h)
  if (type %in% c('uni_G', 'G') & subset == '0'){
    if (type=='uni_G'){
      metaT_data <- readRDS( paste('../data/count_unis_metaG_',clusts0[j],'_',frac,'.rds', sep=''))
      metaT_data <- metaT_data[!(names(metaT_data) %in%
                                c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]
    } else{
      metaT_data <- readRDS(paste('MetaT_',clusts0[j],'_',tax_id,'_', frac,'.rds', sep=''))
      metaT_data <- metaT_data[,!(colnames(metaT_data) %in%
                                c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]
    }
    if (type=='uni_G'){
      data0 <- rbind(data0, metaT_data)
    } else { 
      sums_ue <- apply(metaT_data, 2, sum)
      sums_ug <- apply(data0, 2, sum)
      data0 <- rbind(data0, sums_ue-sums_ug)
    }
    data0[data0<0]<-0
    rownames(data0)[dim(data0)[1]] <- 'Not expressed'
    list_to_plot_in[[j]] <- data0
    if (subset=='0'){
      d_ltp <- dim(data0)[1]
      if (type %in% c('G' , 'uni_G')){
        d_ltp <- d_ltp-1
        root_row=which(rownames(data0)=='root')
        sel=unique(c((d_ltp-os):d_ltp, root_row))
        dt <- apply(data0[sel,], 2, sum)
        anti_sel=1:(d_ltp+1)
        anti_sel=anti_sel[!(anti_sel %in% sel)]
        data0=data0[anti_sel,]
        data0 <- rbind(data0, dt)

      } else{
        root_row=which(rownames(data0)=='root')
        sel=unique(c((d_ltp-os):d_ltp, root_row))
        dt <- apply(data0[sel,], 2, sum)
        anti_sel=1:d_ltp
        anti_sel=anti_sel[!(anti_sel %in% sel)]
        data0=data0[anti_sel,]
        data0 <- rbind(data0, dt)
      }
      rownames(data0)[dim(data0)[1]] <- 'Other'
      if (sum(rownames(data0)=='other')+sum(rownames(data0)=='Other')==2 ){
        inds <- c(which(rownames(data0)=='other'), which(rownames(data0)=='Other'))
        data0[inds[2],] <- data0[inds[1],] + data0[inds[2],]
        data0 <- data0[-inds[1], ]
      }
      rownames(data0)[which(rownames(data0)=='Other')] ='other'
      if (type %in% c('G' , 'uni_G')){
        di <- dim(data0)[1]
        data0 <- data0[c(1:(di-2), di, di-1), ]
      }
      
    }
  }
  if (subset=='0' & type %in% c('uni_T', 'T')){
      h <- apply(data0, 2, sum)

      d_ltp <- dim(data0)[1]
      root_row=which(rownames(data0)=='root')
      sel=unique(c((d_ltp-os):d_ltp, root_row))
      dt <- apply(data0[sel,], 2, sum)
      anti_sel=1:d_ltp
      anti_sel=anti_sel[!(anti_sel %in% sel)]
      data0=data0[anti_sel,]
      data0 <- rbind(data0, dt)
      rownames(data0)[dim(data0)[1]] <- 'Other'
      if (sum(rownames(data0)=='other')+sum(rownames(data0)=='Other')==2 ){
        inds <- c(which(rownames(data0)=='other'), which(rownames(data0)=='Other'))
        data0[inds[2],] <- data0[inds[1],] + data0[inds[2],]
        data0 <- data0[-inds[1], ]
      }
      rownames(data0)[which(rownames(data0)=='Other')] ='other'
  }
  h <- apply(data0, 2, sum)
  print(h)
  data0 <- data0[,!(colnames(data0) %in% c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]
  n_tax <- dim(data0)[1]
  n_to_plot <- min(50,n_tax)
  mxu <- max(h, na.rm=T)
  print(mxu)
  limy <- max(c(1, mxu))
  u<- barplot(as.matrix(data0), las=2,cex.axis=1.3,names.arg=rep('', dim(data0)[2]),ylim=c(0,limy),
          col = as.character(colors_taxon$col[match(rownames(data0), colors_taxon$taxon)]))
  axis(1, at = u[sts %in% transition], col.axis =  'red', labels = sts[sts %in% transition], tick = F,cex.axis=1.3, las=2)
  axis(1, at = u[sts %in% arctic], col.axis =  'darkblue', labels = sts[sts %in% arctic], tick = F, cex.axis=1.3,las=2)
  axis(1, at = u[sts %in% atlantic], col.axis =  'darkorange', labels = sts[sts %in% atlantic], tick = F, cex.axis=1.3,las=2)  
  leg <-append(leg, rownames(data0))
}
leg <- unique(leg)
n_to_plot <- min(50,length(leg))
plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
legend('topleft',legend = leg[1:n_to_plot], 
       fill=as.character(colors_taxon$col[match(leg, colors_taxon$taxon)])[1:n_to_plot], 
       box.lty=0 , ncol=4, cex=0.7)
if (type %in% c('uni_T', 'uni_G')){
  par(mfrow=c(2,2), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 2:4){
    data0 <- list_to_plot[[j]]
    sts <- colnames(data0)
    sts <- sapply(sts, function(x){str_replace(x, 'SUR', 'SRF')})
    sts <- as.character(sts)
    colnames(data0) <- sts
    data0 <- data0[,!(colnames(data0) %in% c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]
    sts<-colnames(data0) 
    if (type=='uni_G' & subset == '0'){
      metaT_data <- readRDS(paste('../data/count_unis_metaG_',clusts0[j],'_',frac,'.rds', sep=''))
      metaT_data <- metaT_data[!(names(metaT_data) %in% 
                                  c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]
      data0 <- rbind(data0, metaT_data)
      data0[data0<0]<-0
      rownames(data0)[dim(data0)[1]] <- 'Not expressed'
      if (subset=='0'){
        d_ltp <- dim(data0)[1]
        if (type %in% c('G' , 'uni_G')){
          d_ltp <- d_ltp-1
          root_row=which(rownames(data0)=='root')
          sel=unique(c((d_ltp-os):d_ltp, root_row))
          dt <- apply(data0[sel,], 2, sum)
          anti_sel=1:(d_ltp+1)
          anti_sel=anti_sel[!(anti_sel %in% sel)]
          data0=data0[anti_sel,]
          data0 <- rbind(data0, dt)
        } else{
          root_row=which(rownames(data0)=='root')
          sel=unique(c((d_ltp-os):d_ltp, root_row))
          dt <- apply(data0[sel,], 2, sum)
          anti_sel=1:(d_ltp)
          anti_sel=anti_sel[!(anti_sel %in% sel)]
          data0=data0[anti_sel,]
          data0 <- rbind(data0, dt)
        }
        rownames(data0)[dim(data0)[1]] <- 'Other'
        if (sum(rownames(data0)=='other')+sum(rownames(data0)=='Other')==2 ){
          inds <- c(which(rownames(data0)=='other'), which(rownames(data0)=='Other'))
          data0[inds[2],] <- data0[inds[1],] + data0[inds[2],]
          data0 <- data0[-inds[1], ]
        }
        rownames(data0)[which(rownames(data0)=='Other')] ='other'
        if (type %in% c('G' , 'uni_G')){
          di <- dim(data0)[1]
          data0 <- data0[c(1:(di-2), di, di-1), ]
        }
      }

    }
    for (i in 1:dim(data0)[2]){
      data0[,i] <- data0[,i]/sums_u[i]
    }
    if (subset=='0' & type=='uni_T'){
      d_ltp <- dim(data0)[1]
      root_row=which(rownames(data0)=='root')
      sel=unique(c((d_ltp-os):d_ltp, root_row))
      dt <- apply(data0[sel,], 2, sum)
      anti_sel=1:(d_ltp)
      anti_sel=anti_sel[!(anti_sel %in% sel)]
      data0=data0[anti_sel,]
      data0 <- rbind(data0, dt)
      rownames(data0)[dim(data0)[1]] <- 'Other'
      if (sum(rownames(data0)=='other')+sum(rownames(data0)=='Other')==2 ){
        inds <- c(which(rownames(data0)=='other'), which(rownames(data0)=='Other'))
        data0[inds[2],] <- data0[inds[1],] + data0[inds[2],]
        data0 <- data0[-inds[1], ]
      }
      rownames(data0)[which(rownames(data0)=='Other')] ='other'
    }
    data0 <- data0[,!(colnames(data0) %in% c('142SRF','201SRF','205SRF', '206SRF', '208SRF', '209SRF', '210SRF'))]
    list_to_plot_in_norm[[j]] <- data0
    n_tax <- dim(data0)[1]
    n_to_plot <- min(50,n_tax)
    u<- barplot(as.matrix(data0), las=2,names.arg=rep('', dim(data0)[2]),cex.axis=1.3,ylim=c(0,1),
            col = as.character(colors_taxon$col[match(rownames(data0), colors_taxon$taxon)]))
    axis(1, at = u[sts %in% transition], col.axis =  'red', labels = sts[sts %in% transition], tick = F, las=2, cex.axis=1.3)
    axis(1, at = u[sts %in% arctic], col.axis =  'darkblue', labels = sts[sts %in% arctic], tick = F, las=2, cex.axis=1.3)
    axis(1, at = u[sts %in% atlantic], col.axis =  'darkorange', labels = sts[sts %in% atlantic], tick = F, las=2, cex.axis=1.3)
    leg <-append(leg, rownames(data0))
  }
}

if (type=='uni_G' | type == 'uni_T' & subset == '0'){
  ty <- paste(type0, 'norm', sep='')  
  if (classi %in% c(5,6)){
    saveRDS(list_to_plot_in_norm, paste('Meta', ty,'_',tax_id,'_', frac,'_',classi,'.rds', sep=''))
  } else{
    saveRDS(list_to_plot_in_norm, paste('Meta', ty,'_',tax_id,'_', frac,'.rds', sep=''))
  }
}
if (type %in% c('uni_G', 'G')){
  ty <- paste(type0, '_nexpr', sep='')
  if (classi %in% c(5,6)){
    saveRDS(list_to_plot_in, paste('Meta',ty ,'_',tax_id,'_', frac,'_',classi,bis,'.rds', sep=''))
  } else{
    saveRDS(list_to_plot_in, paste('Meta',ty ,'_',tax_id,'_', frac,bis,'.rds', sep=''))
  }
}
dev.off()

