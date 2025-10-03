#!/bin/env/usr/env Rscript
library('stringr')

barplot_scaled<-function(t0, exponent, col_vec, title, legend=F, mfrow=F){
  maxy <- apply(t0^(exponent),2, sum)
  maxsumy <- round(max(maxy),1-floor(log(max(maxy))/log(10)))
  maxy0 <- apply(t0,2, sum)
  maxsumy0 <- round(max(maxy0),1-floor(log(max(maxy0))/log(10)))
  if (sum(maxy, na.rm=T) !=0){
  if (is.na(maxsumy) | is.infinite(maxsumy) | is.nan(maxsumy)){
    maxsumy <- max(maxy, na.rm=T)
    maxsumy0 <- max(maxy, na.rm=T)
    print(maxy)
  }
  scaley <- seq(0,maxsumy, maxsumy/6)
  ylabs <- seq(0,maxsumy0, maxsumy0/6)
  ylabs0 <- round(ylabs,digits = 1-floor(log(ylabs)/log(10)))
  if (mfrow==T){
    par(mfrow=c(1,2))
  }
  par(mar=c(7.5, 5.1, 4.1, 4.1))
  sts <- colnames(t0)
  u<- barplot(t0^(exponent),las=2, main=title , ylim=c(0,maxsumy), names.arg=rep('', dim(t0)[2]),
          col = col_vec, yaxt='n', cex.axis=2, cex.lab=2, cex.names=2)
  axis(side = 2, at = scaley, labels = ylabs0, cex.axis=2)
  axis(1, at = u[sts %in% transition], col.axis =  'red', labels = sts[sts %in% transition], tick = F, las=2, cex.axis=1.3)
  axis(1, at = u[sts %in% arctic_stations], col.axis =  'darkblue', labels = sts[sts %in% arctic_stations], tick = F, las=2, cex.axis=1.3)
  axis(1, at = u[sts %in% atlantic_stations], col.axis =  'darkorange', labels = sts[sts %in% atlantic_stations], tick = F, las=2, cex.axis=1.3)
  if (legend==T){
    leg <- rownames(t0)[1:min(50, length(rownames(t0)))]
    plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
    legend('topleft',legend = leg, 
           fill=col_vec[1:min(50, length(rownames(t0)))], box.lty=0 , ncol = 3, cex=0.7)
  }
  }
}
arctic_stations <- c("158SRF","163SRF","168SRF" ,"173SRF", "175SRF", "178SRF","180SRF", "188SRF", "189SRF","191SRF",
                           "193SRF", "194SRF", "196SRF")
atlantic_stations <- c("142SRF","143SRF" ,"144SRF", "145SRF","146SRF" , "147SRF", "148SRF","149SRF" ,"150SRF" ,
                          "151SRF", "152SRF", '155SRF')
transition <- c('')

barplot_GO <- function(GO_list, fraction, taxo, n_clusts, types, name_list, group=NULL, bis){

  arctic_stations <- c("158SRF","163SRF","168SRF" ,"173SRF", "175SRF", "178SRF","180SRF", "188SRF", "189SRF","191SRF",
                       "193SRF", "194SRF", "196SRF")
  atlantic_stations <- c("142SRF","143SRF" ,"144SRF", "145SRF","146SRF" , "147SRF", "148SRF","149SRF" ,"150SRF" ,
                         "151SRF", "152SRF", '155SRF')
  transition <- c('')
  
  n_ty <- length(types)
  GO_data <- rep(list(NULL), n_ty)
  for (i in 1:n_ty){
    GO_data[[i]] <- readRDS(paste('GO_station_table_',fraction,'_',n_clusts,'_',
                                  types[i],'_',taxo,bis,'.rds', sep=''))

  }
  
  GO_data1 <- GO_data[[1]][GO_list]
  if (is.null(group)){
    pdf(paste('GO_meta',types,'_',name_list,'_',fraction,'_',taxo,'_',n_clusts,bis,'.pdf', sep=''),
        width=15, height=15)
  } else{
    pdf(paste('GO_meta',types,'_',name_list,'_',fraction,'_',taxo,'_',
              group,'_',n_clusts,bis,'.pdf', sep=''),
        width=15, height=15)
  }
  for (go in names(GO_data1)){
    title <- paste(GO_table[GO_table[,1]==go,2], go)
    if (!is.null(GO_data1[[go]][[1]])){
      t0 <- as.matrix(GO_data1[[go]][[1]])
      
      d_t0 <- dim(t0)[1]
      
      to_plot <- sum(rownames(t0)!='unknown')
      if (d_t0==1){
        tax <- rownames(t0)
      }
      na <- colnames(t0)
      na <- sapply(na, function(x){str_replace(x, 'SUR', 'SRF')})
      na <- as.character(na)
      colnames(t0) <- na
      t0 <- t0[,!(colnames(t0) %in% c('142SRF','201SRF','205SRF',
                                      '206SRF', '208SRF', '209SRF', '210SRF'))]
      if (d_t0==1){
        na <- na[!(na %in% c('142SRF','201SRF','205SRF',
                             '206SRF', '208SRF', '209SRF', '210SRF'))]
      } else{
        na <- colnames(t0)
      }
      
      if (!is.null(group)){
        t0 <- t0[rownames(t0)==group,]
        t0 <- as.data.frame(t0)
        t0 <-t(t0)
        rownames(t0) <-group
      }
      
      if (to_plot==0 | d_t0==1){
        t0 <- t(t0)
        colnames(t0)<-na
        rownames(t0)<- tax
      }
      other_taxo <- c( 'root','Viruses', 'other Viridiplantae','Insecta', 'Cryptophyta', 'Bacteria', 'Craniata', 'Rhizaria', 'Fungi', 
                       'Streptophyta','Cnidaria', 'Euglenozoa','other', 'Amoebozoa', 'other Alveolata','Archaea', 'Coccosphaerales')

      indices <- which(rownames(t0) %in% other_taxo)
      if (length(indices) > 1){
        indices_good <- which(!(rownames(t0) %in% other_taxo)) 
        to_take <- t0[indices,]
        sum_to_take <- apply(to_take, 2, sum)
        t0_bis <- t0[indices_good, ]
        t0_bis <- rbind(t0_bis, sum_to_take)
        rownames(t0_bis)[dim(t0_bis)[1]] <- 'other'
      } else if(length(indices) ==1){
        indices_good <- which(!(rownames(t0) %in% other_taxo))
        to_take <- t0[indices,]
        t0_bis <- t0[indices_good, ]
        t0_bis <- rbind(t0_bis, to_take)
        rownames(t0_bis)[dim(t0_bis)[1]] <- 'other'
      } else if(length(indices) ==0) {
        t0_bis <- t0
      }
      par(mfrow=c(2,2))
      barplot_scaled(t0 = t0_bis, exponent = 1,title = title,legend = T, mfrow = F,
                     col_vec = as.character(col_taxoS$col[match(rownames(t0_bis),
                                                                col_taxoS$taxon)]))
      if (d_t0==1){
        t1 <- t(t0)
        sum_metaT <- t1
      } else{
        sum_metaT <- apply(t0, 2, sum)
      }
      if (sum(sum_metaT==0)==0){
        cond_atl <- colnames(t0) %in% atlantic_stations
        cond_arc <- colnames(t0) %in% arctic_stations
        vioplot(log(sum_metaT[cond_atl]), log(sum_metaT[cond_arc]), names=c('Atlantic', 'Arctic'),
                col=c('darkorange', 'darkblue'), cex=2.5)
        col1= scales::alpha('firebrick4', 0.8)
        col2= scales::alpha('cyan', 0.8)
        
        stripchart(list(log(sum_metaT[cond_atl]), log(sum_metaT[cond_arc])),cex=4, vertical = TRUE, method = "jitter",
                   pch = 19, add = TRUE, col = c(col1, col2))
        plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
      } else{
        plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
        plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
      }
      
      par(mfrow=c(1,2))
      
      
      
      if (!is.null(group) && group=='unknown'){
        t1 <-t0
        leg <- rownames(t0)
      } else{
        t1 <- t0[rownames(t0)!='unknown',]
        leg <- rownames(t0)[rownames(t0)!='unknown']
      }
      if (to_plot>1){
        if (!is.null(group) && group !='unknown'){
          t1 <- as.data.frame(t1)
          t1 <-t(t1)
          rownames(t1) <-group
        }
        barplot_scaled(t0=t1, exponent = 1, title = title,legend = T, mfrow = T,
                       col_vec = as.character(col_taxoS$col[match(rownames(t1), col_taxoS$taxon)]))
      } else if (to_plot==1){
        par(mar=c(7.5, 5.1, 4.1, 4.1))
        sts <- colnames(t1)
        if (is.null(dim(t1))) {
          lc = length(t1)
        } else{
          lc = dim(t1)[2]
        }
        u<- barplot(t1,las=2, main=title , col = as.character(col_taxoS$col[match(leg, col_taxoS$taxon)]), names.arg=rep('',lc),
                    cex.axis=2, cex.names=2, cex.lab=2)
        axis(1, at = u[sts %in% transition], col.axis =  'red', labels = sts[sts %in% transition], tick = F, las=2, cex.axis=1.3)
        axis(1, at = u[sts %in% arctic_stations], col.axis =  'darkblue', labels = sts[sts %in% arctic_stations], tick = F, las=2, cex.axis=1.3)
        axis(1, at = u[sts %in% atlantic_stations], col.axis =  'darkorange', labels = sts[sts %in% atlantic_stations], tick = F, las=2, cex.axis=1.3)
        plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
        legend('topleft',legend = leg, 
               fill=as.character(col_taxoS$col[match(leg, col_taxoS$taxon)]), box.lty=0 , ncol = 3, cex=0.7)
      }
      exponents <- c(1)
      for (exp in exponents){
        if (n_clusts==4){
          par(mfrow=c(4,2))
        } else if (n_clusts == 2){
          par(mfrow=c(3,2))
        }
        cls <- names(GO_data1[[go]])
        leg <- NULL
        for (cl in cls[2:length(cls)]){
          if (!is.na(cl) & !is.null(GO_data1[[go]][[cl]])){
            t <- GO_data1[[go]][[cl]]
            nas <- colnames(t)
            nas <- sapply(nas, function(x){str_replace(x, 'SUR', 'SRF')})
            nas <- as.character(nas)
            colnames(t) <- nas  
            if (!is.null(group) && group=='unknown'){
              t <- t
              to_plot=1
            } else{
              t <- t[rownames(t)!='unknown', ]
              to_plot <- sum(rownames(t)!='unknown')
            }
            if (!is.null(group)){
              t <- t[rownames(t)==group,]
              t <- as.data.frame(t)
            }
            #to_plot <- sum(rownames(t)!='unknown')
            dimt <- dim(t)[1]
            t <- t[1:(min(20, dimt)),]
            if (to_plot>1){
              t <- t[,!(colnames(t) %in% c('142SRF','201SRF','205SRF',
                                           '206SRF', '208SRF', '209SRF', '210SRF'))]
              barplot_scaled(t0=as.matrix(t), exponent = exp,title = cl, legend = F, mfrow = F,
                             col_vec = as.character(col_taxoS$col[match(rownames(t), col_taxoS$taxon)]))
              leg <- append(leg, rownames(t))
              plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F, main=title)
              legend('topleft',legend = rownames(t),
                     fill=as.character(col_taxoS$col[match(rownames(t), col_taxoS$taxon)]), box.lty=0 , ncol = 3, cex=0.7) 	    
            } else if (to_plot==1){
              t <- t[,!(colnames(t) %in% c('142SRF','201SRF','205SRF',
                                           '206SRF', '208SRF', '209SRF', '210SRF'))]
              leg <-rownames(t)
              par(mar=c(7.5, 5.1, 4.1, 4.1))
              sts <- colnames(t)
              u<-barplot(as.matrix(t),las=2, main=title , cex.axis=2,cex.lab=2,cex.names=2,names.arg=rep('',length(sts)),
                         col = as.character(col_taxoS$col[match(leg, col_taxoS$taxon)]))
              axis(1, at = u[sts %in% transition], col.axis =  'red', labels = sts[sts %in% transition], tick = F, las=2, cex.axis=1.3)
              axis(1, at = u[sts %in% arctic_stations], col.axis =  'darkblue', labels = sts[sts %in% arctic_stations], tick = F, las=2, cex.axis=1.3)
              axis(1, at = u[sts %in% atlantic_stations], col.axis =  'darkorange', labels = sts[sts %in% atlantic_stations], tick = F, las=2, cex.axis=1.3)
              plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
              legend('topleft',legend = leg, 
                     fill=as.character(col_taxoS$col[match(leg, col_taxoS$taxon)]), box.lty=0 , ncol = 3, cex=0.7)
            } else{
              plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
            }
          }
        }
        leg <- unique(leg)
      }
    }
  }
  dev.off()
}

GO_table <- readRDS('GO_table.rds')


tys <- c('T', 'G')
taxo = commandArgs(trailingOnly = T)[1]
n_clusts <- commandArgs(trailingOnly = T)[2]
bis <- commandArgs(trailingOnly = T)[3]
fractions <-c('GGZZ')
if (taxo=='taxo_MGT-v2'){
  col_taxoS <- readRDS('color_table_MGT-v2.rds')
} else if (taxo=='taxo_groups3'){
  col_taxoS <- readRDS('color_table_groups3.rds')
}


for (ty in tys){
for (loc in c('arctic', 'atlantic')){
  for (frac in fractions){
    if (taxo =='taxo_groups3'){
      groups <-c('Hexanauplia', 'Bacillariophyta','Pelagophyceae','Crysophyceae' ,'Coccosphaerales','Phaeocystales' ,
           'Ciliophora', 'Tunicata','Mamiellales', 'Dinophyceae', 'Bacteria', 'Streptophyta', 'Cnidaria','Rhizaria',
           'Insecta', 'Cryptophyta' , 'Euglenozoa', 'Amoebozoa','Craniata',
           'Fungi', 'Viruses','Archaea' ,'Eukaryota (unclassified)','root' ,paste('other ' , c('Haptophyta', 'Opisthokonta', 'Stramenopiles', 'Alveolata', 'Viridiplantae') , sep='' ))
      for (gr in groups){
        GOs <- readRDS(paste('GO_representative_',ty,'_',loc,'_', frac,'_' , 
                             taxo,'_',gr,bis,'.rds', sep=''))
        if (length(GOs)>0){
          barplot_GO(GO_list = GOs, 
                    fraction = frac, taxo = taxo, n_clusts = n_clusts, 
                    types = ty, name_list = loc, group = gr, bis=bis)
        }
      }
      GOs <- readRDS(paste('GO_representative_',ty,'_',loc,'_', frac,'_' , taxo,bis,'.rds', sep=''))
      barplot_GO(GO_list = GOs, 
                 fraction = frac, taxo = taxo, n_clusts = n_clusts, 
                 types = ty, name_list = loc, bis=bis)
    } else{
      GOs <- readRDS(paste('GO_representative_',ty,'_',loc,'_', frac,'_' , taxo,bis,'.rds', sep=''))
      barplot_GO(GO_list = GOs, 
                 fraction = frac, taxo = taxo, n_clusts = n_clusts, 
                 types = ty, name_list = loc, bis=bis)
    }
  }
}
}

