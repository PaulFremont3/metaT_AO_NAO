library("ggplot2")
library('FactoMineR')
library('pvclust')
library('matlab')
library('dplyr')
library('RColorBrewer')
library('stringr')
library('stringdist')
library('tidyr')
library('dendextend')
library('gplots')
library('cluster')
library('parallel')
library('vegan')
library('scatterplot3d')
library('phateR')
library('ecodist')
library('gridExtra')
library('ggrepel')


GO_table <- readRDS('../data/GO_table.rds')
GO_table <- as.data.frame(GO_table)
GO_table <- data.frame(lapply(GO_table, as.character), stringsAsFactors=FALSE)


max_hier <- function(x){
  v <- strsplit(x, ';')[[1]]
  if (length(v)>3){
    u <- paste(v[1:length(v)-1], collapse = ';')
  } else{
    u <- paste(v, collapse = ';')
  }
  return(u)
}
upper_class_go <- function(x){
  v <- strsplit(x, ';')[[1]]
  u <- v[length(v)]
  g <- as.character(GO_table$V2[GO_table$V1==u])
  return(g)
}


prep_data <- function(enriched, taxo, frac, type, vari_s, vari_s0, vari, taxos_sel, count){
  signif_data <- readRDS(paste('../data/Significant',enriched,'_uids_GO_annotation_',taxo,'_1.rds', sep=''))
  sig_dat <- signif_data[[frac]][[type]]
  sig_dat <- sig_dat[sig_dat$vari_s==vari_s | sig_dat$vari_s==vari_s0,]
  sig_dat <- data.frame(lapply(sig_dat, as.character), stringsAsFactors=FALSE)
  sig_dat$GO_hierarchy <- paste(GO_table$V3[match(sig_dat$GO,GO_table$V2)], GO_table$V1[match(sig_dat$GO,GO_table$V2)], sep=';')
  #sig_dat <- sig_dat[sig_dat$GO!='unknown' & sig_dat$GO!='NO_GO',]
  sig_dat$GO_hierarchy1 <- sapply(sig_dat$GO_hierarchy, FUN = max_hier)
  sig_dat$GO_up <- as.character(sapply(sig_dat$GO_hierarchy1, FUN = upper_class_go))
  sig_dat <- sig_dat[with(sig_dat, order(taxo, GO_hierarchy, decreasing = T)),]
  sig_dat <- sig_dat[sig_dat$taxo %in% taxos_sel[[frac]],]
  sig_dat$taxo_var <- paste(sig_dat$taxo, sig_dat$vari_s, sep=' ')
  
  unique_goh <- unique(sig_dat$GO_hierarchy)
  go_unique <- sig_dat$GO[match(unique_goh,sig_dat$GO_hierarchy)]
  
  
  result <- list()
  if (enriched=='_enriched'){
    to_do <- 1
  } else{
    to_do <- 1:3
  }
  for (i in to_do){
    sig_dat_map <- reshape2::acast(sig_dat, sig_dat$GO~sig_dat$taxo_var)
    sig_dat_map0 <- sig_dat_map
    if (count==0){
      sig_dat_map[sig_dat_map>0] <- 1
      sig_dat_map0[sig_dat_map0>0] <- 1
    } else if (count==1){
      print(sum(sig_dat_map0==0.5))
      sig_dat_map0[sig_dat_map0>5] <- 6
      #sig_dat_map0[sig_dat_map0>0 & sig_dat_map0<=5] <- 1
    }
    if (i==2){
      cond <- apply(sig_dat_map0, 1, sum)
      cond <- cond>= length(taxos_sel[[frac]])-1
      sig_dat_map0 <- sig_dat_map0[cond,]
    } else if (i==3){
      cond <- apply(sig_dat_map0, 1, sum)
      cond <- cond<length(taxos_sel[[frac]])
      sig_dat_map0 <- sig_dat_map0[cond,]
    }
    set.seed(1)
    #sig_dat_map0 <- sig_dat_map0[!grepl('NO_GO', rownames(sig_dat_map0)), ]
    if (count==0){
      meth='jaccard'
    } else{
      meth ='bray'
    }
    
    sig_dat_map=sig_dat_map[!(rownames(sig_dat_map) %in% c('unknown', 'NO_GO')), ]  
    data.dist <- vegdist(t(sig_dat_map), method = meth)
    col.clus <- hclust(data.dist, "aver")
    data.dist_r <- vegdist(sig_dat_map, method = meth)
    row.clus <- hclust(data.dist_r, "aver")
    
    data.dist0 <- vegdist(t(sig_dat_map0), method = meth)
    col.clus0 <- hclust(data.dist0, "aver")
    data.dist_r0 <- vegdist(sig_dat_map0, method = meth)
    row.clus0 <- hclust(data.dist_r0, "aver")
    
    if (frac=='GGZZ' & count==1 & type=='physical_clr' & vari=='T' & enriched=='' & i==1){
       write.table(sig_dat_map0, paste('training_set_phate_GO_',frac, '_', type, '_', vari,'.txt', sep=''), row.names=F, col.names=F)
       write.table(sig_dat_map, paste('training_set_phate_GO_all_',frac, '_', type, '_', vari,'.txt', sep=''), row.names=F, col.names=F)
       writeLines(rownames(sig_dat_map0),sep='\t', paste('training_set_phate_GO_rows_', frac, '_', type, '_', vari,'.txt', sep=''))
       writeLines(colnames(sig_dat_map0),sep='\t', paste('training_set_phate_GO_cols_', frac, '_', type, '_', vari,'.txt', sep=''))
    }   
    
    result[[i*6-5]]<-sig_dat_map
    result[[i*6-4]]<-row.clus
    result[[i*6-3]]<-col.clus
    result[[i*6-2]]<-sig_dat_map0
    result[[i*6-1]]<-row.clus0
    result[[i*6]]<-col.clus0
  }
  return(result)
}
plotpca <- function(dat){
  par(mfrow=c(2,1))
  pca.res <- PCA(dat, graph = F, ncp=6)
  if (type=='physical_clr'){
    colo_pca <- c('red', 'blue')[grepl('-', rownames(pca.res$ind$coord))+1]
  } else if (type=='basins_clr'){
    colo_pca <- c('red', 'blue')[grepl('arctic', rownames(pca.res$ind$coord))+1]
  }
  plot.PCA(pca.res ,col.ind=colo_pca, choix='ind')
  plot.PCA(pca.res ,col.ind=colo_pca, choix='ind', axes = c(3,4))
  #plot.PCA(pca.res , choix='var')
}

pltnmds <- function(dat){
  v<- metaMDS(dat, distance = 'bray')
  plot(v)
  ordiplot(v,type="n")
  #orditorp(v,display="species",col="red",air=0.01)
  orditorp(v,display="sites",cex=1.25,air=0.01)
}

plotpcoa <- function(dat){
  varespec <- dat
  if (type=='physical_clr'){
    colo_pca <- c('red', 'blue')[grepl('-', rownames(dat))+1]
  } else if (type=='basins_clr'){
    colo_pca <- c('red', 'blue')[grepl('arctic', rownames(dat))+1]
  }
  varespec.bray <- vegdist(varespec, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan
  
  pcoaVS <- pco(varespec.bray, negvals = "zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999
  
  # y=pcoaVS$vectors[,2]
  # x=pcoaVS$vectors[,1]
  # 
  # wordcloud::textplot(x, y, rownames(dat), col=colo_pca,
  #                     ylim = 1.2*c(min(y),max(y) ), xlim=1.2*c(min(x),max(x)),)
  # points(x, y, col=colo_pca,  pch=19)
  # 
  vals <- round(pcoaVS$values*100/sum(pcoaVS$values)[1],1)
  dt <- data.frame(pcoaVS$vectors)
  par(mfrow=c(2,1))
  v1 <- ggplot(dt, aes(x = X1, y =X2, colour=colo_pca), size=3) +
    geom_point(size=3.5) +
    geom_label_repel(aes(label = rownames(dat), colour=colo_pca), size = 5)+
    scale_color_manual(values = c('blue', 'red'))+
    xlab(paste("PCoA 1 (",vals[1],"%)", sep='')) +
    ylab(paste("PCoA 2 (",vals[2],"%)", sep='')) + theme_bw(base_size=18)+theme(legend.position = "none", axis.text = element_text(size=20), 
                                                                                axis.title = element_text(size=20))
  
  v2 <- ggplot(dt, aes(x = X3, y =X4, colour=colo_pca)) +
          geom_point(size=3.5) +
          geom_label_repel(aes(label = rownames(dat), colour=colo_pca), size = 5)+
          scale_color_manual(values = c('blue', 'red'))+
          xlab(paste("PCoA 3 (",vals[3],"%)", sep='')) +
          ylab(paste("PCoA 4 (",vals[4],"%)", sep='')) + theme_bw(base_size=18)+theme(legend.position = "none", axis.text = element_text(size=20), 
                                                                                      axis.title = element_text(size=20))
  grid.arrange(v1, v2, nrow=2)
}
#u <- phate(data = dat)

taxos_sel <- rep(list(NULL), 4)
names(taxos_sel) <- c('GGZZ', 'SSUU', 'QQSS', 'KKQQ')
taxos_sel[['GGZZ']] <- c(  'Mamiellales','Pelagophyceae',
                          'Phaeocystales', 
                          'Bacillariophyta') #'Bacteria', 'Dinophyceae','unknown','Cryptophyta',
taxos_sel[['SSUU']] <- c( 'unknown', 'Insecta','Hexanauplia')
taxos_sel[['KKQQ']] <- c( 'unknown', 'Mamiellales','Bacillariophyta', 'Bacteria', 
                          'Dinophyceae','Hexanauplia')
taxos_sel[['QQSS']] <- c( 'unknown', 'Bacillariophyta', 
                          'Bacteria', 'Dinophyceae', 'Hexanauplia')

taxo <- 'groups3'
fracs <- c('GGZZ')
types <- c('physical_clr')
vari_l <- c('T_-')
vari_l0 <-c('T_+')
varis <- c('T')

GO_classes <- readRDS('../data/GO_classification_hc.rds')
GO_classes$GO_id <- GO_table[match(GO_classes$GO, GO_table[,1]) ,2]
for (frac in fracs){
  c=1
  
  for (type in types){
    for (count in c(1,0)){
      vari_s <- vari_l[c]
      vari_s0 <- vari_l0[c]
      vari <- varis[c]
      res0 <- prep_data(enriched='', taxo, frac, type, vari_s, vari_s0, vari, taxos_sel, count)
      res_en <- prep_data(enriched='_enriched', taxo, frac, type, vari_s, 
                          vari_s0, vari, taxos_sel, count)
      nopt <- readRDS('../data/n_fit.rds')
      pdf(paste('GO_heatmap_', taxo, '_', frac, '_', type, '_', vari,'_',count,'.pdf' , sep=''),
          height = 15, width=10)
      for (i in 3:1){
        colo <- match(rownames(res0[[i*6-2]]), GO_classes$GO_id)
        colo[is.na(colo)] <- which(GO_classes[[paste('color2_', nopt, sep='')]]=='black')[1]
        clrw <- GO_classes[[paste('color2_', nopt, sep='')]][colo]
        grps <- GO_classes[[paste('Group2_', nopt, sep='')]][colo]
        if (type == 'basins_clr'){
          colr <- c('red','blue')[grepl('arctic', colnames(res0[[1]]))+1]
        } else{
          colr <- c('red','blue')[grepl('-', colnames(res0[[1]]))+1]
        }
        nr=nrow(res0[[i*6-5]])
        nc=ncol(res0[[i*6-5]])
        if (nr<=110){
          fc <- 0.3
        } else{
          fc <- 0
        }
        if (count==0){
          col_vec = heat.colors(15)
        } else if (count==1){
          #col_vec= c(heat.colors(15)[1], heat.colors(15)[15],'green', 'blue',  'purple', 'darkmagenta', 'black')[1:(max(as.matrix(res0[[i*6-2]]))+1)]
          col_vec= c('white', colorRampPalette(c('goldenrod1', 'darkorchid'))(6))[1:(max(as.matrix(res0[[i*6-2]]))+1)]
        }
        heatmap.2(res0[[i*6-2]], trace="none",dendrogram = 'column',col = col_vec,
                  Rowv = as.dendrogram(res0[[i*6-1]]), Colv = as.dendrogram(res0[[i*6]]),
                  keysize=1, margins=c(14,19), cexRow = fc + 1/(1*log10(nr)), cexCol = 0 + 1/log10(nc),
                  colRow = clrw, colCol = colr)
        # if (i %in% c(2,1)){
        
        plotpcoa(t(res0[[i*6-2]]))
        
        col_vec <- c('red', colorRampPalette(c('lightblue', 'blue', 'darkblue'))(max(res0[[i*6-5]])))
        heatmap.2(res0[[i*6-5]], trace="none",dendrogram = 'column',col = col_vec,
                  Rowv = as.dendrogram(res0[[i*6-4]]), Colv = as.dendrogram(res0[[i*6-3]]),
                  keysize=1, margins=c(14,19), cexRow = fc + 1/(1*log10(nr)), cexCol = 0 + 1/log10(nc),
                  colRow = clrw, colCol = colr)
        plotpcoa(t(res0[[i*6-5]]))
        
        # }
      }
      
      ord <- order(grps)
      dat <- res0[[4]][ord,]
      cols <- clrw[ord]
      nr=nrow(dat)
      nc=ncol(dat)
      if (nr<=110){
        fc <- 0.3
      } else{
        fc <- 0
      }
      heatmap.2(dat, trace="none",dendrogram = 'column',col = col_vec,
                Rowv = NULL, Colv = as.dendrogram(res0[[6]]),
                keysize=1, margins=c(14,19), cexRow=fc + 1/(1*log10(nr)), cexCol = 0 + 1/log10(nc),
                colRow = cols, colCol = colr)
      
      colo <- match(rownames(res_en[[4]]), GO_classes$GO_id)
      colo[is.na(colo)] <- which(GO_classes[[paste('color2_', nopt, sep='')]]=='black')[1]
      res_en[[7]] <- GO_classes[[paste('color2_', nopt, sep='')]][colo]
      res_en[[8]] <- GO_classes[[paste('Group2_', nopt, sep='')]][colo]
      if (type == 'basins_clr'){
        colr <- c('red','blue')[grepl('arctic', colnames(res_en[[1]]))+1]
      } else{
        colr <- c('red','blue')[grepl('-', colnames(res_en[[1]]))+1]
      }
      nr=nrow(res_en[[1]])
      nc=ncol(res_en[[1]])
      if (nr<=110){
        fc <- 0.3
      } else{
        fc <- 0
      }
      heatmap.2(res_en[[4]], trace="none",dendrogram = 'column',col = col_vec,
                Rowv = as.dendrogram(res_en[[5]]), Colv = as.dendrogram(res_en[[6]]),
                keysize=1, margins=c(14,19), cexRow = fc + 1/(1*log10(nr)), cexCol = 0 + 1/log10(nc),
                colRow = res_en[[7]], colCol = colr)
      ord <- order(res_en[[7]])
      dat <- res_en[[4]][ord,]
      cols <- res_en[[7]][ord]
      nr=nrow(dat)
      nc=ncol(dat)
      if (nr<=110){
        fc <- 0.3
      } else{
        fc <- 0
      }
      heatmap.2(dat, trace="none",dendrogram = 'column',col = col_vec,
                Rowv = NULL, Colv = as.dendrogram(res_en[[6]]),
                keysize=1, margins=c(14,19), cexRow = fc + 1/(1*log10(nr)), cexCol = 0 + 1/log10(nc),
                colRow = cols, colCol = colr)
      dev.off()
    }
    c=c+1
  }
}

