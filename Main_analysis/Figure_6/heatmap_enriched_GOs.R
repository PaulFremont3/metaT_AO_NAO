#library("rgdal")                                                                                                      
#library("raster")
library("ggplot2")
#library('DT')
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
setwd("/env/cns/scratch_TaraOcean/BioAdvection_II/MetaT_4/Groups_metaT/Expression_plots/")


GO_table <- readRDS('GO_table.rds')
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
  signif_data <- readRDS(paste('Significant',enriched,'_uids_GO_annotation_',taxo,'_1.rds', sep=''))
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
fracs <- c('GGZZ','SSUU', 'QQSS',  'KKQQ')
types <- c('physical_clr', 'basins_clr')
vari_l <- c('T_-', 'arctic_clr_+')
vari_l0 <-c('T_+', 'atlantic_clr_+')
varis <- c('T', 'median_basins')

GO_classes <- readRDS('GO_classification_hc.rds')
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
      nopt <- readRDS('n_fit.rds')
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

# pca.res <- PCA(t(res0[[1]]), graph = F, ncp = 6)
# pdf('3D_pca.pdf')
# par(mfrow=c(2,2))
# for (i in 3:6){
#   zz<- scatterplot3d(x = pca.res$ind$coord[,1], y =pca.res$ind$coord[,2] , 
#                      z= pca.res$ind$coord[,i],
#                      xlab = "1st dim", ylab = "2nd dim", 
#                      zlab = paste(i,"th dim", sep='' ))
#   zz.coords <- zz$xyz.convert(pca.res$ind$coord[,1], pca.res$ind$coord[,2] , 
#                               pca.res$ind$coord[,i]) 
#   
#   text(zz.coords$x, 
#        zz.coords$y,             
#        labels=colnames(res0[[1]]) ,              
#        cex = .5, 
#        pos = 4) 
# }
# dev.off()


# stringdistance <- function(x){
#   stringdist(x, unique_goh,method = 'lv' )
# }
# GO_h_dist <- sapply(unique_goh, FUN = stringdistance)
# n <- length(unique(sig_dat$GO_hierarchy))
# GO_h_dist_map <- matrix(GO_h_dist, nrow = n, ncol = n, byrow = T)
# 
# n <- length(unique_goh)
# GO_h_dist_map <-log(GO_h_dist)
# GO_h_dist_map[GO_h_dist_map==-Inf] <- 0
# 
# # herarchical clustering
# res.hc <- pvclust(data=GO_h_dist_map, method.hclust = "average" ,
#                   method.dist = 'cor', nboot = 1, iseed = 1000)
# clusts <- pvpick(res.hc, alpha=0.95, max.only = T)
# clusts <- clusts$clusters
# grps_hc <- rep(1, length(unique_goh))
# count <- 2
# for (c in clusts){
#   vc <- c
#   for (goh in vc){
#     ind <- which(unique_goh==goh)
#     grps_hc[ind] <- count
#   }
#   count <- count+1
# }
# dend <- res.hc %>% as.dendrogram
# lab_dend <-dendextend::get_leaves_attr(dend, attribute ='label' )
# 
# set.seed(2)
# u <- PCA(GO_h_dist_map, graph = F)
# res.km <- kmeans(u$ind$coord, centers = 20, nstart = 25)
# grp <- as.factor(res.km$cluster)
# means_sil_width <- NULL
# sil_res <- rep(list(NULL), 60)
# for (k in 2:60){
#   #v <- cutree(UPGMA, k=k)
#   res.km <- kmeans(u$ind$coord, centers = k, nstart = 25)
#   grp <- as.factor(res.km$cluster)
#   sil <- silhouette(as.numeric(grp), dmatrix = as.matrix(GO_h_dist_map))
#   sil_res[[k]] <- sil
#   means_sil_width <- append(means_sil_width, mean(sil[,3]))
# }
# plot(2:60, means_sil_width, col='red', pch=19, 
#      xlab='Number of clusters', ylab='Mean silhouette width')
# u <- PCA(GO_h_dist_map, graph = F)
# res.km <- kmeans(u$ind$coord, centers = 16, nstart = 25)
# grp <- as.factor(res.km$cluster)
# 
# x <- u$ind$coord[,1]
# y <- u$ind$coord[,2]
# z <- u$ind$coord[,3]
# u$ind$coord[,1] <- (u$ind$coord[,1]-min(u$ind$coord[,1]))*2/(max(u$ind$coord[,1])-min(u$ind$coord[,1]))-1
# u$ind$coord[,2] <- (u$ind$coord[,2]-min(u$ind$coord[,2]))*2/(max(u$ind$coord[,2])-min(u$ind$coord[,2]))-1
# 
# # plot(x,y)
# e1 <- u$eig[1]
# e2 <-u$eig[2]
# e3 <- u$eig[3]
# print(e2/e1)
# e1_var <- e1/sum(u$eig)
# e2_var <- e2/sum(u$eig)
# e3_var <- e3/sum(u$eig)
# print(sum(c(e1_var, e2_var, e3_var)))
# 
# dimension =3
# x_n <- (x/max(abs(x))+1)*255/2
# print(max(x_n))
# if (dimension==3 | dimension ==2){
#   y_n <- ((e2/e1)*y/max(abs(y))+1)*255/2
# } else{
#   y_n <- rep(0, length(x))
# }
# if (dimension==3){
#   z_n <- ((e3/e1)*z/max(abs(z))+1)*255/2
# } else{
#   z_n <- rep(0, length(x))
# }
# colors <-rgb(x_n, y_n, z_n, maxColorValue = 255)
# 
# cols <- c('red', 'violet', 'lightblue', 'darkblue', 'darkgreen', 'darkred', 
#           'orange', 'brown', 'black', 'yellow', 'darkseagreen1', 'goldenrod',
#           'hotpink', 'lightcoral', 'lightgoldenrod')
# set.seed(3)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col_vector <- rep(as.vector(col_vector), 160)
# col_vector[1] <- 'black'
# while('#FFFFFF' %in% col_vector){
#   inds <- which(col_vector=='#FFFFFF')
#   for(j in inds){
#     col_vector[j]=sample(col_vector,1)
#   }
# }
# pdf('GO_hierarchy_plots.pdf')
# plot.PCA(u, axes=c(1, 2), choix="ind", habillage="none", col.ind = cols[grp])
# plot.PCA(u, axes=c(1, 2), choix="ind", habillage="none", col.ind = colors)
# heatmap.2(GO_h_dist_map, trace="none",
#           Rowv = NA, Colv = NA,
#           dendrogram = "none", keysize=1, cexRow = 0.2)
# heatmap.2(GO_h_dist_map, trace="none",Colv = dend, Rowv = dend,
#           labCol=go_unique, labRow=go_unique,
#           keysize=1,cexRow = 0.2, cexCol = 0.2,
#           colRow = col_vector[grps_hc],
#           colCol = col_vector[grps_hc])
# dev.off()
# 
# sig_dat_map <- reshape2::acast(sig_dat, sig_dat$GO~sig_dat$taxo)
# sig_dat_map0 <- sig_dat_map
# sig_dat_map0[sig_dat_map0>0] <- 1
# data.dist <- vegdist(t(sig_dat_map0), method = "euclidean")
# col.clus <- hclust(data.dist, "aver")
# data.dist_r <- vegdist(sig_dat_map0, method = "euclidean")
# row.clus <- hclust(data.dist_r, "aver")
# pdf('test.pdf')
# #colos <- cols[grp[match(sig_dat$GO_hierarchy[match(rownames(sig_dat_map),sig_dat$GO)], 
# #                        unique_goh)]]
# grps <- grp[match(rownames(sig_dat_map0),go_unique)]
# colos <- col_vector[grps]
# # heatmap.2(sig_dat_map0, trace="none",
# #           Rowv = NA, Colv = as.dendrogram(col.clus),
# #           dendrogram = "none", keysize=1, margins=c(10,15), cexRow = 0.2,
# #           colRow = colos)
# heatmap.2(sig_dat_map0, trace="none",
#           Rowv = match(go_unique[order(grp_)],rownames(sig_dat_map0)), Colv = as.dendrogram(col.clus),
#           dendrogram = "none", keysize=1, margins=c(10,15), cexRow = 0.2,
#           colRow = cols[grp[order(grp)]])
# 
# heatmap.2(sig_dat_map0, trace="none",
#           Rowv = match(go_unique[order(grp)],rownames(sig_dat_map0)), Colv = as.dendrogram(col.clus),
#           dendrogram = "none", keysize=1, margins=c(10,15), cexRow = 0.2,
#           colRow = cols[grp[order(grp)]])
# 
# heatmap.2(sig_dat_map0, trace="none",
#           Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus),
#           dendrogram = "none", keysize=1, margins=c(10,15), cexRow = 0.2, cexCol = 0.5)
# dev.off()
