library('pvclust')
library('dplyr')
library('tidyr')
library('dendextend')
library('gplots')
library('cluster')
library('RColorBrewer')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/simkaMin")

all_data0 <- readRDS('simkaTG_matrix_sur.rds')
all_data_arctic0 <- readRDS('simkaTG_matrix_arctic.rds')
bc_unigenes <- NULL
frac0 <- c('GGZZ')
all_data_arctic0_uni <- rep(list(NULL), 2)
names(all_data_arctic0_uni) <- c('T', 'G')
for (f in frac0){
  fi <- readRDS(paste('MetaTG_bray-curtis_unigenes_',f,'.rds', sep=''))
  all_data_arctic0_uni[['T']][[f]] <- fi[[1]]
  all_data_arctic0_uni[['G']][[f]] <- fi[[2]]
}
fracs <- c('GGZZ')
types <-c('T', 'G')
cols <- colorRampPalette(c('darkred', 'red','orange', 'green','blue' ))(50)
br <- c(seq(0.3, 1, length.out = 51))
clusts <- function(data, name, maxc){
  for (j in 1:2){
    for (i in 1:length(data[[1]])){
      if (grepl(name, 'uni')){
        fracs_b <- frac0
      } else{
        fracs_b <- fracs
      }
      mat <- data[[j]][[i]]
      cond <- grepl('SUR.1', colnames(mat)) | grepl('SUR.2', colnames(mat)) | grepl('SUR.3', colnames(mat)) 
      mat <- mat[!cond,!cond]
      v <-as.integer(unlist(strsplit(colnames(mat), 'SUR')))
      mat <- mat[order(v), order(v)]
      
      if (name %in% c('arctic', 'arctic_sub', 'arctic_uni_bc_sub') & fracs_b[i] %in% c('SSUU', 'QQSS', 'GGZZ', 'KKQQ')){
        stations <- c("142SUR" ,"144SUR", "145SUR", "147SUR", "148SUR", "150SUR" ,"151SUR" ,"152SUR","155SUR", "158SUR", "163SUR", "168SUR", "173SUR", "175SUR", "178SUR", "180SUR","188SUR", "189SUR", "191SUR", "193SUR", "194SUR", "196SUR" ,"201SUR" ,"205SUR","206SUR" ,"208SUR", "209SUR", "210SUR")
        if (fracs_b[i]=='GGZZ'){
          stations <- stations[stations != "191SUR"]
        }
        if (name=='arctic_sub' | name=='arctic_uni_bc_sub'){
          good_stat <- stations
          good_stat <- good_stat[!(good_stat %in% c('142SUR','201SUR','205SUR', '206SUR', '208SUR', '209SUR', '210SUR'))]
          stations <- stations[stations %in% good_stat]
        }
        cond <- colnames(mat) %in% stations
        mat <- mat[cond, cond]
      } else{
        if (fracs_b[i]=='GGZZ'){
          cond <- colnames(mat) != "191SUR"
          mat <- mat[cond, cond]
        }
        if (fracs_b[i] %in% c('GGZZ_GGMM', 'KKQQ_MMQQ') & name=='arctic_sub'){
          cond <- !(colnames(mat) %in% c('142SUR','201SUR','205SUR', '206SUR', '208SUR', '209SUR', '210SUR'))
          mat <- mat[cond, cond]
        }
      }
      
      set.seed(2)
      UPGMA <- pvclust(mat, method.dist="cor", method.hclust="average", nboot=1000)
      if (name %in% c('arctic', 'arctic_sub', 'arctic_uni_bc', 'arctic_uni_bc_sub')){
        csize = 1
      } else if (name=='all'){
        csize=0.6
      }
      dend <- UPGMA %>% as.dendrogram
      pdf(paste('UPGMA_meta',types[j],'_',name,'_', fracs_b[i], '_test.pdf', sep=''), width = 15, height = 12)
      plot(UPGMA, cex=csize)
      pvrect(UPGMA, alpha = 0.85, max.only = F)
      
      if (name %in% c('arctic', 'arctic_sub', 'arctic_uni_bc', 'arctic_uni_bc_sub')){
        UPGMA %>% as.dendrogram %>% 
        set("branches_k_color", k = 3, value = c("green", "blue", 'red')) %>%
        plot(main='k=3 clusters')
        UPGMA %>% text
        UPGMA %>% pvrect(alpha = 0.90, max.only = F)
        
        UPGMA %>% as.dendrogram %>% 
        set("branches_k_color", k = 2, value = c("green", "blue")) %>%
        plot(main='k=2 clusters')
        UPGMA %>% text
        UPGMA %>% pvrect(alpha = 0.85, max.only = F)
      }
      #       means_sil_width <- NULL
      #       sil_res <- rep(list(NULL), 10)
      #       for (k in 2:maxc){
      #         v <- cutree(UPGMA, k=k)
      #         u <- silhouette(v, mat)
      #         sil_res[[k]] <- u
      #         means_sil_width <- append(means_sil_width, mean(u[,3]))
      #       }
      #       plot(2:maxc, means_sil_width, col='red', pch=19, xlab='Number of clusters', ylab='Mean silhouette width')
      #       
      #       opt_n<- which.min(means_sil_width)+1
      #       plot(sil_res[[opt_n]])
      
      #       UPGMA %>% as.dendrogram %>% 
      #       set("branches_k_color", k = opt_n, value = brewer.pal(opt_n, "Dark2")) %>%
      #       plot(main=paste('k=',opt_n ,' clusters', sep=''))
      #       UPGMA %>% text
      #       UPGMA %>% pvrect(alpha = 0.75, max.only = F)
      if (grepl('arctic', name)){
        colos <- rep(NA,length(colnames(mat)))
        colos[colnames(mat)>"155SUR"] <- 'blueviolet'
        colos[colnames(mat)<="155SUR"] <- 'darkblue'
        #colos[colnames(mat)>="155SUR" & colnames(mat)<="163SUR"] <- 'red'
      }
      t1 <- pvpick(UPGMA, alpha = 0.75, max.only = F)
      
      mat <- as.matrix(mat)
      diag(mat) <- NA
      if (name %in% c('arctic', 'arctic_sub', 'arctic_uni_bc', 'arctic_uni_bc_sub')){
        fsize = 1.5
      } else if (name=='all'){
        fsize=0.6
      }
      if (!(grepl('arctic', name))){
        heatmap.2(mat, trace="none", keysize=1,margins=c(10,22), labCol=colnames(mat), labRow=colnames(mat), col= cols,
                  symkey = F, breaks = br, cexRow = fsize, cexCol = fsize,dendrogram = 'none',
                  Colv = NA, Rowv = NA, main='no cluster, no duplicates')
      #     heatmap.2(mat,trace="none",hclustfun = function(x){ set.seed(2);pvclust(mat, method.dist="cor", method.hclust="average", nboot=1000)}, keysize=1,margins=c(10,22), labCol=colnames(mat), labRow=colnames(mat), col= cols,
      #               symkey = F, breaks = br, cexRow = fsize, cexCol = fsize, main='upgma, no duplicates')
        heatmap.2(mat,trace="none",Colv = dend, Rowv = dend,  keysize=1,margins=c(10,22), labCol=colnames(mat), labRow=colnames(mat), col= cols,
                  symkey = F, breaks = br, cexRow = fsize, cexCol = fsize, main='upgma, no duplicates')
      
       } else{
         heatmap.2(mat, trace="none", keysize=1,margins=c(10,22), labCol=colnames(mat), labRow=colnames(mat), col= cols,
                  symkey = F, breaks = br, cexRow = fsize, cexCol = fsize,dendrogram = 'none', colRow=colos, colCol=colos,
                  Colv = NA, Rowv = NA, main='no cluster, no duplicates')
         heatmap.2(mat,trace="none",Colv = dend, Rowv = dend,  keysize=1,margins=c(10,22), labCol=colnames(mat), labRow=colnames(mat), col= cols,colRow=colos, colCol=colos,
                  symkey = F, breaks = br, cexRow = fsize, cexCol = fsize, main='upgma, no duplicates')
       }
       dev.off()
    }
  }
}
clusts(all_data_arctic0_uni, name = 'arctic_uni_bc_sub', maxc = 10)
#clusts(all_data_arctic0_uni, name = 'arctic_uni_bc', maxc = 10)
#clusts(all_data_arctic0, name = 'arctic', maxc = 10)
#clusts(all_data_arctic0, name = 'arctic_sub', maxc = 10)
#clusts(all_data0, name = 'all', maxc = 20)

