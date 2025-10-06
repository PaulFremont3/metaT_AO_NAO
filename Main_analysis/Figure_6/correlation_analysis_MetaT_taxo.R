#!/bin/env/usr/env Rscript
library('treemap')
#library("readxl")
library('gplots')
library('FactoMineR')
library('stringr')
library('RColorBrewer')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")

taxo = commandArgs(trailingOnly = T)[1]
subset = commandArgs(trailingOnly = T)[2]
set.seed(4)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- rep(as.vector(col_vector), 160)


while('#FFFFFF' %in% col_vector){
  inds <- which(col_vector=='#FFFFFF')
  for(j in inds){
    col_vector[j]=sample(col_vector,1)
  }
}

qual_col_pals1 = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector1 = unlist(mapply(brewer.pal, qual_col_pals1$maxcolors, rownames(qual_col_pals1)))
col_vector1 <- rep(col_vector1, 80)

while('#FFFFFF' %in% col_vector1){
  inds <- which(col_vector1=='#FFFFFF')
  for(j in inds){
    col_vector[j]=sample(col_vector1,1)
  }
}

if (taxo=='taxo'){
  data_uni_tax <- readRDS('taxID_uni.rds')
} else {
  data_uni_tax <- readRDS(paste('data_uni_',taxo,'.rds', sep=''))
  colnames(data_uni_tax)<-c('taxName', 'geneID')
  data_uni_tax<- data_uni_tax[!duplicated(data_uni_tax$geneID),]
  if (is.factor(data_uni_tax$taxName)){
    data_uni_tax$taxName <- as.character(levels(data_uni_tax$taxName))[data_uni_tax$taxName]
  }
} 

GO_table <- readRDS('GO_table.rds')
pfam2go <- read.table('pfam2go_27_02_20.txt')
allg <- as.character(unique(pfam2go$V3))
go_names <- GO_table[match(allg, GO_table[,1]),2]
all_gos <- c(go_names, 'unknown', 'NO_GO')
go_colors <- data.frame(go = all_gos, color=col_vector1[1:length(all_gos)])

extract_GO <- function(x){
  ud <- as.character(x[1][[1]])
  pfams <- as.character(x[5][[1]])
  tax <- x[6][[1]]
  vari <- as.character(x[2][[1]])
  cor_sign <- as.character(x[9][[1]])
  if (!is.na(pfams)){
    pfams <- str_split(pfams, ';')[[1]]
    pfams <- pfams[2:length(pfams)]
    GOs <- NULL
    passed <-c(0)
    for (pf in pfams){
      go <- pfam2go$V3[as.character(pfam2go$V1)==pf]
      if (length(go)!=0){
        GOS <- GO_table[GO_table[,1] %in% go,c(1,2,7)]
        if (!is.null(dim(GOS))){
          index <- which.max(GOS[,3])
          GO <- GOS[index,2]
        } else{
          GO <- GOS[2]
        }
        if (!(GO %in% passed)){
          GOs <- append(GOs, c(tax,vari,cor_sign,GO, ud))
        }
        passed <- append(passed, GO)
      }
    }
    if(is.null(GOs)){
      GOs <- c(tax, vari,cor_sign,'NO_GO', ud)
    }
  } else{
    GOs <- c(tax, vari,cor_sign, 'unknown', ud)
  }
  
  return(GOs)
}

barplot_tax <- function(components, name, frac, subset){
  signif_uids_all <- NULL
  for (comp in components){
    signif_uids <- readRDS(paste('Significant_uids_all_',frac,'_0_05_',comp,'.rds', sep=''))
    signif_uids$taxo <- data_uni_tax$taxName[match(signif_uids$uni, data_uni_tax$geneID)]
    signif_uids$taxo[is.na(signif_uids$taxo)] <- 'unknown'
    if (name %in% c('basins', 'basins_clr')){
      signif_uids$vari <- comp
      signif_uids$cor <- 1
    }
    signif_uids_all <- rbind(signif_uids_all, signif_uids)
  }
  if (name %in% c('basins', 'basins_clr')){
    signif_uids_all <- signif_uids_all[,c(1,5,6,2:4)]
  }
  signif_uids_all$id <- paste(signif_uids_all$uni, signif_uids_all$vari, sep='_')
  signif_uids_all <- signif_uids_all[!duplicated(signif_uids_all$id),]
  
  signif_uids_all <- signif_uids_all[!is.na(signif_uids_all$uni),]
  signif_uids_all$count <-1
  
  signif_uids_all$cor_sign <- as.numeric(signif_uids_all$cor>0)
  signif_uids_all$cor_sign[signif_uids_all$cor_sign==0] <-'-'
  signif_uids_all$cor_sign[signif_uids_all$cor_sign==1] <-'+'
  signif_uids_all$cor_sign<-as.factor(signif_uids_all$cor_sign)
  
  test0 <- apply(signif_uids_all,1, extract_GO)
  test0 <- matrix(unlist(test0), ncol=5, byrow = T)
  test0 <- as.data.frame(test0)
  
  test <- reshape2::acast(signif_uids_all, value.var='count', 
                          signif_uids_all$taxo~signif_uids_all$vari+signif_uids_all$cor_sign)
  sums_tax <- apply(test, 1, sum)
  test <- test[order(sums_tax, decreasing = T),]
  signif_enriched_all <- NULL
  signif_all <- NULL
  if (subset=='0'){ 
    pdf(paste('correlation_',name,'_taxo_',taxo,'_', frac, '.pdf', sep=''), width = 15, height = 7)
  } else if (subset=='1'){
    pdf(paste('correlation_',name,'_taxo_',taxo,'_', frac, '_subset.pdf', sep=''), width = 15, height = 7)
  }
  par(mfrow=c(1,2))
  if (frac == 'GGZZ'){ 
    d_ltp <- dim(test)[1]
    os <- 17
    dt <- apply(test[(d_ltp-os):d_ltp,], 2, sum)
    data0 <- test[1:(d_ltp-(os+1)),]
    #rownames(data0) <- rownames(test)[1:(d_ltp-(os+1))]
    data0 <- rbind(data0, dt)
    rownames(data0) <- c(rownames(test)[1:(d_ltp-(os+1))] ,'other')
  } else{
    data0 <- test
  }
  #print( rownames(data0) )
  #print( col_list$taxon )
  #print( match(rownames(data0), col_list$taxon) )
  #print( col_list$col[match(rownames(data0), col_list$taxon)] )
  #print( as.character(col_list$col[match(rownames(data0), col_list$taxon)]))
  if ('Mamiellales' %in% rownames(data0) & 'Bacillariophyta' %in% rownames(data0) & 'Pelagophyceae' %in% rownames(data0) & 'Phaeocystales' %in% rownames(data0)){
    data1=data0[rownames(data0) %in% c('Mamiellales', 'Bacillariophyta', 'Pelagophyceae', 'Phaeocystales'), ]
    barplot(data1, col=col_list$col[match(rownames(data1), col_list$taxon)], las=1 )
    plot(15,15, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
  
    #print(rownames(data1))
    legend('topleft',legend = rownames(data1),
           fill=as.character(col_list$col[match(rownames(data1), col_list$taxon)]),
           box.lty=0 ,ncol=3, cex=0.7)
  }

  barplot(data0, col=col_list$col[match(rownames(data0), col_list$taxon)], las=1 )
  plot(15,15, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
  legend('topleft',legend = rownames(data0), 
         fill=as.character(col_list$col[match(rownames(data0), col_list$taxon)]), 
         box.lty=0 ,ncol=3, cex=0.7)
  
  
  
  tot <- min(11,dim(test)[1])
  par(mfrow=c(1,2))
  barplot(test[2:tot,], col=as.character(col_list$col[match(rownames(test)[2:tot], col_list$taxon)]), las=1 )
  plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
  legend('topleft',legend = rownames(test)[2:tot], 
         fill=as.character(col_list$col[match(rownames(test)[2:tot], col_list$taxon)]), box.lty=0, cex=0.7 )

  if (subset=='1'){
    baseline_matou <- readRDS(paste('baseline_GO_matou_',taxo,'s_',frac,'_subset.rds', sep=''))
  } else{
    baseline_matou <- readRDS(paste('baseline_GO_matou_',taxo,'s_',frac,'.rds', sep=''))
  }
  if (taxo %in% c('MGT-v2','MGT-v2-phylum')){
    M=15
  } else if (taxo=='MGT'){
    M=5
  } else if (taxo=='groups'){
    M=9
  } else if (taxo=='groups2'){
    M=19
  } else if (taxo=='groups3'){
    M=29
  } else if (taxo=='taxo'){
    M=8
  }

  for (tax in rownames(test)[1:M]){
    test1 <- test0[test0$V1==tax,]
    colnames(test1)<-c('taxo', 'vari', 'cor_sign', 'GO', 'uid')
    test2 <- reshape2::acast(test1,  
                             test1$GO~test1$vari+test1$cor_sign, 
                             fun.aggregate = length )
    title_vec <- colnames(test2)
    if (dim(test2)[1]>1){
      test2 <- test2[order(apply(test2, 1, sum), decreasing = T),]
      test2 <-as.matrix(test2)
      colnames(test2) <-title_vec
      par(mfrow=c(1,2))
      barplot(test2, las=1, 
              col=as.character(go_colors$color[match(rownames(test2), go_colors$go)]), main=tax)
      plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
      legend('topleft',legend = rownames(test2), 
             fill=as.character(go_colors$color[match(rownames(test2), go_colors$go)]), 
             box.lty=0 , ncol=3, cex=0.5)
    
    
      par(mfrow=c(1,2))
      to_plot <- !(rownames(test2) %in% c('unknown', 'NO_GO'))
      if (sum(to_plot)>1){
        test2_1 <- as.data.frame(test2[to_plot,])
        m_ind <- min(30, dim(test2_1)[1])
        test2_2 <- as.matrix(test2_1[1:m_ind,])
        colnames(test2_2) <- title_vec
        barplot(test2_2, las=1,
                col=as.character(go_colors$color[match(rownames(test2_1)[1:m_ind], go_colors$go)]), main=tax)
        plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
        legend('topleft',legend = rownames(test2_1)[1:m_ind], 
               fill=as.character(go_colors$color[match(rownames(test2_1)[1:m_ind], go_colors$go)]),
               box.lty=0 , ncol=2, cex=0.7)
      }
      baseline_matou_tax <- baseline_matou[baseline_matou$Group.1==tax, c(2,3)]
      baseline_matou_tax$Group.2 <- as.character(levels(baseline_matou_tax$Group.2))[baseline_matou_tax$Group.2]
      baseline_matou_tax$Group.2[baseline_matou_tax$Group.2=='Unknown']<-'unknown'
      colnames(baseline_matou_tax) <- c('subgroup', 'value')
      test3 <- test2
      p_vals <- NULL
      for (i in 1:dim(test2)[2]){
        data <- data.frame(value=test2[,i], subgroup=rownames(test2))
        p_values<- sapply(1:dim(data)[1], FUN = hypergeometric_test, data=data, baseline=baseline_matou_tax)
        p_values_adj <- p.adjust(p_values, method = 'hochberg')
        #p_values_adj[p_values_adj==0]<- min(p_values_adj[p_values_adj !=0])
        if (min(p_values_adj[p_values_adj !=0])==Inf){
          p_values_adj[p_values_adj==0]=0.00001
        } else if (min(p_values_adj[p_values_adj !=0])==-Inf){
	        p_values_adj[p_values_adj==0]= 0.00001
         # p_values_adj[p_values_adj==0]= min(p_values_adj[p_values_adj !=0])
        } else{
	        p_values_adj[p_values_adj==0]= 0.00001
	      }
        p_values_adj<- -log10(p_values_adj)
        test3[!(p_values_adj>-log10(0.05)), i] <- 0
        p_vals <- cbind(p_vals, p_values_adj)
      }
      selec <- apply(test3,1,sum)!=0
      if (sum(selec!=0)){
        test3 <- test3[selec,]
        p_vals <- p_vals[selec,]
        p_vals[test3==0]<-0
        p_vals <- as.data.frame(p_vals)
        if (sum(selec)>1){
          colnames(p_vals)<- colnames(test2)
          rownames(p_vals)<- rownames(test2)[selec]
        } else{
          
        }
      } else{
        test3 <-NULL
      }
      
     
      
      if (!is.null(test3)){
        par(mfrow=c(1,2))
        legs <- rownames(test2)[which(selec==T)]
        barplot(test3, las=1, 
                col=as.character(go_colors$color[match(legs, go_colors$go)]), main=paste(tax, 'Enriched GO'))
        plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F)
        legend('topleft',legend = legs[1:min(30, sum(selec))], 
               fill=as.character(go_colors$color[match(legs[1:min(30, sum(selec))], go_colors$go)]), 
               box.lty=0 , ncol=2, cex=0.7)
        par(mfrow=c(1,1))
        for (var in colnames(test3)){
          selection <- test3[,which(colnames(test3)==var)]!=0
          if (sum(selection)>0){
            plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F, main=var)
            legend('topleft',legend = legs[selection],
                   fill=as.character(go_colors$color[match(legs[selection], go_colors$go)]),
                   box.lty=0 , ncol=2, cex=0.7)
          }
        }
        
        test1$vari_s <- paste(test1$vari, test1$cor_sign, sep='_')
	co <- 1
        for (vs in colnames(test3)){
	  GEs <- rownames(test3)[test3[,co] != 0]
          signif_enriched <- test1[test1$vari_s==vs & test1$GO %in% GEs,]
	  signif <- test1[test1$vari_s==vs,]
	  signif_enriched_all <- rbind(signif_enriched_all, signif_enriched)
	  signif_all <- rbind(signif_all, signif)
	  co <- co+1
        }
      }
      extract_parent_GO <- function(g){
        zero_level <- c('molecular_function', 'cellular_component', 'biological_process')
        if (!(g %in% c('unknown', 'NO_GO'))){
          prts <- strsplit(GO_table[GO_table[,2]==g,4], ';')[[1]]
          u <- prts[1]
          v <- prts[1]
          c=1
          while (v %in% zero_level){
            v <- strsplit(GO_table[GO_table[,2]==g,4], ';')[[1]][c+1]
            c=c+1
          }
          if (is.na(v)){
            v <- prts[1]
          }
        }else{
          u <- g
          v <- g
        }
        return(c(u,v))
      }
      if (sum(selec)>1){
        go_parents <- t(sapply(rownames(p_vals), extract_parent_GO))
        for (i in 1:dim(test2)[2]){
          dat <- data.frame(parent=go_parents[,1],group=go_parents[,2], subgroup=rownames(p_vals), value =p_vals[,i])
          if (sum(dat$value)!=0){
            treemap(dtf = dat,vSize = 'value',index = c('parent','group', 'subgroup') ,
                    title=title_vec[i], overlap.labels=0.8, 
                    align.labels = list(c('left', 'top'),c('right', 'top'), c('center', 'center')))
          }
        }
      }
    }
  }
  dev.off()
  return(list(signif_enriched_all, signif_all))
}

hypergeometric_test <- function(i, data, baseline){
  q0=data$value[i]
  if (q0!=0){
    # m0=baseline$value[baseline$subgroup==as.character(data$subgroup[i])]
    m0=sum(baseline$value[baseline$subgroup==as.character(data$subgroup[i])])
    # n0=sum(baseline$value)-baseline$value[baseline$subgroup==as.character(data$subgroup[i])]
    n0=sum(baseline$value)-sum(baseline$value[baseline$subgroup==as.character(data$subgroup[i])])
    k0=sum(data$value)
    proba<-phyper(q=q0, m=m0, n=n0,k=k0, lower.tail = F, log.p = F)
  }else{
    proba <-1
  }
  hypergeometric_test<-proba
}
comp_list <- list(c('T_clr', 'Sal_clr', 'SSD_clr'), c('Phos_clr', 'Si_clr', 'NO3_clr'),c('Si._clr', 'N._clr', 'Fe_clr'),
                  c('arctic_clr', 'atlantic_clr'),
                  c('pca1_pca', 'pca2_pca', 'pca3_pca'), c('T_log', 'Sal_log', 'SSD_log'),
                  c('Phos_log', 'Si_log', 'NO3_log'), c('Si._log', 'N._log', 'Fe_log'),
                  c('Environment', 'Travel_time'),c('arctic', 'atlantic'))
names0 <- c('physical_clr', 'nutrients1_clr', 'nutrients2_clr','basins_clr', 'pca','physical', 'nutrients1', 'nutrients2', 'env_travel-time', 'basins')
taxons <- c(unique(data_uni_tax$taxName), 'unknown')

col_list <- data.frame(col=col_vector[1:length(taxons)], taxon=sort(taxons))
#if (taxo=='groups3'){
#  col_list$col[col_list$taxon=='Mamiellales']='#000AA88'
#  col_list$col[col_list$taxon=='Pelagophyceae']='#000A000'
#}
#print(col_list$col)
#print(' ')
if (taxo=='groups3'){
  #col_list$col <- as.character(levels(col_list$col))[col_list$col]
  #print(col_list$col)
  #print(' ')
  col_list$col[col_list$taxon=='Cnidaria'] =col_vector[31]
  col_list$col[col_list$taxon=='Tunicata'] =col_vector[34]
  col_list$col[col_list$taxon=='other'] =col_vector[38]
  col_list$col[col_list$taxon=='Mamiellales']='#00AA88'
  col_list$col[col_list$taxon=='Pelagophyceae']='#008000'   
  col_list$col[col_list$taxon=='other Stramenopiles']='#FFC000'
#col_list$col <- as.factor(col_list$col)
}
print(col_list$col)

if (taxo=='MGT-v2'){
  col_list$col[col_list$taxon=='42_Pelagomonadales'] = col_vector[13]
}


saveRDS(col_list, paste('color_table_',taxo,'.rds',sep=''))
saveRDS(go_colors,  paste('color_table_go.rds',sep=''))

fractions <-c('GGZZ','SSUU',  'QQSS', 'KKQQ')
results_enriched <- rep(list(NULL), length(fractions))
results_not_enriched <- rep(list(NULL), length(fractions))
for(frac in fractions){
  results_enriched_int <- rep(list(NULL), length(names0))
  results_not_enriched_int <- rep(list(NULL), length(names0))
  names(results_enriched_int) <- names0
  for (i in 1:length(comp_list)){
    components <- comp_list[[i]]
    name <- names0[i]
    result <- barplot_tax(components = components, name=name, frac = frac, subset=subset)
    results_enriched_int[[name]] <- result[[1]]
    results_not_enriched_int[[name]] <- result[[2]]
  }
  results_enriched[[frac]] <- results_enriched_int
  results_not_enriched[[frac]] <- results_not_enriched_int
  #saveRDS(results_enriched, paste('Significant_uids_GO_annotation_', taxo, '_',frac, '.rds', sep=''))
}
saveRDS(results_enriched, paste('Significant_enriched_uids_GO_annotation_', taxo,'_',subset,'.rds', sep=''))
saveRDS(results_not_enriched, paste('Significant_uids_GO_annotation_', taxo,'_',subset,'.rds', sep=''))
