#!/bin/env/usr/env Rscript
library("gbm")
library("dismo")
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library('scales')
library('parallel')
library('bestglm')
library('FactoMineR')

fraction = commandArgs(trailingOnly = T)[2] # 'GGZZ', 'SSUU', 'QQSS', 'KKQQ'
arg = commandArgs(trailingOnly = T)[1]    
df <- readRDS(paste("subset_metat_",fraction,"/",fraction,"metaTnMetaG_",arg,".rds", sep=''))
uids <- readRDS(paste('subset_metat_', fraction,"/pre_unilist_",fraction,'_',arg,'.rds', sep=''))
df <- df[df$UID %in% uids,]

df_clr <- readRDS(paste("subset_metat_",fraction,"/",fraction,"metaTnMetaG_",arg,"_clr.rds", sep=''))
df_clr <- df_clr[df_clr$UID %in% uids,]

df$MetaT[is.na(df$MetaT)]<-0
df$MetaG[is.na(df$MetaG)]<-0

matou <- readRDS(paste("subset_metat_",fraction,"/matou", fraction, '_', arg,'.rds', sep=''))
colnames(matou)<- c('geneID','pfam','Evalue')
matou <-as.data.frame(matou)

Lagr <- readRDS('../../Main_analysis/data/Lagr_distances.rds')
Lagr <- as.matrix(Lagr)

env_a <- read.table('../../Main_analysis/data/env_arctic_v3.txt', header = T)
env_a <- env_a[!(env_a$Station %in% c(142,201,205, 206, 208, 209, 210)),]
variables1 <- colnames(env_a)[4:13]
v <- env_a[,c(4:13)]
set.seed(1)
u <- PCA(v[,1:9], graph = F)
variables_pca <- u$ind$coord[, 1:3]

env_distances <- readRDS('../../Main_analysis/data/environmental_distances_arctic_0.rds')
env_distances <- cbind(env_a$Station, env_distances)

# rho from Erb and Notredame 2016
rho_erb <- function(x, y){
  rho = 1-var(x-y)/(var(x)+var(y))
  return(rho)
}

df$Station <- as.numeric(unlist(strsplit(df$Station, split = 'SUR')))
df_clr$Station <- as.numeric(unlist(strsplit(df_clr$Station, split = 'SUR')))
stations <- unique(df$Station)
if (fraction=='GGZZ'){
	good_stat <- stations[stations!=191]
} else{
	good_stat <- stations
}
good_stat <- good_stat[!(good_stat %in% c(142,201,205, 206, 
                                          208, 209, 210))]
df <- df[df$Station %in% good_stat,]
df_clr <- df_clr[df_clr$Station %in% good_stat,]

correlations <- function(uni){
  metat1 <- df$MetaT[df$UID==uni  & df$Station %in% good_stat]
  stat0 <- df$Station[df$UID==uni & df$Station %in% good_stat]
  metag0 <- df$MetaG[df$UID==uni  & df$Station %in% good_stat]
  condi_g <- metag0 != 0
  stat <- stat0[condi_g]
  stat2 <- paste(stat, collapse = '_')
  metatg_clr2 <- NA
  if (length(metat1[metat1!=0])>4 & max(stat)>155 & min(stat)<155){
    correls <- NULL
    correls_log <- NULL
    correls_pca <- NULL
    correls_clr <- NULL
    metat0 <- df$MetaT[df$UID==uni & df$Station %in% good_stat]
    
    metat <- metat1[condi_g]/metag0[condi_g]
    metag_c <- metag0[condi_g]
    metat_c <- metat1[condi_g]
    
    condit = metat1 !=0 | metag0 != 0
    metat2 = metat1[condit]
    metag = metag0[condit]
    dataf <- env_a[env_a$Station %in% stat,4:13]
    dataf <- cbind(dataf, metat)
    
    if (dim(dataf)[1]>4 & sum(metat, na.rm=T) != 0){
      nst = dim(dataf)[1]
      a= metag[stat<=155]
      b= metat2[stat<=155]
      c= metat2[stat>155]
      d= metag[stat>155]
      if ( median(a)>0 &  median(b)>0 &  median(c)>0 & median(d)>0){
        tr <- F
        count <- 0
        while (tr != T){
          sple <- sample(stat<=155)
          e= metag[sple]
          f= metat2[sple]
          g= metat2[!sple]
          h= metag[!sple]
          tr <- median(e)>0 &  median(f)>0 &  median(g)>0 & median(h)>0
          count <- count+1
          if (count >= 20){
            tr <- T
          }
        }
        ratio_t <- median(c)/median(b)
        ratio_g <- median(d)/median(a)
        if (count >= 20){
          ratio_t_r <- NA
          ratio_g_r <- NA
        } else{
          ratio_t_r <- median(g)/median(f)
          ratio_g_r <- median(h)/median(e)
        }
      } else{
        ratio_t <- NA
        ratio_g <- NA
        ratio_t_r <- NA
        ratio_g_r <- NA
      }
      
      a= metag_c[stat<=155]
      b= metat_c[stat<=155]
      c= metat_c[stat>155]
      d= metag_c[stat>155]
      if ( median(a)>0 &  median(b)>0 &  median(c)>0 & median(d)>0){
        tr <- F
        count <- 0
        while (tr != T){
          sple <- sample(stat<=155)
          e= metag_c[sple]
          f= metat_c[sple]
          g= metat_c[!sple]
          h= metag_c[!sple]
          tr <- median(e)>0 &  median(f)>0 &  median(g)>0 & median(h)>0
          count <- count+1
          if (count >= 20){
            tr <- T
          }
        }
        ratio_tc <- median(c)/median(b)
        ratio_gc <- median(d)/median(a)
        if (count >= 20){
          ratio_tc_r <- NA
          ratio_gc_r <- NA
        } else{
          ratio_tc_r <- median(g)/median(f)
          ratio_gc_r <- median(h)/median(e)
        }
      } else{
        ratio_tc <- NA
        ratio_gc <- NA
        ratio_tc_r <- NA
        ratio_gc_r <- NA
      }
      
      metag_s <- sample(metag0)  
      cor_g <- cor(metat0, metag0, method = 'pearson')
      cor_g_s <- cor(metat0, metag_s, method = 'pearson')
      correls <- append(correls, cor_g)
      correls <- append(correls, cor_g_s)
      
      Infs <- c(Inf, -Inf)
      metat_clr <- df_clr$MetaT[df_clr$UID==uni  & df_clr$Station %in% good_stat]
      metag_clr <- df_clr$MetaG[df_clr$UID==uni  & df_clr$Station %in% good_stat]
      conda <- !is.na(metat_clr) & !is.na(metag_clr) & !is.nan(metat_clr) & !is.nan(metag_clr) & !(metat_clr %in% Infs) & !(metag_clr %in% Infs)
      
      if (sum(conda)>4){
        stat_clr <- stat0[conda]
        stat2_clr <- paste(stat_clr, collapse='_')
        nst_clr=length(stat_clr)
        metat_clr <- metat_clr[conda]
        metag_clr <- metag_clr[conda]
        metag_clr_s <- sample(metag_clr)
        cor_g_clr <- rho_erb(metat_clr, metag_clr)
        cor_g_clr_s <- rho_erb(metat_clr, metag_clr_s)
      }else{
        cor_g_clr <- NA
        cor_g_clr_s <- NA
        stat_clr <- stat0[conda]
        stat2_clr <- paste(stat_clr, collapse='_')
        nst_clr=length(stat_clr)
      }
      
      metag_s_log <- sample(log10(metag0+1))  
      cor_g_log <- cor(log10(metat0+1), log10(metag0+1), method = 'spearman')
      cor_g_s_log <- cor(log10(metat0+1), metag_s_log, method = 'spearman')
      correls_log <- append(correls_log, cor_g_log)
      correls_log <- append(correls_log, cor_g_s_log)
      
      
      lgr <- NULL
      ev <- NULL
      mt <- NULL
      mt_log <- NULL
      # log_mt <- NULL
      for (st in stat){
        inde <- match(st, env_distances[,1])
        indl <- match(st, Lagr[,1])
        indt <- match(st, stat)
        for (st1 in stat){
          if (as.integer(st) < as.integer(st1)){
            inde1 <- match(st1, env_distances[,1])
            indl1 <- match(st1, Lagr[,1])
            indt1 <- match(st1, stat)
            if (!(is.na(Lagr[indl, indl1+1])) | !(is.na(Lagr[indl1, indl+1]))){
              ev <- append(ev, as.numeric(env_distances[inde, inde1+1]))
              min_lgr <- min(Lagr[indl, indl1+1], Lagr[indl1, indl+1], na.rm = T)
              lgr <- append(lgr, min_lgr)
              ok <- match(min_lgr, c(Lagr[indl, indl1+1], Lagr[indl1, indl+1]))
              if (ok ==1){
                mt <- append(mt, metat[indt1]-metat[indt])
                mt_log <- append(mt_log, log10(metat[indt1]+1)-log10(metat[indt]+1))
              } else{
                mt <- append(mt, metat[indt]-metat[indt1])
                mt_log <- append(mt_log, log10(metat[indt]+1)-log10(metat[indt1]+1))
              }
            }
          }
        }
      }
      
      cor_e <- cor(mt, ev, method = 'pearson')
      cor_e_s <- cor(mt, sample(ev), method = 'pearson')
      
      correls <- append(correls, cor_e)
      correls <- append(correls, cor_e_s)
      
      cor_l <- cor(mt, lgr, method ='pearson')
      cor_l_s <- cor(mt, sample(lgr), method ='pearson')
      
      correls <- append(correls, cor_l)
      correls <- append(correls, cor_l_s)
      
      cor_e_log <- cor(mt_log, ev, method = 'spearman')
      cor_e_s_log <- cor(mt_log, sample(ev), method = 'spearman')
      
      correls_log <- append(correls_log, cor_e_log)
      correls_log <- append(correls_log, cor_e_s_log)
      
      cor_l_log <- cor(mt_log, lgr, method ='spearman')
      cor_l_s_log <- cor(mt_log, sample(lgr), method ='spearman')
      
      correls_log <- append(correls_log, cor_l_log)
      correls_log <- append(correls_log, cor_l_s_log)
      
      correls_clr <- NULL
      metatg_clr <- metat_clr-metag_clr
      for (i in 4:13){
        value <- NULL
        for (l in stat){
          value <- append(value, as.numeric(env_a[env_a[,1]==l, i]))
        }
        cr <- cor(metat, value)
        cr_s <- cor(metat, sample(value))
        correls <- append(correls , cr)
        correls <- append(correls, cr_s)
        
        cr_log <- cor(log10(metat+1), value, method ='spearman')
        cr_s_log <- cor(log10(metat+1), sample(value), method ='spearman')
        correls_log <- append(correls_log , cr_log)
        correls_log <- append(correls_log, cr_s_log)
        
        if (length(stat_clr)>4){
          value_ <- NULL
          for (l in stat_clr){
            value_<- append(value_, as.numeric(env_a[env_a[,1]==l, i]))
          }
          cr_clr <- cor(metatg_clr, value_)
          cr_clr_s <- cor(metatg_clr, sample(value_))
          correls_clr <- append(correls_clr , cr_clr)
          correls_clr <- append(correls_clr, cr_clr_s)
        } else{
          correls_clr <- rep(NA, 20)
        }
      }
      
      for (i in 1:3){
        value <- NULL
        for (l in stat){
          value <- append(value, as.numeric(variables_pca[env_a[,1]==l, i]))
        }
        cr_pca <- cor(log10(metat+1), value, method ='spearman')
        cr_s_pca <- cor(log10(metat+1), sample(value), method ='spearman')
        correls_pca <- append(correls_pca , cr_pca)
        correls_pca <- append(correls_pca, cr_s_pca)
      }

      
      
      p_val <- NULL
      p_val_r <- NULL
      p_val_1 <- NULL
      if (sum(metat[stat>155]) > 0 & sum(metat[stat<=155]) > 0){
        p_val <- wilcox.test(metat[stat>155], metat[stat<=155], exact=FALSE)$p.value
        spl <- sample(stat>155)
        p_val_r <- wilcox.test(metat[spl], metat[!spl], exact=FALSE)$p.value
        ratio_tg_arc <- median(metat[stat>155])
        ratio_tg_atl <- median(metat[stat<=155])
        if (median(metat[stat>155]) > median(metat[stat<=155])){
          p_val_1 <- 'arctic'
        } else{
          p_val_1 <- 'atlantic'
        }
      } else if (sum(metat[stat>155])==0 & sum(metat[stat<=155]) > 0){
        p_val <- 'atlantic'
        p_val_r <- NA
        p_val_1 <- NA
        ratio_tg_arc <- NA
        ratio_tg_atl <- NA
      } else if (sum(metat[stat>155])>0 & sum(metat[stat<=155])==0){
        p_val <- 'arctic'
        p_val_r <- NA
        p_val_1 <- NA
        ratio_tg_arc <- NA
        ratio_tg_atl <- NA
      }
      if (is.null(p_val)){
        p_val <- NA
        p_val_r <- NA
        p_val_1 <- NA
        ratio_tg_arc <- NA
        ratio_tg_atl <- NA
      }
      
      p_val_clr <- NULL
      p_val_r_clr <- NULL
      p_val_1_clr <- NULL
      metatg_clr1 <- metatg_clr[!is.na(metatg_clr)]
      metatg_clr2 <- paste(metatg_clr1, collapse='_')
      if (sum(metatg_clr1[stat_clr>155]!=0, na.rm = T) > 0 & sum(metatg_clr1[stat_clr<=155]!=0, na.rm = T)  > 0){
        p_val_clr <- wilcox.test(metatg_clr1[stat_clr>155], metatg_clr1[stat_clr<=155], exact=FALSE)$p.value
        spl <- sample(stat_clr>155)
        p_val_r_clr <- wilcox.test(metatg_clr1[spl], metatg_clr1[!spl], exact=FALSE)$p.value
        ratio_tg_arc_clr <- median(metatg_clr1[stat_clr>155], na.rm=T)
        ratio_tg_atl_clr <- median(metatg_clr1[stat_clr<=155], na.rm=T)
        if (median(metatg_clr1[stat_clr>155], na.rm=T) > median(metatg_clr1[stat_clr<=155], na.rm=T)){
          p_val_1_clr <- 'arctic'
        } else{
          p_val_1_clr <- 'atlantic'
        }
      } else if (sum(metatg_clr1[stat_clr>155]!=0, na.rm=T)==0 & sum(metatg_clr1[stat_clr<=155]!=0, na.rm=T)>0){
        p_val_clr <- 'atlantic'
        p_val_r_clr <- NA
        p_val_1_clr <- NA
        ratio_tg_arc_clr <- NA
        ratio_tg_atl_clr <- NA
      } else if (sum(metatg_clr1[stat_clr>155]!=0, na.rm=T)>0 & sum(metatg_clr1[stat_clr<=155]!=0, na.rm=T)==0){
        p_val_clr <- 'arctic'
        p_val_r_clr <- NA
        p_val_1_clr <- NA
        ratio_tg_arc_clr <- NA
        ratio_tg_atl_clr <- NA
      }
      if (is.null(p_val_clr)){
        p_val_clr <- NA
        p_val_r_clr <- NA
        p_val_1_clr <- NA
        ratio_tg_arc_clr <- NA
        ratio_tg_atl_clr <- NA
      }
      
      if (dim(dataf)[1] > 11){
        flag <- TRUE
        possibleError <- tryCatch(
          b_lm_model <- bestglm::bestglm(dataf, IC='AIC'),
          error=function(e) flag <<- FALSE
        )
        if (flag==T){
          lm_model <- b_lm_model$BestModel
          p_values <- summary(lm_model)$coefficients[,4]
          lm_mod_t <- T
        } else{
          lm_model <- list()
          lm_model$call = 'lm(formula = y ~ 1)' 
          lm_mod_t <- F
        }
      } else{
        lm_model <- list()
        lm_model$call = 'lm(formula = y ~ 1)' 
        lm_mod_t <- F
      }
      if (lm_model$call != 'lm(formula = y ~ 1)' & lm_mod_t == T ){
        beta_coefs <- QuantPsyc::lm.beta(lm_model)
        j = 0
        signif <- NULL
        coefs <- NULL
        # print(lm_model$call)
        variables <- rownames(as.data.frame(b_lm_model$BestModel$coefficients))[2:length(b_lm_model$BestModel$coefficients)]
        for ( p in p_values){
          if (p<0.001 & j!=0){
            signif <- append(signif, variables[j])
            coefs <- append(coefs, beta_coefs[j][[1]])
          }
          j =j+1
        }
        
        variables0<-NULL
        for (var in variables){
          variables0 <- paste(variables0, var, sep='_')
        }
        if (is.null(variables0)){
          variables0 <- NA
        }
        signif0 <- NA
        for (i in 1:length(signif)){
          signif0 <- paste(signif0, signif[i], sep='_')
        }
        if (is.null(signif0)){
          signif0 <- NA
        }
        coefs0 <- NA
        for (i in 1:length(signif)){
          coefs0 <- paste(coefs0, coefs[i], sep='_')
        }
        if (is.null(coefs0)){
          coefs0 <- NA
        }
        
        aic <- extractAIC(lm_model, scale=0)
        if (is.null(aic)){
          aic <- NA
        }
        r_2<- summary(lm_model)$adj.r.squared
        effect_size <- NULL
        c<-0
        for (u in variables1){
          if (u %in% signif0){
            if (correls[7+c] > 0.5 | correls[7+c] < -0.5){
              effect_size <- paste(effect_size, u, sep='_')
            }
          }
          c<-c+2
        }
        if (is.null(effect_size)){
          effect_size <- 'no'
        }
      } else{
        signif0 <- NA
        coefs0 <- NA
        aic <- c(NA,NA)
        effect_size <- NA
        variables0 <- NA
        r_2 <-  NA
	      
	#c=1
        #for (st in stat){
        #  if (metat[stat==st] != 0){
        #    stat0 <- paste(stat0, st, sep='_')
        #  }
        #  c=c+1
        #}
        #if (is.null(stat0)){
        #  stat0 <- NA
        #}	        
      }
      index_min <-NULL
      index_min <- which.min(matou$Evalue[matou$geneID==uni])
      if (length(index_min)>0){
        pfam <- as.character(matou$pfam[matou$geneID==uni][index_min])
        pfams <- unique(matou$pfam[matou$geneID==uni])
        pfam_list <- NULL
        for (pf in pfams){
          pfam_list <- paste(pfam_list, pf, sep=';' )
        }
      } else{
        pfam<- NA
        pfam_list <- NA
      }
      to_write <- c(uni,correls, signif0, coefs0,aic[2], r_2, effect_size,variables0,stat2, 
                    p_val, p_val_1, pfam, ratio_t, ratio_g, ratio_tc ,ratio_gc, 
                    ratio_tg_arc, ratio_tg_atl, p_val_r, ratio_t_r, ratio_g_r, ratio_tc_r ,ratio_gc_r, pfam_list, 
                    nst, correls_log, correls_pca, cor_g_clr, cor_g_clr_s, correls_clr, 
                    p_val_clr, p_val_1_clr, ratio_tg_arc_clr, ratio_tg_atl_clr, p_val_r_clr,
                    stat2_clr, nst_clr, metatg_clr2)
      #print(length(to_write))
      if (length(to_write)==112){
        write(length(to_write), 'erreur_metat_an.txt', append=T)
        write(arg, 'erreur_metat_an.txt', append=T)	
      }
      if (length(to_write)==112){
        correlations<-to_write
      } else{
        correlations<-NULL
      }
    } else{
      correlations<-NULL
    }
  } else{
    correlations<-NULL
  }
}

uids <- unique(df$UID)

#no_cores <-detectCores()
# Initiate cluster
#cl <- makeCluster(no_cores)
#clusterExport(cl=cl, varlist=c("df_clr","rho_erb","df", 'Lagr',  'env_distances','variables_pca', 
#'env_a', 'variables1', 'good_stat', 'arg', 'matou'))
cor_list <- lapply(uids,FUN= correlations)
#stopCluster(cl) 
if (!is.null(unlist(cor_list))){
  cor_list_final <- NULL
  for (u in cor_list){
    if (!is.null(u) & length(u)==112){
      cor_list_final<- append(cor_list_final, u)
    }
  }
  cor_list_final <- matrix(cor_list_final, ncol=112, byrow = T)
  write.table(cor_list_final,paste('metat_analysis_',fraction,'_',arg,'.txt', sep=''), col.names = F, row.names = F)
} else {
  write.table(" ",paste('no_metat_analysis_',fraction,'_',arg,'.txt', sep=''), col.names = F, row.names = F)
}
