#!/bin/env/usr/env Rscript
#library("readxl")
library('scales')
library('parallel')
library('dplyr')
library('stringr')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")
pfam2go <- read.table('pfam2go_27_02_20.txt', header = F)
GO_table <- readRDS('GO_table.rds')

fraction = commandArgs(trailingOnly = T)[1]
n_clusts = commandArgs(trailingOnly = T)[2]
taxo = commandArgs(trailingOnly = T)[3]
bis = commandArgs(trailingOnly = T)[4]

if (bis != '_bis'){
  bis <- ''
}
if (taxo != "0" ){
  tax_id = taxo
}else{
  tax_id = ''
}

n_clusts <- as.numeric(n_clusts)
if (n_clusts==1){
  clusts=''
} else if (n_clusts %in% c(2,5,6) & taxo != "0"){
  clusts <- paste(c('com', 'atl', 'arc'), tax_id, sep='_')
} else if (n_clusts %in% c(3,4) & taxo != "0"){
  clusts <- paste(c('trans', 'atl', 'arc', 'com'), tax_id, sep='_')
} else if (n_clusts %in% c(2,5,6) & taxo == "0"){
  clusts <- c('com', 'atl', 'arc')
} else if (n_clusts %in% c(3,4) & taxo == "0"){
  clusts <- c('trans', 'atl', 'arc', 'com')
} 

types <- c('T', 'uni_T', 'G', 'uni_G')
pfam_data <- rep(list(NULL), n_clusts+2)
for (clust in clusts){
  if (n_clusts!=1 & (n_clusts!=4) & (n_clusts!=6)){
    pfam_station_table <- readRDS(paste('pfams_station_table_T_',clust,'_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_uni <- readRDS(paste('pfams_station_table_uni_T_',clust ,'_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_G <- readRDS(paste('pfams_station_table_G_',clust,'_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_uni_G<- readRDS(paste('pfams_station_table_uni_G_',clust,'_allpf_',fraction,bis,'.rds', sep=''))
    pfam_data[[clust]] <- list(pfam_station_table , pfam_station_table_uni, pfam_station_table_G, pfam_station_table_uni_G)
    names(pfam_data[[clust]]) <- types
  } else if (n_clusts==1 & taxo=='0'){
    pfam_station_table <- readRDS(paste('pfams_station_table_T_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_uni <- readRDS(paste('pfams_station_table_uni_T_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_G <- readRDS(paste('pfams_station_table_G_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_uni_G <- readRDS(paste('pfams_station_table_uni_G_allpf_',fraction,bis,'.rds', sep=''))
    pfam_data[[1]] <- list(pfam_station_table , pfam_station_table_uni, pfam_station_table_G, pfam_station_table_uni_G)
    names(pfam_data[[1]]) <- types
  } else if (n_clusts==1 & taxo!='0'){
    pfam_station_table <- readRDS(paste('pfams_station_table_T_',taxo,'_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_uni <- readRDS(paste('pfams_station_table_uni_T_',taxo,'_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_G <- readRDS(paste('pfams_station_table_G_',taxo,'_allpf_',fraction,bis,'.rds', sep=''))
    pfam_station_table_uni_G <- readRDS(paste('pfams_station_table_uni_G_',taxo,'_allpf_',fraction,bis,'.rds', sep=''))
    pfam_data[[1]] <- list(pfam_station_table , pfam_station_table_uni, pfam_station_table_G, pfam_station_table_uni_G)
    names(pfam_data[[1]]) <- types
  } else if (n_clusts %in% c(4,6) & taxo!='0'){
    pfam_station_table <- readRDS(paste('pfams_station_table_T_',clust,'_allpf_',fraction,'_',n_clusts,bis,'.rds', sep=''))
    pfam_station_table_uni <- readRDS(paste('pfams_station_table_uni_T_',clust,'_allpf_',fraction,'_',n_clusts,bis,'.rds', sep=''))
    pfam_station_table_G <- readRDS(paste('pfams_station_table_G_',clust,'_allpf_',fraction,'_',n_clusts,bis,'.rds', sep=''))
    pfam_station_table_uni_G <- readRDS(paste('pfams_station_table_uni_G_',clust,'_allpf_',fraction,'_',n_clusts,bis,'.rds', sep=''))
    pfam_data[[clust]] <- list(pfam_station_table , pfam_station_table_uni, pfam_station_table_G, pfam_station_table_uni_G)
    names(pfam_data[[clust]]) <- types
  }
}
if (taxo=='0'){
  all_pfams <- rownames(pfam_station_table)
} else if (taxo!='0'){
  all_pfams <- names(pfam_station_table)
}
rm(pfam_station_table , pfam_station_table_uni, pfam_station_table_G, pfam_station_table_uni_G)

create_GO_list <- function(type){
  library('dplyr')
  library('stringr')
  if (n_clusts==4){
    n_clusts0 =  3
  } else{
    n_clusts0 = n_clusts
  }
  if (n_clusts0>1){
    single_GO_table <- rep(list(NULL), n_clusts0+2)
    names(single_GO_table) <- c('all', clusts)
    GO_data <- rep(list(single_GO_table), dim(GO_table)[1]+1)
    names(GO_data)<- c(GO_table[,1], 'Unknown')
  } else{
    single_GO_table <- list(NULL)
    names(single_GO_table) <- 'all'
    GO_data <- rep(single_GO_table, dim(GO_table)[1]+1)
    names(GO_data)<- c(GO_table[,1], 'Unknown')
  }
  count= 1
  for (pfam in all_pfams){
    GOs <- pfam2go$V3[pfam2go$V1==pfam]
    passed <- NULL
    if (length(GOs)>0){
      for (go in GOs){
        passed <- append(passed, go)
        for (clust in clusts){
          if (n_clusts0>1 & taxo != '0'){
            to_bind <- pfam_data[[clust]][[type]][[pfam]]
            if (!is.null(to_bind)){
              GO_data[[go]][[clust]] <- rbind(GO_data[[go]][[clust]], to_bind)
              GO_data[[go]][['all']] <- rbind(GO_data[[go]][['all']], to_bind)
            }
          } else if (n_clusts0==1 & taxo=='0'){
            GO_data[[go]] <- rbind(GO_data[[go]], pfam_data[[1]][[type]][which(rownames(pfam_data[[1]][[type]])==pfam),])
          } else if (n_clusts0==1 & taxo!='0'){
            to_bind <- pfam_data[[1]][[type]][[pfam]]
            if (!is.null(to_bind)){
              GO_data[[go]] <- rbind(GO_data[[go]], to_bind)
            }
          } else if (n_clusts0>1 & taxo == '0'){
            GO_data[[go]][[clust]] <- rbind(GO_data[[go]][[clust]], pfam_data[[clust]][[type]][which(rownames(pfam_data[[clust]][[type]])==pfam),])
            GO_data[[go]][['all']] <- rbind(GO_data[[go]][['all']], pfam_data[[clust]][[type]][which(rownames(pfam_data[[clust]][[type]])==pfam),])
          }
        }
        hierarchy <- strsplit(GO_table[GO_table[,1]==go, 3], ';')[[1]]
        for (go_h in hierarchy){
          if (go_h %in% passed | go_h =='none'){
            next
          } else{
            for (clust in clusts){
              if (n_clusts0>1 & taxo != '0'){
                to_bind <- pfam_data[[clust]][[type]][[pfam]]
                if (!is.null(to_bind)){
                  GO_data[[go_h]][[clust]] <- rbind(GO_data[[go_h]][[clust]], to_bind)
                  GO_data[[go_h]][['all']] <- rbind(GO_data[[go_h]][['all']], to_bind)
                }
              } else if (n_clusts0==1 & taxo=='0'){
                GO_data[[go_h]]<- rbind(GO_data[[go_h]], pfam_data[[1]][[type]][which(rownames(pfam_data[[1]][[type]])==pfam),])
              } else if (n_clusts0==1 & taxo!='0'){
                to_bind <- pfam_data[[1]][[type]][[pfam]]
                if (!is.null(to_bind)){
                  GO_data[[go_h]] <- rbind(GO_data[[go_h]], to_bind)
                }
              } else if (n_clusts0>1 & taxo == '0'){
                GO_data[[go_h]][[clust]] <- rbind(GO_data[[go_h]][[clust]], pfam_data[[clust]][[type]][which(rownames(pfam_data[[clust]][[type]])==pfam),])
                GO_data[[go_h]][['all']] <- rbind(GO_data[[go_h]][['all']], pfam_data[[clust]][[type]][which(rownames(pfam_data[[clust]][[type]])==pfam),])
              }
            }
          }
        }
      } 
    } else{
      for (clust in clusts){
        if (n_clusts0>1 & taxo != '0'){
          to_bind <- pfam_data[[clust]][[type]][[pfam]]
          if (!is.null(to_bind)){
            GO_data[['Unknown']][[clust]] <- rbind(GO_data[['Unknown']][[clust]], to_bind)
            GO_data[['Unknown']][['all']] <- rbind(GO_data[['Unknown']][['all']], to_bind)
          }
        } else if (n_clusts0==1 & taxo=='0'){
          GO_data[['Unknown']] <- rbind(GO_data[['Unknown']], pfam_data[[1]][[type]][which(rownames(pfam_data[[1]][[type]])==pfam),])
        } else if (n_clusts0==1 & taxo!='0'){
          to_bind <- pfam_data[[1]][[type]][[pfam]]
          if (!is.null(to_bind)){
            GO_data[['Unknown']] <- rbind(GO_data[['Unknown']],to_bind)
          }
        } else if (n_clusts0>1 & taxo == '0'){
          GO_data[['Unknown']][[clust]] <- rbind(GO_data[['Unknown']][[clust]], pfam_data[[clust]][[type]][which(rownames(pfam_data[[clust]][[type]])==pfam),])
          GO_data[['Unknown']][['all']] <- rbind(GO_data[['Unknown']][['all']], pfam_data[[clust]][[type]][which(rownames(pfam_data[[clust]][[type]])==pfam),])
        }
      }
    }
    write(count*100/length(all_pfams),paste('follow_GO_', fraction,'_',n_clusts,'_' ,type,'_' ,taxo,bis,'.txt', sep=''), append = F)
    count=count+1
  }
  saveRDS(GO_data, paste('GO_station_table_interm_',fraction,'_',n_clusts,'_' ,type,'_' ,taxo,bis,'.rds', sep=''))
  #GO_data <- readRDS(paste('GO_station_table_interm_',fraction,'_',n_clusts,'_' ,type,'_' ,taxo,'.rds', sep=''))
  go_names  <- names(GO_data)
  for (go in go_names){
    data <- GO_data[[go]]
    if (n_clusts0 != 1 & taxo != '0'){
      if (!is.null(data[[1]])){
        for (clust in names(single_GO_table)){
          data1 <- data[[clust]]
          if (!is.null(data1)){
            data0 <- as.data.frame(data[[clust]])
            v <- rownames(data0)
            for (i in 1:4){
              v <- stringr::str_replace("[:digit:]$", string = v, replace='')
            }
            data0$taxo <- v
            data0 <- data0 %>% group_by(taxo) %>% summarise_all(funs(sum))
            data0 <- as.data.frame(data0)
            row.names(data0) <- data0[,1]
            data0 <- data0[,-1]
            sums_func <- apply(data0, 1, sum)
            order_func <- order(sums_func, decreasing = T)
            data0 <- data0[order_func,]
            GO_data[[go]][[clust]] <- data0
          }
        }
      }
      write(which(go_names==go)*100/length(go_names),paste('follow_GO1_', fraction,'_',n_clusts,'_' ,type,'_' ,taxo,bis,'.txt', sep=''), append = F)
    } else if (n_clusts0==1 & taxo != '0'){
      if (!is.null(data)){
        data0 <- as.data.frame(data)
        v <- rownames(data0)
        for (i in 1:4){
          v <- stringr::str_replace("[:digit:]$", string = v, replace='')
        }
        data0$taxo <- v
        data0 <- data0 %>% group_by(taxo) %>% summarise_all(funs(sum))
        data0 <- as.data.frame(data0)
        row.names(data0) <- data0[,1]
        data0 <- data0[,-1]
        sums_func <- apply(data0, 1, sum)
        order_func <- order(sums_func, decreasing = T)
        data0 <- data0[order_func,]
        GO_data[[go]] <- data0
      }
      write(which(go_names==go)*100/length(go_names),paste('follow_GO1_', fraction,'_',n_clusts,'_' ,type,'_' ,taxo,bis,'.txt', sep=''), append = F)
    } else if (n_clusts0==1 & taxo == '0'){
      if (!is.null(data)){
        data0 <- as.data.frame(data)
        data0 <- apply(data, 2, sum)
        GO_data[[go]] <- data0
      }
      write(which(go_names==go)*100/length(go_names),paste('follow_GO1_', fraction,'_',n_clusts,'_' ,type,'_' ,taxo,bis,'.txt', sep=''), append = F)
    } else if (n_clusts0!=1 & taxo == '0'){
      if (!is.null(data[[1]])){
        for (clust in names(single_GO_table)){
          data1 <- data[[clust]]
          if (!is.null(data1)){
            data0 <- as.data.frame(data[[clust]])
            data0 <- apply(data0, 2, sum)
            GO_data[[go]][[clust]] <- data0
          }
        }
      }
      write(which(go_names==go)*100/length(go_names),paste('follow_GO1_', fraction,'_',n_clusts,'_' ,type,'_' ,taxo,bis,'.txt', sep=''), append = F)
    }
  }
  return(GO_data)
}
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c('fraction',   'n_clusts',  'types', 'taxo', 'pfam_data', 'all_pfams', 'clusts', 'GO_table', 'pfam2go','bis' ))
Go_data_all_stations <- parLapply(cl, X = types,fun = create_GO_list)
stopCluster(cl)


save_rds <- function(data, k){
  saveRDS(data, paste('GO_station_table_',fraction,'_',n_clusts,'_' ,types[k],'_' ,taxo,bis,'.rds', sep=''))
}
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c('fraction',   'n_clusts',  'types', 'taxo', 'bis'))
v <- clusterMap(cl, save_rds,Go_data_all_stations, c(1:4))
stopCluster(cl)

