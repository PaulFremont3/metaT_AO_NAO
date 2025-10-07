library('tidyr')
library('dplyr')
library('RColorBrewer')
#setwd("~/Groups_metaT/Distances")
source('aldex_bis.R')
fracs =c('GGZZ', 'SSUU', 'QQSS', 'KKQQ')
arctic_stations <- paste(156:196, 'SUR', sep='')
atlantic_stations <- paste(143:155, 'SUR', sep='')
types =c('', 'uni_')
for (frac in fracs){
  taxo <- '_taxo_groups3'
  for (ty in types){
    mT <- readRDS(paste('pfams_station_table_',ty,'T',taxo,'_allpf_',frac,'_bis.rds', sep=''))
    mG <-  readRDS(paste('pfams_station_table_',ty,'G',taxo,'_allpf_',frac,'_bis.rds', sep=''))
    
    na <- match(names(mG), names(mT))
    cond <- !is.na(na) & na!=1
    na <- na[cond]
    na_G <- names(mG)
    na_G <- na_G[cond]
    
    if (frac%in% c('GGZZ', 'KKQQ')){
      taxis <- c('groups3','Mamiellales', 'Dinophyceae', 'Bacillariophyta', 'Phaeocystales','Pelagophyceae', 'Bacteria','Ciliophora', 'unknown')
    }  else if (frac %in% c('QQSS', 'SSUU')){
      taxis <- c('groups3','Hexanauplia','Bacillariophyta', 'unknown', 'Insecta', 'other Opisthokonta', 'Cnidaria')
    }
    for (tax in taxis){
      dat_G <- NULL
      dat_T <- NULL
      pfs <- NULL
      c=1
      for (i in na){
        Tm <- mT[[i]]
        Gm <- mG[[na_G[c]]]
        if (tax != 'groups3'){
          if (tax %in% rownames(Tm) & tax %in% rownames(Gm)){
            dat_T <- rbind(dat_T, Tm[rownames(Tm)==tax,])
            dat_G <- rbind(dat_G, Gm[rownames(Gm)==tax,])
            pf <- na_G[c]
            pfs <- append(pfs, pf)
          }
        } else{
          if (!is.null(Tm)){
            dat_T <- rbind(dat_T, apply(Tm, 2, sum))
            dat_G <- rbind(dat_G, apply(Gm, 2, sum))
            pf <- na_G[c]
            pfs <- append(pfs, pf)
          }
        }
        c=c+1
      }
      rownames(dat_T) <- pfs
      rownames(dat_G) <- pfs
      saveRDS(dat_T,paste('pfams_station_table_',ty,'T_',tax,'_allpf_',frac,'.rds', sep=''))
      saveRDS(dat_G,paste('pfams_station_table_',ty,'G_',tax,'_allpf_',frac,'.rds', sep=''))
    }
  }
}
