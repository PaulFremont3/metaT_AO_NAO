#setwd("~/Groups_metaT/Distances")
fracs =c('GGZZ')#, 'QQSS', 'KKQQ','SSUU')
arctic_stations <- paste(156:196, 'SUR', sep='')
atlantic_stations <- paste(143:155, 'SUR', sep='')
for (frac in fracs){
  if (frac%in% c('GGZZ', 'KKQQ')){
    taxis <- c('groups3','Mamiellales', 'Ciliophora','Dinophyceae', 'Bacillariophyta', 'Phaeocystales','Pelagophyceae', 'Bacteria', 'unknown')
  }  else if (frac %in% c('QQSS', 'SSUU')){
    taxis <- c('groups3','Hexanauplia', 'unknown', 'Insecta', 'other Opisthokonta', 'Cnidaria')
  }
  
  for (tax in taxis){
    test_uniT <- readRDS(paste('pfams_station_table_uni_T_',tax,'_allpf_',frac,'.rds', sep=''))
    test_uniG <- readRDS(paste('pfams_station_table_uni_G_',tax,'_allpf_',frac,'.rds', sep=''))
    
    w.test <- function(x){
      me_at <- median(as.numeric(x[names(x) %in% atlantic_stations]))
      me_ar <- median(as.numeric(x[names(x) %in% arctic_stations]))
      p<-wilcox.test(as.numeric(x[names(x) %in% atlantic_stations]), as.numeric(x[names(x) %in% arctic_stations]))$p.value
      return(c(me_at, me_ar, p))
    }
    
    tests_T <- t(apply(test_uniT, 1, FUN = w.test))
    tests_T[,3] <- p.adjust(tests_T[,3], method = 'BH')
    tests_T <- cbind(tests_T, tests_T[,1]-tests_T[,2])
    colnames(tests_T) <- c('median_atlantic_T', 'median_arctic_T', 'p.val_w_BH_T', 'diff.btw_T')
    
    tests_G <- t(apply(test_uniG, 1, FUN = w.test))
    tests_G[,3] <- p.adjust(tests_G[,3], method = 'BH')
    tests_G <- cbind(tests_G, tests_G[,1]-tests_G[,2])
    colnames(tests_G) <- c('median_atlantic_G', 'median_arctic_G', 'p.val_w_BH_G', 'diff.btw_G')
    
    tests <- cbind(tests_T, tests_G)
    
    saveRDS(tests, paste("unigenes_",tax,"_mTG_",frac,".rds", sep=''))
  }
}
