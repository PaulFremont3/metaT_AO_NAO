library('data.table')

data_uni <- readRDS('data_uni_groups3.rds')
data_uni<- data_uni[!duplicated(data_uni$geneID),]
fraction <- commandArgs(trailingOnly = T)[1]
groups <-c('Hexanauplia', 'Bacillariophyta','Pelagophyceae','Crysophyceae' ,'Coccosphaerales','Phaeocystales' ,
           'Ciliophora', 'Tunicata','Mamiellales', 'Dinophyceae', 'Bacteria', 'Streptophyta', 'Cnidaria','Rhizaria',
           'Insecta', 'Cryptophyta' , 'Euglenozoa', 'Amoebozoa','Craniata',
           'Fungi', 'Viruses','Archaea' ,'Eukaryota (unclassified)','root' ,paste('other ' , c('Haptophyta', 'Opisthokonta', 'Stramenopiles', 'Alveolata', 'Viridiplantae') , sep='' ))

for (arg in 1:792){
  expr_uds <- readRDS(paste("subset_metat_",fraction,'/expressed_unilist_',fraction,'_',arg,'.rds', sep=''))
  df <- readRDS(paste("subset_metat_",fraction,"/",fraction,"metaTnMetaG_",arg,".rds", sep=''))
  uids <- unique(df$UID)
  uids <- uids[expr_uds]
  data_uni0 <- data_uni[data_uni$geneID %in% uids,]
  uk_uids <- uids[is.na(match(uids, data_uni0$geneID))]
  data_uni0 <- rbind(data_uni0, data.frame(geneID=uk_uids, taxName='unknown'))
  data_uni0 <- data_uni0[order(data_uni0$geneID),]
  df <- df[df$UID %in% uids,]
  if (fraction=='GGZZ'){
    df <- df[df$Station != '191SUR',]
  }
  df$MetaT[is.na(df$MetaT)]<-0
  df$MetaG[is.na(df$MetaG)]<-0
  v <- rep(NA, length(uids))
  for (g in groups){
    v[data_uni0$taxName==g] <- g
  }
  v[is.na(v) & data_uni0$taxName!='unknown']<-'other'
  v[is.na(v)] <- 'unknown'
  saveRDS(v, paste('TaxID_groups3_', fraction,'_' ,arg, '.rds',sep=''))
}
