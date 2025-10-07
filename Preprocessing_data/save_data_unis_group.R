library('data.table')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_2")

group <- commandArgs(trailingOnly = T)[1]
data_uni <- readRDS('taxID_uni_all.rds')
if (group == 'groups'){
  groups <-c('Hexanauplia', 'Bacillariophyta', 'Coccosphaerales', 'Ciliophora', 'Tunicata',
             'Mamiellales', 'Dinophyceae')
} else if (group =='groups3'){
  groups <-c('Hexanauplia', 'Bacillariophyta','Pelagophyceae','Crysophyceae' ,'Coccosphaerales','Phaeocystales' ,
           'Ciliophora', 'Tunicata','Mamiellales', 'Dinophyceae', 'Bacteria', 'Streptophyta', 'Cnidaria','Rhizaria',
           'Insecta', 'Cryptophyta' , 'Euglenozoa', 'Amoebozoa','Craniata',
           'Fungi', 'Viruses','Archaea' ,'Eukaryota (unclassified)','root' ,paste('other ' , c('Haptophyta', 'Opisthokonta', 'Stramenopiles', 'Alveolata', 'Viridiplantae') , sep='' ))
} else if (group == 'groups2'){
  groups <-c('Hexanauplia', 'Bacillariophyta','Pelagophyceae','Crysophyceae' ,'Coccosphaerales','Phaeocystales' ,
           'Ciliophora', 'Tunicata','Mamiellales', 'Dinophyceae', 'Bacteria', 'Streptophyta', 'Cnidaria','Rhizaria',
           'Insecta', 'Cryptophyta' , 'Euglenozoa', 'Amoebozoa','Craniata',
           'Fungi', 'Viruses')
}

v <- rep(NA, dim(data_uni)[1])
for (g in groups){
  if (grepl('root', g) | grepl('Eukaryota', g)){
    if (grepl('root', g)){
      v[data_uni$taxLineage==g & is.na(v)] <- g
    } else{
      gu <- 'root;cellular organisms;Eukaryota'
      v[data_uni$taxLineage==gu & is.na(v)] <- g
    }
  } else if (grepl('other', g) ){
    gu <- strsplit(g, split=' ')[[1]][2]
    v[grepl(gu, data_uni$taxLineage) & is.na(v)] <- g 
  } else{
    v[grep(g, data_uni$taxLineage)] <- g
  }
}
cond <- is.na(v)
u <-  data_uni$taxLineage[cond]
unique_u <- unique(u)
saveRDS(u, paste('other_lineages_',group ,'.rds', sep=''))
v[cond]<-'other'
write(sum(v=='other'), paste('not_in_groups_',group,'.txt', sep=''))
dat <- data.frame(taxName=v, geneID=data_uni$geneID)
dat$taxName <- as.character(levels(dat$taxName))[dat$taxName]
saveRDS(dat, paste('data_uni_',group,'.rds', sep=''))
