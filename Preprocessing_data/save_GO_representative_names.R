
types <-c('T', 'G')
basins <- c('atlantic', 'arctic')
frcs <- c('GGZZ', 'QQSS', 'SSUU')
bis <- c('', '_bis')
GT <- readRDS('GO_table.rds')
for (b in bis){
for (frac in frcs){
  for (ty in types){
    for (basin in basins){
      v <- readRDS(paste('GO_representative_',ty,'_',basin,'_',frac,'_taxo_groups3',b,'.rds', sep=''))
      nams <- GT[match(v, GT[,1]),2]
      saveRDS(nams,paste('GO_representative_names_',ty,'_',basin,'_',frac,'_taxo_groups3',b,'.rds', sep=''))
    }
  }
}
}
