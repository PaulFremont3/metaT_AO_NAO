library('data.table')
data_uni <- fread('MATOU-v2.taxonomy.tsv',  select = c(1,3),header = T, sep='\t') # taxonomy file of MATOU
saveRDS(data_uni, 'taxID_uni.rds')
