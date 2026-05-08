
taxo=commandArgs(trailingOnly = T)[1]
u=readRDS(paste('../data/Significant_enriched_uids_GO_annotation_',taxo,'_1.rds', sep=''))
# 1. Filter dataset
df <- u$GGZZ$physical_clr[
  u$GGZZ$physical_clr$taxo == "MBPP" &
  u$GGZZ$physical_clr$cor_sign == "-" &
  u$GGZZ$physical_clr$vari == "T",
]

df$uid=as.numeric(df$uid)
# 2. Store unique GO terms
go_terms <- unique(df$GO)

data_uni_tax=readRDS(paste('data_uni_',taxo,'.rds', sep=''))
#print(colnames(data_uni_tax))
if (taxo=='MGT-v2'){
  colnames(data_uni_tax)=c("taxName","geneID")
}# 3. Map taxonomy (uid ↔ geneID)
# keep only relevant columns to avoid duplicates explosion
tax_map <- data_uni_tax[, c("geneID", "taxName")]

df_tax <- merge(df, tax_map, by.x = "uid", by.y = "geneID", all.x = TRUE)

saveRDS(df_tax, paste('core_set_T_neg_',taxo,'.rds', sep=''))
print(dim(df_tax))
print(dim(df))
# 4. Count unique instances (uid) per GO
go_counts <- aggregate(uid ~ GO, data = df_tax, function(x) length(unique(x)))
colnames(go_counts)[2] <- "n_unique"

# 5. Count unique instances per GO per taxon
go_taxo_counts <- aggregate(uid ~ taxName + GO, data = df_tax,
                            function(x) length(unique(x)))
colnames(go_taxo_counts)[3] <- "n_unique"
saveRDS(list(go_counts, go_taxo_counts),paste('core_set_T_neg_functions_',taxo,'.rds', sep=''))
