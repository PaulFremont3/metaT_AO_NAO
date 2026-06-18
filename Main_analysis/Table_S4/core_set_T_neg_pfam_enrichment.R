fraction='GGZZ'
taxo=commandArgs(trailingOnly = T)[1]
data_uni_tax=readRDS(paste('data_uni_',taxo,'.rds', sep=''))
#print(colnames(data_uni_tax))
if (taxo=='MGT-v2'){
  colnames(data_uni_tax)=c("taxName","geneID")
}
# 3. Map taxonomy (uid ↔ geneID)
tax_map <- data_uni_tax[, c("geneID", "taxName")]
df_tax=readRDS(paste('core_set_T_neg_',taxo,'.rds', sep=''))

if (taxo=='groups3'){
  subts <- c("Mamiellales", "Bacillariophyta", "Pelagophyceae", "Phaeocystales")
} else{
  subts <- c('3_Phaeocystales','439_Bathycoccaceae', '79_Bathycoccaceae', '195_Cymatosiraceae', '76_root', '4_Haptista',
                '42_Pelagomonadales')
}

# keep only relevant mapping
tax_map_sub <- data_uni_tax[data_uni_tax$taxName %in% subts, c("geneID", "taxName")]

all_pfam <- list()

for (arg in 1:10) {
  file <- paste("subset_metat_", fraction, "/", "matouUK", fraction, "_", arg, ".rds", sep = "")
  
  if (file.exists(file)) {
    matou <- readRDS(file)
    colnames(matou) <- c("geneID", "pfam", "Evalue")
    matou=data.frame(matou)
    # keep only genes with taxonomy of interest
    idx <- match(matou$geneID, tax_map_sub$geneID)
    keep <- !is.na(idx)
    
    if (any(keep)) {
      tmp <- matou[keep, c("geneID", "pfam")]
      all_pfam[[length(all_pfam) + 1]] <- tmp
    }
  }
}

# combine
pfam_all <- do.call(rbind, all_pfam)

pfam_all <- pfam_all[!duplicated(pfam_all[, c("geneID", "pfam")]), ]
pfam_all$taxName <- tax_map_sub$taxName[match(pfam_all$geneID, tax_map_sub$geneID)]
# background genes
bg_genes <- unique(pfam_all$geneID)

# test genes (your df_tax)
test_genes <- intersect(df_tax$uid, bg_genes)

pfam_bg <- aggregate(geneID ~ pfam, data = pfam_all,
                     function(x) length(unique(x)))
colnames(pfam_bg)[2] <- "n_bg"

pfam_test <- pfam_all[pfam_all$geneID %in% test_genes, ]

pfam_test_counts <- aggregate(geneID ~ pfam, data = pfam_test,
                              function(x) length(unique(x)))
colnames(pfam_test_counts)[2] <- "n_test"

pfam_stats <- merge(pfam_bg, pfam_test_counts, by = "pfam", all.x = TRUE)
pfam_stats$n_test[is.na(pfam_stats$n_test)] <- 0

N_bg <- length(bg_genes)
N_test <- length(test_genes)

pfam_stats$pval <- sapply(1:nrow(pfam_stats), function(i) {
  a <- pfam_stats$n_test[i]
  b <- N_test - a
  c <- pfam_stats$n_bg[i] - a
  d <- N_bg - (a + b + c)
  
  fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
})

pfam_stats$padj <- p.adjust(pfam_stats$pval, method = "hochberg")

saveRDS(pfam_stats, paste('pfam_enrichment_core_set_T_neg_',taxo,'.rds', sep=''))
saveRDS(pfam_all, paste('pfam_baseline_BMPP_',taxo,'.rds', sep=''))

