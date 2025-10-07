#!/bin/env/usr/env Rscript
#library("readxl")

GO_table <- read_excel('GO_table.xlsx', col_names = F)

count_GO_hierarchy <- function(u){
  return(length(strsplit(u, split=";")[[1]])-1)
}

max_GO_hierarchy <- function(u, GO_tab){
  GOS <- strsplit(GO_tab$X__4[GO_tab$X__1==u], split=";")[[1]]
  GOS <- GOS[GOS!=u]
  vec <- NULL
  vec0 <- NULL
  for (go in GOS){
    vec <- append(vec,GO_tab$X__5[GO_tab$X__1==go] )
    vec0 <- append(vec0, GO_tab$X__2[GO_tab$X__1==go])
  }
  GOs_hierarchy <- GOS[order(vec, decreasing = F)]
  GO_hierarc <- paste(GOs_hierarchy, collapse = ';')
  GOs_hierarchy0 <- vec0[order(vec, decreasing = F)]
  GO_hierarc0 <- paste(GOs_hierarchy0, collapse = ';')
  return(c(u,GO_tab$X__2[GO_tab$X__1==u], GO_hierarc,  GO_hierarc0, paste(sort(vec), collapse = ';'),max(vec), max(vec)+1 ))
}

GO_table$X__5 <- sapply(GO_table$X__4, count_GO_hierarchy)

GO_table1 <- GO_table[order(GO_table$X__5, decreasing = F),]
GO_table2 <- cbind(GO_table1$X__1[1:3], GO_table1$X__2[1:3] , rep('none', 3), rep('none', 3), rep(0,3), rep(0,3), rep(0,3))

unique_lengths = unique(GO_table1$X__5)
for (i in unique_lengths[2:length(unique_lengths)]){
  test <- lapply(GO_table1$X__1[GO_table1$X__5==i], max_GO_hierarchy, GO_tab=GO_table1)
  test1 <- unlist(test)
  test2 <- matrix(test1, ncol=7, byrow = T)
  
  GO_table1$X__5[GO_table1$X__5==i]= test2[,7]
  GO_table1$X__5 <- as.numeric(GO_table1$X__5 )
  GO_table2 <- rbind(GO_table2, test2)
}
GO_table2[,7]<-as.numeric(GO_table2[,7])

pfam2go <- read.table('pfam2go_27_02_20.txt', header = F)

missing_GO <- unique(pfam2go$V3)[!(unique(pfam2go$V3) %in% GO_table2[,1])]
update <- c('GO:1902600', 'GO:0038023', 'GO:0031012', 'GO:0036376', 'GO:0065003', 'obsolete', 
            'GO:1904680', 'GO:0044403', 'GO:0044877', 'GO:0007052', 'GO:1990904', 'GO:0016197')
pfam2go$V3 <- as.character(levels(pfam2go$V3))[pfam2go$V3]
pfam2go_update <- pfam2go
for (i in 1:length(missing_GO)){
  pfam2go_update$V3[pfam2go_update$V3==missing_GO[i]]=update[i]
}
GO_table3 <- GO_table2[GO_table2[,1] %in% unique(pfam2go_update$V3),]
saveRDS(GO_table2, 'GO_table.rds')
saveRDS(GO_table3, 'GO_table_pfam.rds')
colnames(GO_table2) <- c('ID', 'NAME', 'PARENTS_ID', 'PARENTS_NAME', 'PARENTS_LEVELS',  'MAX_LEVEL_PARENT', 'LEVEL')
write.table(GO_table2,'GO_table_hierarchy.txt', row.names = F)
#write.table(pfam2go_update, 'pfam2go_update.tab', row.names = F, col.names = F)
