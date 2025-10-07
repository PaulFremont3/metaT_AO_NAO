#!/bin/env/usr/env Rscript
library('treemap')
library('gplots')
library('FactoMineR')
library('stringr')


frac <-commandArgs(trailingOnly = T)[1]
n <- 792

matou <- readRDS(paste("subset_metat_",frac,"/",'matou', frac, '_', 1,'.rds', sep=''))
colnames(matou)<- c('geneID','pfam','Evalue')
matou <-as.data.frame(matou)
matou <- matou[with(matou, order(geneID, Evalue)),]
matou <- matou[!duplicated(matou$geneID),]
matou_baseline <- as.data.frame(table(matou$pfam, useNA = 'ifany'))

matou_allpf <- readRDS(paste("subset_metat_",frac,"/",'matouUK', frac, '_', 1,'.rds', sep=''))
colnames(matou_allpf)<- c('geneID','pfam','Evalue')
matou_allpf <-as.data.frame(matou_allpf)
matou_baseline_allpf <- as.data.frame(table(matou_allpf$pfam, useNA = 'ifany'))

matou_allpf0 <- readRDS(paste("subset_metat_",frac,"/",'matou', frac, '_', 1,'.rds', sep=''))
colnames(matou_allpf0)<- c('geneID','pfam','Evalue')
matou_allpf0 <-as.data.frame(matou_allpf0)
matou_baseline_allpf0 <- as.data.frame(table(matou_allpf0$pfam, useNA = 'ifany'))
local(
for (i in 2:n){
  matou <- readRDS(paste("subset_metat_",frac,"/",'matou', frac, '_', i,'.rds', sep=''))
  colnames(matou)<- c('geneID','pfam','Evalue')
  matou <-as.data.frame(matou)
  matou <- matou[with(matou, order(geneID, Evalue)),]
  matou <- matou[!duplicated(matou$geneID),]
  baseline <- as.data.frame(table(matou$pfam, useNA = 'ifany'))
  pfams_to_add <- baseline$Freq[match(matou_baseline$Var1, baseline$Var1)]
  pfams_to_add[is.na(pfams_to_add)]<-0
  matou_baseline$Freq <<- matou_baseline$Freq + pfams_to_add
  new_pfams <- baseline[!(baseline$Var1 %in% matou_baseline$Var1),]
  matou_baseline<<-rbind(matou_baseline, new_pfams)
  
  matou_allpf <- readRDS(paste("subset_metat_",frac,"/",'matouUK', frac, '_', i,'.rds', sep=''))
  colnames(matou_allpf)<- c('geneID','pfam','Evalue')
  matou_allpf <-as.data.frame(matou_allpf)
  baseline_allpf <- as.data.frame(table(matou_allpf$pfam, useNA = 'ifany'))
  pfams_to_add_allpf <- baseline_allpf$Freq[match(matou_baseline_allpf$Var1, baseline_allpf$Var1)]
  pfams_to_add_allpf[is.na(pfams_to_add_allpf)]<-0
  matou_baseline_allpf$Freq <<- matou_baseline_allpf$Freq + pfams_to_add_allpf
  new_pfams_allpf <- baseline_allpf[!(baseline_allpf$Var1 %in% matou_baseline_allpf$Var1),]
  matou_baseline_allpf<<-rbind(matou_baseline_allpf, new_pfams_allpf)
  
  matou_allpf0 <- readRDS(paste("subset_metat_",frac,"/",'matou', frac, '_', i,'.rds', sep=''))
  colnames(matou_allpf0)<- c('geneID','pfam','Evalue')
  matou_allpf0 <-as.data.frame(matou_allpf0)
  baseline_allpf0 <- as.data.frame(table(matou_allpf0$pfam, useNA = 'ifany'))
  pfams_to_add_allpf0 <- baseline_allpf0$Freq[match(matou_baseline_allpf0$Var1, baseline_allpf0$Var1)]
  pfams_to_add_allpf0[is.na(pfams_to_add_allpf0)]<-0
  matou_baseline_allpf0$Freq <<- matou_baseline_allpf0$Freq + pfams_to_add_allpf0
  new_pfams_allpf0 <- baseline_allpf0[!(baseline_allpf0$Var1 %in% matou_baseline_allpf0$Var1),]
  matou_baseline_allpf0<<-rbind(matou_baseline_allpf0, new_pfams_allpf0)
})
rm(matou)
saveRDS(matou_baseline, paste('matou_baseline_', frac, '.rds', sep=''))
saveRDS(matou_baseline_allpf, paste('matou_baseline_allpf_', frac, '.rds', sep=''))
saveRDS(matou_baseline_allpf0, paste('matou_baseline_allpf_noUK_', frac, '.rds', sep=''))
