library('tidyr')
library('dplyr')
library('RColorBrewer')
library('ALDEx2')
#setwd("~/Groups_metaT/Distances")
source('aldex_bis.R')
fracs =commandArgs(trailingOnly = T)[1]
arctic_stations <- paste(156:196, 'SUR', sep='')
atlantic_stations <- paste(143:155, 'SUR', sep='')

for(frac in fracs){
  if (frac%in% c('GGZZ', 'KKQQ')){
    #taxis <- c('','_Mamiellales', '_Dinophyceae', '_Bacillariophyta', '_Phaeocystales','_Pelagophyceae', '_Bacteria', '_unknown', '_Ciliophora')
    taxis <- c('_Ciliophora')
  }  else if (frac %in% c('QQSS', 'SSUU')){
    taxis <- c('','_Hexanauplia', '_unknown', '_Insecta', '_other Opisthokonta', '_Cnidaria')
  }
  for (tax in taxis){
    mT <- readRDS(paste('pfams_station_table_T',tax,'_allpf_',frac,'.rds', sep=''))
    mG <- readRDS(paste('pfams_station_table_G',tax,'_allpf_',frac,'.rds', sep=''))
    mT <- mT[,colnames(mT) %in% arctic_stations | colnames(mT) %in% atlantic_stations ]
    mG <- mG[,colnames(mG) %in% arctic_stations | colnames(mG) %in% atlantic_stations ]
    
    mG_reads <- read.table('metaG_590_filters_stats.mappedReads')
    mT_reads <- read.table('metaT_581_filters_stats')
    
    f <- function(x){
      paste(strsplit(x, 'SUR')[[1]][1], 'SUR', sep='')
    }
    
    if (frac != 'KKQQ'){
      mT_reads <- mT_reads[grepl(frac, mT_reads$V1) & grepl('SUR', mT_reads$V1),1:2]
    } else{
      mT_reads <- mT_reads[grepl(frac, mT_reads$V1) & grepl('SUR', mT_reads$V1) | grepl('MMQQ', mT_reads$V1) & grepl('SUR', mT_reads$V1) ,1:2]
    }
    if (is.factor(mT_reads$V1)){
      mT_reads$V1 <- as.character(levels(mT_reads$V1))[mT_reads$V1]
    }
    mT_reads$V1 <- sapply(mT_reads$V1, FUN = f)
    mT_reads <- mT_reads[mT_reads$V1 %in% arctic_stations | mT_reads$V1 %in% atlantic_stations,]
    
    if (frac != 'KKQQ'){
      mG_reads <- mG_reads[grepl(frac, mG_reads$V1) & grepl('SUR', mG_reads$V1),1:2]
    }else{
      mG_reads <- mG_reads[grepl(frac, mG_reads$V1) & grepl('SUR', mG_reads$V1) | grepl('MMQQ', mG_reads$V1) & grepl('SUR', mG_reads$V1),1:2]
    }
    if (is.factor(mG_reads$V1)){
      mG_reads$V1 <- as.character(levels(mG_reads$V1))[mG_reads$V1]
    }
    mG_reads$V1 <- sapply(mG_reads$V1, FUN = f)
    mG_reads <- mG_reads[mG_reads$V1 %in% arctic_stations | mG_reads$V1 %in% atlantic_stations,]
    
    g <- function(x, n){
      round(x*n)
    }
    mT0 <- t(apply(mT, 1, g, n=mT_reads$V2[match(colnames(mT), mT_reads$V1)]))
    mG0 <- t(apply(mG, 1, g, n=mG_reads$V2[match(colnames(mG), mG_reads$V1)]))
    
    
    condis <- c('atlantic', 'arctic')
    conds_T <- condis[as.numeric(colnames(mT0) %in% arctic_stations)+1]
    
    set.seed(1)
    mt0_aldex_clr <- aldex.clr(mT0, conds = conds_T)
    set.seed(1)
    mt0_aldex <- aldex(mT0, conditions = conds_T)
    
    condis <- c('atlantic', 'arctic')
    conds_G <- condis[as.numeric(colnames(mG0) %in% arctic_stations)+1]
    print(conds_G)
    set.seed(1)
    mg0_aldex_clr <- aldex.clr(mG0, conds = conds_G)
    set.seed(1)
    mg0_aldex <- aldex(mG0, conditions = conds_G)
    
    
    to_co_g<- match(rownames(mt0_aldex_clr@analysisData[[1]]), rownames(mg0_aldex_clr@analysisData[[1]]))
    to_co_t <- match(rownames(mg0_aldex_clr@analysisData[[1]]), rownames(mt0_aldex_clr@analysisData[[1]]))
    na <- names(mg0_aldex_clr@analysisData)
    me0_aldex_clr <-  mg0_aldex_clr
    e_clr <- list()
    for (st in na){
      e_clr[[st]] <- mt0_aldex_clr@analysisData[[st]][to_co_t[!is.na(to_co_t)],]-mg0_aldex_clr@analysisData[[st]][to_co_g[!is.na(to_co_g)],]
    }
    
    me0_aldex_clr@analysisData <- e_clr
    set.seed(1)
    me0_aldex <- aldex_bis(me0_aldex_clr, conditions=conds_G)
    
    
    saveRDS(list(mt0_aldex_clr, mt0_aldex, mg0_aldex_clr, mg0_aldex, me0_aldex_clr, me0_aldex),paste('aldex',tax,'_mTG_',frac,'.rds', sep=''))
  }
}

# understanding ALDEx
mg0_aldex <- readRDS(paste('aldex_mTG_',frac,'.rds', sep=''))
mg0_aldex_clr <- mg0_aldex[[1]]
mg0_aldex0 <- mg0_aldex[[2]]
na <- names(mg0_aldex_clr@analysisData)
arc <- 1:12
atl <- 13:20
wij_arc <- NULL
dwij_arc <- NULL
for (i in arc){
  wij_arc = append(wij_arc, mg0_aldex_clr@analysisData[[na[i]]][1,])
}
wij_arc_m <- median(wij_arc)



wij_atl <- NULL
dwij_atl <- NULL
for (i in atl){
  wij_atl = append(wij_atl, mg0_aldex_clr@analysisData[[na[i]]][1,])
}
wij_atl_m <- median(wij_atl)

nw <- min(c(length(wij_arc), length(wij_atl)))
sp1 <- sample(wij_arc, nw)
sp2 <- sample(wij_arc, nw)
dwij_arc_m <- abs(sp1-sp2)

sp1_at <- sample(wij_atl, nw)
sp2_at <- sample(wij_atl, nw)
dwij_atl_m <- abs(sp1_at-sp2_at)

win.max <- apply(  rbind( dwij_arc_m, dwij_atl_m )  , 2 , max )
delta_w <- median(win.max)

diff_b <- sp1_at-sp1
delta_A <- median(diff_b)
effec <- (sp1_at-sp1)/win.max
delta_R <- median(effec)
# 
# mt0_aldex_clr <- mg0_aldex[[1]]
# to_co_g<- match(rownames(mt0_aldex_clr@analysisData[[1]]), rownames(mg0_aldex_clr@analysisData[[1]]))
# to_co_t <- match(rownames(mg0_aldex_clr@analysisData[[1]]), rownames(mt0_aldex_clr@analysisData[[1]]))
# me0_aldex_clr <-  mg0_aldex[[3]]
# e_clr <- list()
# for (st in na){
#     e_clr[[st]] <- mt0_aldex_clr@analysisData[[st]][to_co_t[!is.na(to_co_t)],]-mg0_aldex_clr@analysisData[[st]][to_co_g[!is.na(to_co_g)],]
# }
# me0_aldex_clr@analysisData <- e_clr
# 
# me0_aldex <- aldex_bis(me0_aldex_clr)

