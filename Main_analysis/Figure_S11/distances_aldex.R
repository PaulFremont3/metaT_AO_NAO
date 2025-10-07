library('tidyr')
library('dplyr')
library('RColorBrewer')
setwd("~/Groups_metaT/Distances")
fracs =c('GGZZ', 'SSUU', 'QQSS', 'KKQQ')
arctic_stations <- paste(155:196, 'SUR', sep='')
atlantic_stations <- paste(143:156, 'SUR', sep='')

for (frac in fracs){
  aldex_data <- readRDS(paste('aldex_mTG_',frac,'.rds', sep=''))
  
  clr_T0 <- aldex_data[[1]]@analysisData
  clr_G0 <- aldex_data[[3]]@analysisData
  clr_E0 <- aldex_data[[5]]@analysisData
  
  compute_dists <- function(clr_D0){
    nmonte <- length(clr_D0$X158SUR[1,])
    clr_D <- rep(list(NULL), nmonte)
    new_na <- as.character(sapply(names(clr_D0), FUN = function(x){strsplit(x, 'X')[[1]][2]}))
    for (i in 1:nmonte){
      for (st in names(clr_D0)){
        clr_D[[i]] <- cbind(clr_D[[i]], clr_D0[[st]][,i])
      }
      colnames(clr_D[[i]]) <- new_na
    }
    
    deuc <- function(x,y){
      v <- sqrt(sum((x-y)^2))
      return(v)
    }
    dist_euc <- function(x, d){
      d <- apply(d, MARGIN = 2, FUN = deuc, x=x)
    }
    
    data_d <- list()
    for (i in 1:nmonte){
      dmat <- apply(clr_D[[i]], MARGIN = 2, FUN = dist_euc, d=clr_D[[i]])
      data_d[[i]] <- dmat
    }
    return(data_d)
  }
  data_G <- compute_dists(clr_G0)
  data_E <- compute_dists(clr_E0)
  data_T <- compute_dists(clr_T0)
  
  arc_dists=NULL
  atl_dists=NULL
  arc_dists_E=NULL
  atl_dists_E=NULL
  arc_dists_T=NULL
  atl_dists_T=NULL
  for (j in 1:length(data_G)){
    arc_mat=data_G[[j]][9:20, 9:20]
    arc_dists = c(arc_dists, as.vector(arc_mat[upper.tri(arc_mat)]))
    atl_mat=data_G[[j]][1:8, 1:8]
    atl_dists = c(atl_dists, as.vector(atl_mat[upper.tri(atl_mat)]))
    
    arc_mate=data_E[[j]][9:20, 9:20]
    arc_dists_E = c(arc_dists_E, as.vector(arc_mate[upper.tri(arc_mat)]))
    atl_mate=data_E[[j]][1:8, 1:8]
    atl_dists_E = c(atl_dists_E, as.vector(atl_mate[upper.tri(atl_mat)]))
    
    arc_matt=data_T[[j]][9:20, 9:20]
    arc_dists_T = c(arc_dists_T, as.vector(arc_matt[upper.tri(arc_mat)]))
    atl_matt=data_T[[j]][1:8, 1:8]
    atl_dists_T = c(atl_dists_T, as.vector(atl_matt[upper.tri(atl_mat)]))
  }
  source('vioplot1.R')
  pdf('violinplot_dists.pdf')
  vioplot1(list(atl_dists,arc_dists), col=c('orange','darkblue') , names = c('Atlantic', 'Artcic'))
  vioplot1(list(atl_dists_E,arc_dists_E), col=c('orange','darkblue') , names= c('Atlantic', 'Artcic'))
  vioplot1(list(atl_dists_T,arc_dists_T), col=c('orange','darkblue') , names= c('Atlantic', 'Artcic'))
  dev.off()

  ratio <- function(a,b,i){
    a[[i]]/b[[i]]
  }
  
  data_r <- lapply(1:length(data_G), ratio, a=data_G,b=data_E)
  
  mean_r<- Reduce('+',data_r)/length(data_r)
  sd_r0 <- lapply(data_r, function(x){(x-mean_r)^2})
  sd_r <- sqrt(Reduce('+',sd_r0))*sqrt(1/length(data_r))
  
  env <- read.table('env_arctic_3.txt')
  env <- env[order(env$T, decreasing = F),]
  env$Station <- paste(env$Station, 'SUR', sep='')
  
  to_keep <- colnames(mean_r) %in% env$Station  
  mean_r <- mean_r[to_keep,to_keep]
  env <- env[env$Station %in% colnames(mean_r), ]
  
  d_ratios<- as_tibble(mean_r, rownames = "Station_1") %>% 
    gather(-Station_1, key = "Station_2", value = "ratiodist") %>% 
    ## Keep only upper triangular part
    arrange(Station_1, Station_2) %>% filter(Station_1 < Station_2)
  
  
  data <- NULL
  data_T <- NULL
  sigs <- NULL
  perc_d <- NULL
  arcs <- paste(155:210, 'SUR', sep='')
  for (i in 1:(dim(env)[1]-4)){
    sts <- env$Station[i:(i+4)]
    vals <- d_ratios$ratiodist[d_ratios$Station_1 %in% sts & d_ratios$Station_2 %in% sts]
    test <- wilcox.test(vals-1)
    sigs <- append(sigs, test$p.value)
    percs <- sum(sts %in% arcs)/length(sts)
    perc_d <- append(perc_d, percs)
    data <- append(data, median(vals))
    data_T <- append(data_T, median(env$T[env$Station %in% sts]))
  }
  sigs <- p.adjust(sigs, method = 'holm')
  cols <- colorRampPalette(c('blue', 'white', 'red'))(11)
  pdf(paste('distances_ratio_plot_pfams_aldex_', frac, '.pdf', sep=''))
  plot(data_T[sigs<=0.05], data[sigs<=0.05], xlab = 'Median temperature of the bin (°C)', ylab='Abundance-based distance / Expression-based distance',
       col=cols[perc_d*10+1][sigs<=0.05], pch=20, ylim = c(0.6, 1.1), cex=3, xlim = c(min(data_T), max(data_T)), cex.lab=1.5, cex.axis=1.5)
  points(data_T[sigs>0.05], data[sigs>0.05],  pch=23, ylim = c(0.6, 1.2), cex=3, bg=cols[perc_d*10+1][sigs>0.05], col=cols[perc_d*10+1][sigs>0.05])
  abline(h=1, lty=5, lwd=3)
  dev.off()
}

