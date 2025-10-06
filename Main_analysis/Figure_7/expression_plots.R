library("ggplot2")
library('gridExtra')
library('tidyr')
library('dplyr')
library('RColorBrewer')




custom_circle_part <- function(x,y,r,d1,d2,nsteps=nsteps,...){
  rs <- seq(d1,d2,len=nsteps)
  xc <- x+r*cos(rs)
  yc <- y+r*sin(rs)
  data <- data.frame(x=c(x,xc,x), y=c(y,yc,y))
  return(data)
}



make_expression_plots <- function(frac, name, taxo, vari, type){
  expr_data <- readRDS(paste('../data/metat_analysis_',frac,'.rds', sep=''))
  signif_data <- readRDS(paste('../data/Significant_',type,'uids_GO_annotation_',taxo,'_1.rds', sep=''))
  color_table <- readRDS(paste('../data/color_table_',taxo,'.rds', sep=''))
  color_table_go <- readRDS('../data/color_table_go.rds')
  colnames(color_table) <- c('col', 'item')
  colnames(color_table_go) <- c('item', 'col')
  env_arctic <- read.table('../data/env_arctic_3.txt', header=T)
  env_arctic <- env_arctic[order(env_arctic$T, decreasing = T),]
  env_arctic <- env_arctic[!(env_arctic$Station %in% c('143', '146', '149', '191')),]
  
  sig_dat <- signif_data[[frac]][[name]]
  sig_dat <- sig_dat[sig_dat$vari_s==vari,]
  rem_spaces <- function(x){
    gsub(" ", "", x, fixed = TRUE)
  }
  if (is.factor(sig_dat$uid)){
    sig_dat$uid <- as.character(levels(sig_dat$uid))[sig_dat$uid] 
  }
  sig_dat$uid <- sapply(sig_dat$uid, rem_spaces)
  expressions <- expr_data[expr_data[,1] %in% sig_dat$uid, c(1,49, 110, 112)]
  expr_data_ <- expr_data[expr_data[,1] %in% sig_dat$uid,c(1,49, 110, 112)]
  n1 <- dim(expressions)[1]
  dup_uids0 <- sig_dat$uid[duplicated(sig_dat$uid)]
  dup_uids <- sig_dat$uid[duplicated(sig_dat$uid)]
  while (length(dup_uids)>0){
    expressions <- rbind(expressions, expr_data_[expr_data_[,1] %in% dup_uids,])
    dup_uids <- dup_uids[duplicated(dup_uids)]
  }
  n2 <- dim(expressions)[1]
  expressions$taxo <- sig_dat$taxo[match(expressions$V1, sig_dat$uid)]
  if (is.factor(expressions$taxo)){
    expressions$taxo <- as.character(levels(expressions$taxo))[expressions$taxo]
  }
  expressions$go <- NULL
  expressions$go <- sig_dat$GO[match(expressions$V1, sig_dat$uid)]
  if (sum(duplicated(sig_dat$uid))>0){
    expressions$go[(n1+1):n2] <- sig_dat$GO[duplicated(sig_dat$uid)]
  }
  count_st <- function(x){
    length(strsplit(as.character(x[2][[1]]), split = '_')[[1]])
  }
  expressions$nst <- apply(expressions, 1, count_st)
  expressions <- expressions[with(expressions, order(taxo, go, nst, decreasing = T)),]
  
  
  data_to_plot <- matrix(0,nrow = dim(expressions)[1], ncol = 54)
  colnames(data_to_plot) <- 143:196
  for (i in 1:dim(expressions)[1]){
    stats <- strsplit(as.character(expressions[i, 3]), split = '_')[[1]]
    exprs <- as.numeric(strsplit(as.character(expressions[i, 4]), split = '_')[[1]])
    data_to_plot[i, match(stats,colnames(data_to_plot))] <- exprs
  }
  data_to_plot   
  
  sum_f <- function(x){
    return(sum(x!=0)>0)
  }
  data_to_plot0 <- data_to_plot[, apply(data_to_plot, 2, FUN = sum_f)]
  
  data_to_plot1 <- data_to_plot0
  for (j in 1:dim(expressions)[1]){
    max_ab <- max(data_to_plot0[j, data_to_plot0[j,]!=0], na.rm=T)
    min_ab <- min(data_to_plot0[j, data_to_plot0[j,]!=0], na.rm=T)
    
    OldRange = (max_ab - min_ab)
    newmax <- 1
    newmin <- 0.1
    NewRange = (newmax - newmin)
    data_to_plot1[j, data_to_plot0[j,]!=0] <- (data_to_plot1[j, data_to_plot0[j,]!=0]-min_ab)*NewRange/OldRange+newmin
  }
  
  
  extract_polygons <- function(data_to_plot, expressions_tab, type, title){
    pol_data <- list()
    pol_data_b <- list()
    pol_data_t <- list()
    #pol_data_4 <- list()
    functions <- NULL
    count=1
    colors_pol <- NULL
    items <- NULL
    if (is.null(dim(data_to_plot))){
      n_sp <- 1
      for (st in names(dat)){
        xo <- env_arctic$Station[env_arctic$Station==st]
        yo <- env_arctic$T[env_arctic$Station==st]
        
        xo_i <- which(env_arctic$Station==xo)
        
        if (data_to_plot[[st]] >0){
          polg <- custom_circle_part(x = xo_i, y = yo, d1=2*pi*j/n_sp,d2=2*pi*(j+1)/n_sp,r = data_to_plot[[st]],
                                     nsteps = 1000)
          polg_b <- custom_circle_part(x = xo_i*2, y = 1, d1=2*pi*j/n_sp,d2=2*pi*(j+1)/n_sp,r = data_to_plot[[st]],
                                       nsteps = 1000)
          
          polg_t <- c(j, xo_i, data_to_plot[[st]])
          
          pol_data[[count]] <- polg
          pol_data_b[[count]] <- polg_b
          pol_data_t[[count]] <- polg_t
          if (type=='taxo'){
            colo <- as.character(color_table$col[color_table$item== expressions_tab$taxo[1]])
            items <- append(items, as.character(expressions_tab$taxo[1]))
          } else if (type=='go'){
            colo <- as.character(color_table_go$col[color_table_go$item== as.character(expressions_tab$go[1])])
            items <- append(items, as.character(expressions_tab$go[1]))
          }
          colors_pol <- append(colors_pol, colo)
          
          count=count+1
        }
      }
      functions <- expressions_tab$go[1]
    } else{
      n_sp <- dim(data_to_plot)[1]
      for (j in 1:n_sp){
         for (i in 1:dim(data_to_plot)[2]){   
          xo <- env_arctic$Station[env_arctic$Station==colnames(data_to_plot)[i]]
          yo <- env_arctic$T[env_arctic$Station==colnames(data_to_plot)[i]]
          
          xo_i <- which(env_arctic$Station==xo)
          if (data_to_plot[j,i] >0){
            polg <- custom_circle_part(x = xo_i, y = yo, d1=2*pi*j/n_sp,d2=2*pi*(j+1)/n_sp,r = data_to_plot[j,i],
                                       nsteps = 1000)
            polg_b <- custom_circle_part(x = xo_i*2, y = 1, d1=2*pi*j/n_sp,d2=2*pi*(j+1)/n_sp,r = data_to_plot[j,i],
                                       nsteps = 1000)
            polg_t <- c(j, xo_i, data_to_plot[j,i])
            pol_data[[count]] <- polg
            pol_data_b[[count]] <- polg_b
            pol_data_t[[count]] <- polg_t
            if (type=='taxo'){
              colo <- as.character(color_table$col[color_table$item== expressions_tab$taxo[j]])
              items <- append(items, as.character(expressions_tab$taxo[j]))
            } else if (type=='go'){
              colo <- as.character(color_table_go$col[color_table_go$item== as.character(expressions_tab$go[j])])
              items <- append(items, as.character(expressions_tab$go[j]))
            }
            colors_pol <- append(colors_pol, colo)
            
            count=count+1
          }
        }
        functions <- append(functions, expressions_tab$go[j])
      }
    }
    data_polygons <- pol_data
    return(list(data_polygons, colors_pol, items, title, pol_data_b, pol_data_t, n_sp, functions))
  }
  
  conditions <- list(rep(T, dim(expressions)[1]),
                     !(expressions$taxo %in% c('unknown', 'other')))
  titles <- c('all', 'known')
  c=3
  for (t in unique(expressions$taxo)){
    if (t=='79_Bathycoccaceae'){
      condi <- expressions$taxo==t | expressions$taxo=='439_Bathycoccaceae'
    } else{
      condi <- expressions$taxo==t
    }
    if (sum(condi)>0){
      conditions[[c]] <- condi
      titles <- append(titles, t)
      if ('unknown' %in% expressions$go[condi]){
        condi1 = condi & expressions$go != 'unknown'
        if (sum(condi1)>0){
          conditions[[c+1]] <-condi1
          titles <- append(titles, paste(t, ' (without unknown)', sep=''))
          c=c+1
        }
      }
      if ('ribosome' %in% expressions$go[condi]){
        condi1 = condi & expressions$go != 'ribosome' & expressions$go != 'unknown'
        if (sum(condi1)>0){
          conditions[[c+1]] <-condi1
          titles <- append(titles, paste(t, ' (without ribosome)', sep=''))
          c=c+1
        }
      }
      if ('nucleus' %in% expressions$go[condi]){
        condi1 = condi & expressions$go != 'nucleus' & expressions$go != 'unknown'
        if (sum(condi1)>0){
          conditions[[c+1]] <-condi1
          titles <- append(titles, paste(t, ' (without nucleus)', sep=''))
          c=c+1
        }
      }
      c=c+1
    }
  }
  all_data <- list()
  count=1
  for (cond in conditions){
    dat <- data_to_plot1[cond, ]
    expre <- expressions[cond,]


    data_pol <- extract_polygons(data_to_plot = dat,expressions_tab = expre,
                                   type = "taxo",title =  titles[count])

    all_data[[count]] <- data_pol
    count=count+1
  }

  pdf(paste('Expression_map1_bis_',type,'uids_',frac,'_',
            name,'_', taxo,'_', vari,'.pdf', sep=''))
  for (dt in all_data){
    v <- scatter_plot_expr(data_c = dt[[1]], colors_pol = dt[[2]],
                    title = paste(dt[[4]],' ',vari ,sep=''), envo=env_arctic)
    print(v)
    plot(0,0, col='white', xaxt = 'n', yaxt='n', xlab = '', ylab = '', axes=F,
         main=paste('Taxonomy: ' ,dt[[4]],' ',vari, sep=''), cex.main=1)
    passed <- NULL
    leg <- NULL
    col_leg <- NULL
    counti=1
    for (it in dt[[3]]){
      if (!(it %in% passed)){
        leg <-append(leg, it)
        col_leg <- append(col_leg, dt[[2]][counti])
        passed <- append(passed, it)
      }
      counti=counti+1
    }
    if (length(leg)>30 & length(leg)<90){
      legend('topleft',legend = leg,
             fill=col_leg,
             box.lty=0 , ncol=4, cex=0.5)
    } else if (length(leg)>90){
      legend('topleft',legend = leg,
             fill=col_leg,
             box.lty=0 , ncol=5, cex=0.5)
    } else{
      legend('topleft',legend = leg,
             fill=col_leg,
             box.lty=0 ,ncol=2, cex=0.5)
    }
  }
  dev.off()
  pdf(paste('Expression_map1_ter_',type,'uids_',frac,'_',
            name,'_', taxo,'_', vari,'.pdf', sep=''), width=140)
  for (dt in all_data){
    v <- scatter_plot_expr_bis(data_c = dt[[5]], colors_pol = dt[[2]],
                           title = paste(dt[[4]],' ',vari ,sep=''), envo=env_arctic)
    print(v)

  }
  dev.off()
  pdf(paste('Expression_map1_bubble_',type,'uids_',frac,'_',
            name,'_', taxo,'_', vari,'.pdf', sep=''))
  for (dt in all_data){
    v <- scatter_plot_expr_bubble(data_c = dt[[6]], colors_pol = dt[[2]],
                           title = paste(dt[[4]],' ',vari ,sep=''), envo=env_arctic, n_genes=dt[[7]])
    print(v)
  }
  dev.off()
  pdf(paste('Expression_map1_bubble_bis_',type,'uids_',frac,'_',
            name,'_', taxo,'_', vari,'.pdf', sep=''))
  for (dt in all_data){
    v <- scatter_plot_expr_bubble_bis(data_c = dt[[6]], colors_pol = dt[[2]],
                           title = paste(dt[[4]],' ',vari ,sep=''), envo=env_arctic, n_genes=dt[[7]], funcs=dt[[8]])
    print(v)
  }
  dev.off()
}

scatter_plot_expr <- function(data_c, colors_pol, title, envo){
  plot(0,0, xlim=c(0, length(envo$Station)), ylim=c(min(envo$T), max(envo$T)), col=scales::alpha('white', 0),
       xlab='', ylab='Temperature', xaxt='n', main=title)
  axis(side = 1, at = c(1:length(envo$Station)), labels = paste(envo$Station, 'SRF', sep=''), las=2)
  i=1
  for (po in data_c){
    polygon(po$x, po$y, col=colors_pol[i], border = colors_pol[i])
    i=i+1
  }
}

scatter_plot_expr_bubble <- function(data_c, colors_pol, title, envo, n_genes){
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(0,0, ylim=c(1,n_genes), xlim=c(1, length(envo$Station)), col=scales::alpha('white', 0),
       xlab='', ylab='', xaxt='n', main=title, yaxt='n')
  axis(side = 1, at = c(1:length(envo$Station)), labels = paste(envo$Station, 'SRF', sep=''), las=2, cex=10)
  i=1
  for (po in data_c){
    points(po[2], po[1], col=scales::alpha(colors_pol[i], 0.5),pch=19, cex=(po[3]-0.1)*(3-0.5)/(1-0.1)+0.5)
    i=i+1
  }
  
}

scatter_plot_expr_bubble_bis <- function(data_c, colors_pol, title, envo, n_genes, funcs){
  par(mar=c(5.1, 1.1, 4.1, 15.1))
  plot(0,0, ylim=c(1,n_genes), xlim=c(1, length(envo$Station)), col=scales::alpha('white', 0),
       xlab='', ylab='', xaxt='n', main=title, yaxt='n')
  axis(side = 1, at = c(1:length(envo$Station)), labels = paste(envo$Station, 'SRF', sep=''), las=2, cex=10)
  i=1
  for (po in data_c){
    points(po[2], po[1], col=scales::alpha(colors_pol[i], 0.5),pch=19, cex=(po[3]-0.1)*(3-0.5)/(1-0.1)+0.5)
    i=i+1
  }
  axis(side=4, at = c(1:n_genes), labels =funcs, cex.axis=0.4, las=2)
}

scatter_plot_expr_bis <- function(data_c, colors_pol, title, envo){
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(0,0, xlim=c(1, 2*length(envo$Station)), ylim=c(0,2), col=scales::alpha('white', 0),
       xlab='', ylab='Temperature', xaxt='n', main=title)
  axis(side = 1, at = 2*c(1:length(envo$Station)), labels = paste(envo$Station, 'SRF', sep=''), las=2, cex=10)
  i=1
  for (po in data_c){
    polygon(po$x, po$y, col=colors_pol[i], border = colors_pol[i])
    i=i+1
  }
}


fractions <- c('GGZZ')
names <- c('physical_clr')
taxos <- c('groups3','MGT-v2')
variables <- list(c( 'T_-', 'T_+'))
types <- c('enriched_')
for (frac in fractions){
  for (name in names){
    varis <- variables[[which(names==name)]]
    for (taxo in taxos){
      for (v in varis){
        for (t in types){
          print(frac)
          make_expression_plots(frac = frac, name = name, taxo = taxo, vari = v, type = t)
          print(frac)
          print(name)
          print(taxo)
          print(v)
          print(t)
        }
      }
    }
  }
}


