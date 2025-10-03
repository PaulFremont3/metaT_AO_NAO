find_lines<-function(cond0, typo){
  lines <- list()
  count <- 1
  passed <- NULL
  if (length(cond0)>0){
    for (i in 1:length(cond0)){
      if (!(cond0[i] %in% passed)){
        if (typo == 2 & i != 1){
          new_line <- c(cond0[i]-1, cond0[i])
        } else{
          new_line <- c(cond0[i])
        }
        passed <- append(passed, cond0[i])
        j=which(cond0==cond0[i])
        while(cond0[j]==cond0[j+1]-1 & j < length(cond0)){
          new_line <- append(new_line, cond0[j+1])
          passed <- append(passed, cond0[j+1])
          j=j+1
        }
        if (typo == 2){
          new_line <- append(new_line, cond0[j]+1)
        }
        lines[[count]]=new_line
        count <- count+1
      } else{
        next
      }
    }
  }
  return(lines)
}

extract_data <- function(type, frac){
  data <- NULL
  count <-1
  simk <- simka[[type]][[frac]]
  simk_g <- simka[['G']][[frac]]
  bc <- bc_unigenes[[frac]][[type]]
  bc_g <- bc_unigenes[[frac]][['G']]
  for (i in env_arctic$Station){
    for (j in env_arctic$Station){
      # if (i<j & !is.nan(Tmin1$Tmin[Tmin1$Station_1==i & Tmin1$Station_2==j])){
      sk <- simk$Distance[simk$Station_1==i & simk$Station_2==j]
      if (length(sk)==0){
        sk <- NA
      }
      sk_g <- simk_g$Distance[simk_g$Station_1==i & simk_g$Station_2==j]
      if (length(sk_g)==0){
        sk_g <- NA
      }
      ed <- env_dist[rownames(env_dist)==i ,rownames(env_dist)==j]
      if (length(ed)==0){
        ed <- NA
      }
      br_c <- bc[which(rownames(bc)==i), which(colnames(bc)==j)]
      if (length(br_c)==0){
        br_c <- NA
      }
      br_c_g <- bc_g[which(rownames(bc_g)==i), which(colnames(bc_g)==j)]
      if (length(br_c_g)==0){
        br_c_g <- NA
      }
      if (i<="155SUR" & j <= "155SUR"){
        co <-'blue'
      } else if (i<="155SUR" & j > "155SUR"){
        co <- 'red'
      } else if (i>"155SUR" & j <= "155SUR"){
        co <- 'red'
      }else if (i>"155SUR" & j > "155SUR"){
        co <- 'green'
      }
      dot <- c(env_arctic$Latitude[env_arctic$Station==i],
               env_arctic$Longitude[env_arctic$Station==i],
               sk,sk_g,ed,
               Tmin1$Tmin[Tmin1$Station_1==i & Tmin1$Station_2==j], i,j, co, br_c, br_c_g)
      if (length(dot)!=11){
        print( i)
        print(j)
      }
      data <- rbind(data,dot)
    }
  }
  data <- as.data.frame(data)
  colnames(data) <- c('la', 'lo',  'sk','sk_g','ed','tmin', 'st1', 'st2', 'col', 'bc', 'bc_g')
  data$sk <- as.numeric(data$sk)
  data$sk_g <- as.numeric(data$sk_g)
  data$ed <- as.numeric(data$ed)
  data$tmin <- as.numeric(data$tmin)
  data$la <- as.numeric(data$la)
  data$lo <- as.numeric(data$lo)
  data$st1 <- as.character(data$st1)
  data$st2 <- as.character(data$st2)
  data$col<- as.character(data$col)
  data$bc <- as.numeric(data$bc)
  data$bc_g <- as.numeric(data$bc_g)
  
  
  return(data)
}

cum_cor <- function(X,  range = c(0.03, 20), var) {
  source('mantel1.R')
  source('getPermuteMatrix.R')
  Y <- X %>% arrange(st1, st2) %>% filter(st1 < st2)
  data <- Y[!is.na(Y$ed) & !is.na(Y$tmin),]
  
  
  simka <- data[[var]]
  tmin  <- data$tmin
  
  ## Sequences of Tmin thresholds
  threshold_seq <- 10^(seq(from = log10(range[1]), to = log10(range[2]), 
                           length.out = n_bins))
  ## Correlation between X and Y for Y-distances below a given threshold
  treshold_cor <- function(thresh) { 
    imin <- tmin <= thresh
    if (sum(imin)>5){
      corel <- cor(simka[imin], tmin[imin], method = "spearman", use = 'pairwise.complete.obs')
    }else{
      corel <- NA
    }
    return(corel)
  }
  threshold_mantel <- function(thresh){
    u <- matrix(unlist(X$tmin), ncol=sqrt(dim(X)[1]))
    v <- matrix(unlist(X[[var]]), ncol=sqrt(dim(X)[1]))
    v[u>thresh] <- NA
    if(sum(!is.na(u) & !is.na(v), na.rm=T)<20){
      sig <- 1
    } else{
      x <- mantel1(u, v, na.rm = T, method = 'spearman')
      sig <- x$signif
    }
    return(sig)
  }
  treshold_sd <- function(thresh) { 
    imin <- tmin <= thresh
    cors <- NULL
    for (i in 1:30){
      sub <- sample(1:length(tmin[imin]), 0.5*length(tmin[imin]))
      s <- simka[imin]
      t <- tmin[imin]
      c <- cor(s[sub], t[sub], method = "spearman")
      cors <- append(cors, c)
    }
    return(sd(cors, na.rm = T)/sqrt(30))
  }
  ## Compute correlation for every threshold value and store result in a data.frame
  tu <-tibble(thresh = threshold_seq,
         cor    = vapply(threshold_seq, treshold_cor, numeric(1)), 
         sd     = vapply(threshold_seq, treshold_sd, numeric(1)),
         mt     = vapply(threshold_seq, threshold_mantel, numeric(1)),
  )
  return(tu)
}

cum_cor_ed <- function(X,  range = c(0.03, 20), var) {
  source('mantel1.R')
  source('getPermuteMatrix.R')
  Y <- X %>% arrange(st1, st2) %>% filter(st1 < st2)
  data <- Y[!is.na(Y$ed) & !is.na(Y$tmin),]
  
  
  simka <- data[[var]]
  tmin  <- data$tmin
  ed <- data$ed
  
  ## Sequences of Tmin thresholds
  threshold_seq <- 10^(seq(from = log10(range[1]), to = log10(range[2]), 
                           length.out = n_bins))
  ## Correlation between X and Y for Y-distances below a given threshold
  treshold_cor <- function(thresh) { 
    imin <- tmin <= thresh
    if (sum(imin)>5){
      corel <- cor(simka[imin], ed[imin], method = "spearman", use = 'pairwise.complete.obs')
    }else{
      corel <- NA
    }
    return(corel)
  }
  threshold_mantel <- function(thresh){
    tm <- matrix(unlist(X$tmin), ncol=sqrt(dim(X)[1]))
    u <- matrix(unlist(X$ed), ncol=sqrt(dim(X)[1]))
    v <- matrix(unlist(X[[var]]), ncol=sqrt(dim(X)[1]))
    v[tm>thresh] <- NA
    if(sum(!is.na(u) & !is.na(v), na.rm=T)<20){
      sig <- 1
    } else{
      x <- mantel1(u, v, na.rm = T, method = 'spearman')
      sig <- x$signif
    }
    return(sig)
  }
  treshold_sd <- function(thresh) { 
    imin <- tmin <= thresh
    cors <- NULL
    for (i in 1:30){
      sub <- sample(1:length(tmin[imin]), 0.5*length(tmin[imin]))
      s <- simka[imin]
      t <- ed[imin]
      c <- cor(s[sub], t[sub], method = "spearman")
      cors <- append(cors, c)
    }
    return(sd(cors, na.rm = T)/sqrt(30))
  }
  ## Compute correlation for every threshold value and store result in a data.frame
  tu <-tibble(thresh = threshold_seq,
              cor    = vapply(threshold_seq, treshold_cor, numeric(1)), 
              sd     = vapply(threshold_seq, treshold_sd, numeric(1)),
              mt     = vapply(threshold_seq, threshold_mantel, numeric(1)),
  )
  return(tu)
}
