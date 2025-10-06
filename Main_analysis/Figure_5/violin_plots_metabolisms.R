library('vioplot')
library('gplots')
dat=readRDS('../data/metat_analysis_GGZZ.rds')
unigenes=readLines('../data/unigenes_phate_GGZZ_0_0_0.txt')
data_expression=readRDS('../data/to_plot1_phate.rds')
names(data_expression)=c('delta0', 'phate-1', 'phate-2', 'mean0','atlanticity', 'mean1', 'delta1' )
data_expression[['delta1']]=-data_expression[['delta1']]
data_expression[['atlanticity']]=-data_expression[['atlanticity']]
new_min=-1
new_max=1
mi=min(data_expression[['atlanticity']])
mx=max(data_expression[['atlanticity']])
data_expression[['atlanticity']]=(data_expression[['atlanticity']]-mi)/(mx-mi)*(new_max-new_min)+new_min

raw_expr=read.table('../data/phate_training_expression_GGZZ_0_0_0.txt')



photosynthesis_pfams <- c('PF00504' ,'PF00127','PF11623','PF02276' , 'PF01716','PF00556' ,'PF00124', 'PF03967', 'PF00223', 'PF02605', 'PF07465', 'PF05479', 'PF00737', 'PF02532', 'PF02533', 'PF05151'
                          , 'PF02468', 'PF04725', 'PF01405', 'PF06514', 'PF07123', 'PF06596', 'PF06298' , 'PF00421', 'PF05969', 'PF13326', 'PF00796',
                          'PF02427', 'PF02507', 'PF03244', 'PF01701', 'PF01241', 'PF10657', 'PF11947')
nitrogen_pfams <- c('PF02665', 'PF01077','PF03460','PF08376'
                    ,'PF07732', 'PF00384',
                    'PF01568', 'PF04324', 'PF00174', 'PF07992', 'PF00355')
carbon_fixation_pfams <- c('PF00194', 'PF00016', 'PF00101', 'PF00485')#)#
iron_functions <- c( 'PF00111',  'PF13085', 'PF00258')
other_functions <- c('PF01061','PF00005', 'PF00854','PF07690', 'PF02421', 'PF04145', 'PF16983',
                     'PF02535', 'PF03842', 'PF01384','PF16974', 'PF00909')
temperature_pfams <- c('PF00360','PF03952','PF01786','PF06415',
                       'PF05971', 'PF01180', 'PF05351', 'PF05035', 'PF01346', 'PF03104', 'PF11999',
                       'PF00249', 'PF03259', 'PF02209', 'PF00241', 'PF00626', 'PF01119', 'PF01979',
                       'PF01964', 'PF02358','PF14249', 'PF02649',
                       'PF00313', 'PF05562', 'PF14169')
CET_photosynthesis = c('PF00124', 'PF00421', 'PF03967', 'PF11623')
other_photosynthesis = c('PF00127','PF02276' , 'PF01716',  'PF00223', 'PF02605', 'PF07465', 'PF05479', 'PF00737', 'PF02532', 'PF02533', 'PF05151'
                         , 'PF02468', 'PF04725', 'PF01405', 'PF06514', 'PF07123', 'PF06596', 'PF06298' , 'PF05969', 'PF13326', 'PF00796',
                         'PF02427', 'PF02507', 'PF03244', 'PF01701', 'PF01241', 'PF10657', 'PF11947')
light_harvest = c('PF00504','PF00556')
specific_nitrogen = c('PF02665', 'PF01077','PF03460')
nitrate_sensing = c('PF08376')
others_nitrogen = c('PF07732', 'PF00384',
                    'PF01568', 'PF04324', 'PF00174', 'PF07992', 'PF00355')
ferredoxin = c('PF13085', 'PF00111')
flavodoxin = c('PF00258')
general_transporters = c('PF01061','PF00005','PF07690' , 'PF00854') #
iron_transport = c('PF02421')
silicon_transporter = c('PF03842')
phosphate_transp = c('PF01384')
nitrate_transp = c('PF16974')
ammonium_transp = c('PF00909')
DMSP=c('PF16867')
cell_div=c('PF01111', 'PF00134','PF04042')
hydrolases=c('PF00703', 'PF03663')
glycolysis=c('PF00349', 'PF00365', 'PF00274', 'PF02887')
iron_storage=c('PF00210', 'PF01475')
copper_transport=c('PF04145')
molybdate_transport=c('PF16983')
zinc_transport=c('PF02535', 'PF16916')


list_all_pfams <- list(photosynthesis_pfams, nitrogen_pfams, carbon_fixation_pfams, iron_functions, other_functions, 
                        temperature_pfams)
names_pfams_list <- c('Photosynthesis', 'Nitrogen metabolism', 'Carbon fixation', 'Iron metabolism', 'Nutrient transporters',
                       'Temperature sensitive functions', 'all')

list_all_pfams1 <- list(CET_photosynthesis, other_photosynthesis,light_harvest,others_nitrogen,nitrate_sensing ,
                       specific_nitrogen, ferredoxin, flavodoxin, general_transporters,iron_transport,  silicon_transporter, 
                       phosphate_transp, nitrate_transp, ammonium_transp, DMSP, cell_div, hydrolases, glycolysis, iron_storage, 
                       copper_transport,molybdate_transport, zinc_transport )
names_pfams_list1 <- c('CET_photosynthesis', 'other_photosynthesis',
                      'light_harvest','others_nitrogen','n_sensing', 
                      'n_specific', 'ferredoxin', 'flavodoxin', 'general_transporters',
                      'iron_transporter' ,
                      'silicon_transporter', 'phosphate_transporter', 
                      'nitrate_transporter', 'ammonium_transporter', 'DMSP', 'cell division', 'hydrolase', 'glycolysis',
                      'iron storage', 'copper_transport', 'molybdate_transport', 'zinc_transport')


dat <- dat[dat[,111]>4,]
unigenes <- strsplit(unigenes, ' ')[[1]]
dat <- dat[dat[,1] %in% unigenes, ]


taxo='MGT-v2'
data_uni_tax <- readRDS(paste('../data/data_uni_',taxo,'.rds', sep=''))

tax_unis=readRDS('../data/taxID_uni_all.rds')
colnames(tax_unis)<-c('geneID', 'taxName')
tax_unis<- tax_unis[!duplicated(tax_unis$geneID),]
tax_unis <- tax_unis[tax_unis$geneID %in% unigenes,]


colnames(data_uni_tax)<-c('taxName', 'geneID')
data_uni_tax<- data_uni_tax[!duplicated(data_uni_tax$geneID),]
data_uni_tax <- data_uni_tax[data_uni_tax$geneID %in% unigenes,]

c_no_pfam <- 'no-Pfam'
col_not <- c_no_pfam
metabo <- rep(col_not, length(dat[,1]))
metabo_bis= rep(col_not, length(dat[,1]))
taxo_vec = data_uni_tax$taxName[match(dat[,1], data_uni_tax$geneID)]
taxo_vec[is.na(taxo_vec)]='no-taxo'
h=1
for (pfam_types in list_all_pfams){
  met=names_pfams_list[[h]]
  for (pf in pfam_types){
    metabo[grep(pf, dat[,49])] <-met
  }
  h=h+1
}

h=1
for (pfam_types in list_all_pfams1){
  met=names_pfams_list1[[h]]
  for (pf in pfam_types){
    metabo_bis[grep(pf, dat[,49])] <-met
  }
  h=h+1
}


if (taxo=='groups3'){
  taxos=c('Bacillariophyta', 'Mamiellales', 'Phaeocystales', 'Pelagophyceae')
} else{
  taxos=c("195_Cymatosiraceae", '76_root', "79_Bathycoccaceae", 
          "3_Phaeocystales","4_Haptista", "42_Pelagomonadales")
}

dat=rep(list(NULL), length(taxos))
names(dat)=taxos
for (tx in taxos){
  
  dt=tax_unis$taxName[match( data_uni_tax$geneID[data_uni_tax$taxName==tx], tax_unis$geneID )]
  unique_taxos=unique(dt)
  for (t in unique_taxos){
    co=sum(dt==t, na.rm=T)
    dat[[tx]]=rbind(dat[[tx]], c(t, co, co*100/length(data_uni_tax$geneID[data_uni_tax$taxName==tx])))
  }
  dat[[tx]]=dat[[tx]][order(as.numeric(dat[[tx]][,2]), decreasing = T),]
}


if (taxo=='groups3'){
  sep_vals=c(-0.2, 0, 0.2)
  cex_taxo=1.75
} else{
  sep_vals=c(-0.25, -0.125, 0, 0.125, 0.25)
  cex_taxo=1.25
}

vals_to_plot=c('mean1', 'atlanticity', 'delta1')

col_table=readRDS(paste('../data/color_table_',taxo,'.rds', sep=''))



pdf(paste('violin_plots_metabolisms_taxons_medians_',taxo,'.pdf', sep=''), height=6)
for (func in names_pfams_list){
  for (v in vals_to_plot){
    to_plot=rep(list(NULL), length(taxos))
    names(to_plot)=taxos
    for (tx in taxos){
      if (func !='all'){
        if (tx=="79_Bathycoccaceae"){
          dt=data_expression[[v]][metabo==func & taxo_vec==tx | taxo_vec=="439_Bathycoccaceae"]
        } else{
          dt=data_expression[[v]][metabo==func & taxo_vec==tx ]
        }
        
      } else{
        if (tx=="79_Bathycoccaceae"){
          dt=data_expression[[v]][taxo_vec==tx | taxo_vec=="439_Bathycoccaceae"]
        } else{
          dt=data_expression[[v]][taxo_vec==tx]
        }
        
      }
      
      dt=dt[!is.na(dt)]
      if (length(dt)==0){
        to_plot[[tx]]=0
      } else{
        to_plot[[tx]]=dt
      }
      
    }
    vec_to_plot=unlist(to_plot)
    g_to_plot=NULL
    lowercase_alphabet <- letters
    for (tx in taxos){
      if (length(to_plot[[tx]])>0){
        g_to_plot=c(g_to_plot, rep(tx, length(to_plot[[tx]])))
      }
      
    }
    
    t_test=pairwise.wilcox.test(x = vec_to_plot,g=g_to_plot, p.adjust.method = 'bonferroni')
    
    if (v=='atlanticity'){
      y_max=1.5
    } else if (v=='mean1'){
      y_max=0.88*max(data_expression[[v]], na.rm=T)
      if (func=='all'){
        y_max=max(data_expression[[v]], na.rm=T)
      }
    } else if (v=='delta1'){

      y_max=max(data_expression[[v]], na.rm=T)

    }

    
    vioplot::vioplot(to_plot, las=2, ylim=c(-y_max, y_max), main=paste(func, v, sep=' '), cex.axis=2, axes=F)

    stripchart(to_plot,cex=0.1, vertical = TRUE, method = "jitter",
               pch = 19, add = TRUE, col = rep('black', length(taxos)))
    
    
    
    
    lowercase_alphabet <- letters
    sigs_x=NULL
    sigs_text=NULL
    c=1
    cs=rep(1,length(taxos))
    
    if (length(t_test$p.value)!=0){
      lowercase_alphabet <- letters
      sigs_x=NULL
      sigs_text=NULL
      c=1
      cs=rep(1,length(taxos))
      
      for (i in 1:dim(t_test$p.value)[1]){
        for (j in 1:dim(t_test$p.value)[1]){
          if (j<=i){
            i_t=match(rownames(t_test$p.value)[i],taxos )
            j_t=match(colnames(t_test$p.value)[j],taxos )
            ci=cs[i_t]
            cj=cs[j_t]
            
            if (!is.na(t_test$p.value[i,j]) & t_test$p.value[i,j]<0.01){
              
              sigs_x=c(sigs_x, i_t+sep_vals[ci])
              sigs_x=c(sigs_x, j_t+sep_vals[cj])
              sigs_text=c(sigs_text, letters[c])
              sigs_text=c(sigs_text, letters[c])
              c=c+1
              cs[i_t]=cs[i_t]+1
              cs[j_t]=cs[j_t]+1
            }
          }
        }
        
      }
      if (length(sigs_x)>0){
        text(x=sigs_x, y=1.12*max(vec_to_plot), sigs_text, cex=cex_taxo)
      }
    }
    
    abline(h=0, lwd=2)
    j=1
    for (tx in taxos){
      segments(x0 = j-0.5, x1 = j+0.5, y0 = median(data_expression[[v]][taxo_vec==tx], na.rm=T),y1= median(data_expression[[v]][taxo_vec==tx], na.rm=T), lwd=4, col=col_table$col[col_table$taxon==tx])
      j=j+1
    }
    
    p_vals=NULL
    signs_m=NULL
    j=1
    for (dt in to_plot){
      tx=taxos[j]
      med_val=median(data_expression[[v]][taxo_vec==tx], na.rm=T)
      if (length(dt)>2){
        p_val=wilcox.test(dt, mu = med_val)
        p_vals=c(p_vals,p_val$p.value)
      } else{
        p_vals=c(p_vals,NA)
      }
      signs_m=c(signs_m, as.integer(median(dt)>med_val)+1)
      j=j+1
    }
    p_vals_adj=p.adjust(p_vals, method = 'bonferroni')
    threshs=c(0.001)
    signs_1=c('-', '+')
    signs_2=c('-','+')
    signs_3=c('-', '+')

    co=1
    done=NULL

    colos=NULL
    for (tx in taxos){
      colos=c(colos, col_table$col[col_table$taxon==tx])
    }
    
    for (t in threshs){
      p_sigs=p_vals_adj<t
      signifs=which(p_sigs==T)
      signs=signs_m[which(p_sigs==T)]
      cols=colos[which(p_sigs==T)]
      for (d in done){
        cols=cols[signifs!=d]
        signs=signs[signifs!=d]
        signifs=signifs[signifs!=d]
        #print(signs)
      }
      #print(cols)
      #print(signifs)
      if (t==0.05 & length(signifs)>0){
        text(signifs, y=0.95*y_max, signs_3[signs] , cex=3, col=cols)
      } else if (t==0.01 & length(signifs)>0){
        text(signifs, y=0.95*y_max, signs_2[signs] , cex=3, col=cols)
      } else if (t==0.001 & length(signifs)>0){
        text(signifs, y=0.95*y_max, signs_1[signs] , cex=3, col=cols)
      }
      
      for (sig in signifs){
        done=c(done, sig)
      }
      co=co+1
    }
    
    
    
  }
  # plot phate
  to_plot1=rep(list(NULL), length(taxos))
  to_plot2=rep(list(NULL), length(taxos))
  names(to_plot)=taxos
  for (tx in taxos){
    dt1=data_expression[['phate-1']][metabo==func & taxo_vec==tx]
    dt1=dt1[!is.na(dt1)]
    dt2=data_expression[['phate-2']][metabo==func & taxo_vec==tx]
    dt2=dt2[!is.na(dt2)]
    if (length(dt)==0){
      to_plot1[[tx]]=0
      to_plot2[[tx]]=0
    } else{
      to_plot1[[tx]]=dt1
      to_plot2[[tx]]=dt2
    }
    
  }
  
  for (i in 1:length(taxos)){
    tx=taxos[i]
    if (i==1){
      plot(to_plot1[[tx]], to_plot2[[tx]], xlim=c(min(data_expression[['phate-1']]), max(data_expression[['phate-1']])), ylim=c(min(data_expression[['phate-2']]), max(data_expression[['phate-2']])), col=col_table$col[col_table$taxon==tx],
           xlab='PHATE 1', ylab='PHATE 2', pch=19)
    } else{
      points(to_plot1[[tx]], to_plot2[[tx]], col=col_table$col[col_table$taxon==tx], pch=19)
    }
    
  }
}
dev.off()

pdf(paste('violin_plots_metabolisms_taxon_medians_',taxo,'_bis.pdf', sep=''), height = 6)
for (func in names_pfams_list1){
  for (v in vals_to_plot){
    to_plot=rep(list(NULL), length(taxos))
    names(to_plot)=taxos
    for (tx in taxos){
      if (func !='all'){
        if (tx=="79_Bathycoccaceae"){
          dt=data_expression[[v]][metabo_bis==func & taxo_vec==tx | taxo_vec=="439_Bathycoccaceae"]
        } else{
          dt=data_expression[[v]][metabo_bis==func & taxo_vec==tx ]
        }
        
      } else{
        if (tx=="79_Bathycoccaceae"){
          dt=data_expression[[v]][taxo_vec==tx | taxo_vec=="439_Bathycoccaceae"]
        } else{
          dt=data_expression[[v]][taxo_vec==tx]
        }
        
      }
      dt=dt[!is.na(dt)]
      if (length(dt)==0){
        to_plot[[tx]]=0
      } else{
        to_plot[[tx]]=dt
      }
      
    }
    vec_to_plot=unlist(to_plot)
    g_to_plot=NULL
    lowercase_alphabet <- letters
    for (tx in taxos){
      if (length(to_plot[[tx]])>0){
        g_to_plot=c(g_to_plot, rep(tx, length(to_plot[[tx]])))
      }
      
    }
    
    t_test=pairwise.wilcox.test(x = vec_to_plot,g=g_to_plot, p.adjust.method = 'bonferroni')
    
    #y_max=1.3*max(vec_to_plot, na.rm = T)
    if (v=='atlanticity'){
      y_max=1.5
    } else if (v=='mean1'){
      y_max=0.88*max(data_expression[[v]], na.rm=T)
    } else if (v=='delta1'){
      y_max=0.8*max(data_expression[[v]], na.rm=T)
    }
    #y_max=1.15*max(
    
    vioplot::vioplot(to_plot, las=2,  main=paste(func, v, sep=' '), cex.axis=2, ylim=c(-y_max, y_max), axes=F)
    stripchart(to_plot,cex=0.1, vertical = TRUE, method = "jitter",
               pch = 19, add = TRUE, col = rep('black', length(taxos)))
    
    
    
    if (length(t_test$p.value)!=0){
      lowercase_alphabet <- letters
      sigs_x=NULL
      sigs_text=NULL
      c=1
      cs=rep(1,length(taxos))
      
      for (i in 1:dim(t_test$p.value)[1]){
        for (j in 1:dim(t_test$p.value)[1]){
          if (j<=i){
            i_t=match(rownames(t_test$p.value)[i],taxos )
            j_t=match(colnames(t_test$p.value)[j],taxos )
            ci=cs[i_t]
            cj=cs[j_t]
            
            if (!is.na(t_test$p.value[i,j]) & t_test$p.value[i,j]<0.01){
              sigs_x=c(sigs_x, i_t+sep_vals[ci])
              sigs_x=c(sigs_x, j_t+sep_vals[cj])
              sigs_text=c(sigs_text, letters[c])
              sigs_text=c(sigs_text, letters[c])
              c=c+1
              cs[i_t]=cs[i_t]+1
              cs[j_t]=cs[j_t]+1
            }
          }
        }
        
      }
      if (length(sigs_x)>0){
        text(x=sigs_x, y=1.12*max(vec_to_plot), sigs_text, cex=cex_taxo)
      }
      
      
      j=1
      for (tx in taxos){
        segments(x0 = j-0.5, x1 = j+0.5, y0 = median(data_expression[[v]][taxo_vec==tx], na.rm=T),y1= median(data_expression[[v]][taxo_vec==tx], na.rm=T), lwd=4, col=col_table$col[col_table$taxon==tx])
        j=j+1
      }
      
      
      p_vals=NULL
      signs_m=NULL
      j=1
      for (dt in to_plot){
        tx=taxos[j]
        med_val=median(data_expression[[v]][taxo_vec==tx], na.rm=T)
        if (length(dt)>2){
          p_val=wilcox.test(dt, mu = med_val)
          p_vals=c(p_vals,p_val$p.value)
        } else{
          p_vals=c(p_vals,NA)
        }
        signs_m=c(signs_m, as.integer(median(dt)>med_val)+1)
        j=j+1
      }
      
      p_vals_adj=p.adjust(p_vals, method = 'bonferroni')
      threshs=c(0.001)
      signs_1=c('-', '+')
      signs_2=c('-','+')
      signs_3=c('-', '+')
      co=1
      done=NULL
      colos=NULL
      for (tx in taxos){
        colos=c(colos, col_table$col[col_table$taxon==tx])
      }
      
      for (t in threshs){
        p_sigs=p_vals_adj<t
        signifs=which(p_sigs==T)
        signs=signs_m[which(p_sigs==T)]
        cols=colos[which(p_sigs==T)]
        for (d in done){
          cols=cols[signifs!=d]
          signs=signs[signifs!=d]
          signifs=signifs[signifs!=d]
          #print(signs)
        }
        if (t==0.05 & length(signifs)>0){
          text(signifs, y=0.95*y_max, signs_3[signs] , cex=3, col=cols)
        } else if (t==0.01 & length(signifs)>0){
          text(signifs, y=0.95*y_max, signs_2[signs] , cex=3, col=cols)
        } else if (t==0.001 & length(signifs)>0){
          text(signifs, y=0.95*y_max, signs_1[signs] , cex=3, col=cols)
        }
         
        for (sig in signifs){
          done=c(done, sig)
        }
        co=co+1
      }
    }
  }
  
  # plot phate
  to_plot1=rep(list(NULL), length(taxos))
  to_plot2=rep(list(NULL), length(taxos))
  names(to_plot)=taxos
  for (tx in taxos){
    dt1=data_expression[['phate-1']][metabo_bis==func & taxo_vec==tx]
    dt1=dt1[!is.na(dt1)]
    dt2=data_expression[['phate-2']][metabo_bis==func & taxo_vec==tx]
    dt2=dt2[!is.na(dt2)]
    if (length(dt)==0){
      to_plot1[[tx]]=0
      to_plot2[[tx]]=0
    } else{
      to_plot1[[tx]]=dt1
      to_plot2[[tx]]=dt2
    }
    
  }
  
  for (i in 1:length(taxos)){
    tx=taxos[i]
    if (i==1){
      plot(to_plot1[[tx]], to_plot2[[tx]], xlim=c(min(data_expression[['phate-1']]), max(data_expression[['phate-1']])), ylim=c(min(data_expression[['phate-2']]), max(data_expression[['phate-2']])), col=col_table$col[col_table$taxon==tx],
           xlab='PHATE 1', ylab='PHATE 2', pch=19)
    } else{
      points(to_plot1[[tx]], to_plot2[[tx]], col=col_table$col[col_table$taxon==tx], pch=19)
    }
    
  }
  
}
dev.off()


pdf(paste('violin_plots_metabolisms_taxons_and_community_medians_',taxo,'.pdf', sep=''), height=6)
for (func in names_pfams_list){
  for (v in vals_to_plot){
    to_plot=rep(list(NULL), length(taxos))
    names(to_plot)=taxos
    for (tx in taxos){
      if (func !='all'){
        if (tx=="79_Bathycoccaceae"){
          dt=data_expression[[v]][metabo==func & taxo_vec==tx | taxo_vec=="439_Bathycoccaceae"]
        } else{
          dt=data_expression[[v]][metabo==func & taxo_vec==tx ]
        }
        
      } else{
        if (tx=="79_Bathycoccaceae"){
          dt=data_expression[[v]][taxo_vec==tx | taxo_vec=="439_Bathycoccaceae"]
        } else{
          dt=data_expression[[v]][taxo_vec==tx]
        }
        
      }
      
      dt=dt[!is.na(dt)]
      if (length(dt)==0){
        to_plot[[tx]]=0
      } else{
        to_plot[[tx]]=dt
      }
      
    }
    vec_to_plot=unlist(to_plot)
    g_to_plot=NULL
    lowercase_alphabet <- letters
    for (tx in taxos){
      if (length(to_plot[[tx]])>0){
        g_to_plot=c(g_to_plot, rep(tx, length(to_plot[[tx]])))
      }
      
    }
    
    t_test=pairwise.wilcox.test(x = vec_to_plot,g=g_to_plot, p.adjust.method = 'bonferroni')
    
    if (v=='atlanticity'){
      y_max=1.5
    } else if (v=='mean1'){
      y_max=max(data_expression[[v]], na.rm=T)
      if (func=='all'){
        y_max=max(data_expression[[v]], na.rm=T)
      }
    } else if (v=='delta1'){
      y_max=max(data_expression[[v]], na.rm=T)
      #}
    }
    
    vioplot::vioplot(to_plot, las=2, ylim=c(-y_max, y_max), main=paste(func, v, sep=' '), cex.axis=2, axes=F)
    stripchart(to_plot,cex=0.1, vertical = TRUE, method = "jitter",
                 pch = 19, add = TRUE, col = rep('black', length(taxos)))
    #}
    
    
    
    
    
    lowercase_alphabet <- letters
    sigs_x=NULL
    sigs_text=NULL
    c=1
    cs=rep(1,length(taxos))
    
    if (length(t_test$p.value)!=0){
      lowercase_alphabet <- letters
      sigs_x=NULL
      sigs_text=NULL
      c=1
      cs=rep(1,length(taxos))
      
      for (i in 1:dim(t_test$p.value)[1]){
        for (j in 1:dim(t_test$p.value)[1]){
          if (j<=i){
            i_t=match(rownames(t_test$p.value)[i],taxos )
            j_t=match(colnames(t_test$p.value)[j],taxos )
            ci=cs[i_t]
            cj=cs[j_t]
            
            if (!is.na(t_test$p.value[i,j]) & t_test$p.value[i,j]<0.01){
              
              sigs_x=c(sigs_x, i_t+sep_vals[ci])
              sigs_x=c(sigs_x, j_t+sep_vals[cj])
              sigs_text=c(sigs_text, letters[c])
              sigs_text=c(sigs_text, letters[c])
              c=c+1
              cs[i_t]=cs[i_t]+1
              cs[j_t]=cs[j_t]+1
            }
          }
        }
        
      }
      if (length(sigs_x)>0){
        text(x=sigs_x, y=1.07*min(vec_to_plot), sigs_text, cex=cex_taxo)
      }
    }
    
    abline(h=0, lwd=2)
    j=1
    for (tx in taxos){
      segments(x0 = j-0.5, x1 = j+0.5, y0 = median(data_expression[[v]][taxo_vec==tx], na.rm=T),y1= median(data_expression[[v]][taxo_vec==tx], na.rm=T), lwd=4, col=col_table$col[col_table$taxon==tx])
      j=j+1
    }
    
    p_vals=NULL
    signs_m=NULL
    j=1
    for (dt in to_plot){
      tx=taxos[j]
      med_val=median(data_expression[[v]][taxo_vec==tx], na.rm=T)
      if (length(dt)>2){
        p_val=wilcox.test(dt, mu = med_val)
        p_vals=c(p_vals,p_val$p.value)
      } else{
        p_vals=c(p_vals,NA)
      }
      signs_m=c(signs_m, as.integer(median(dt)>med_val)+1)
      j=j+1
    }
    p_vals_adj=p.adjust(p_vals, method = 'bonferroni')
    threshs=c(0.001)
    signs_1=c('-', '+')
    signs_2=c('-','+')
    signs_3=c('-', '+')
    co=1
    done=NULL
    colos=NULL
    for (tx in taxos){
      colos=c(colos, col_table$col[col_table$taxon==tx])
    }
    
    for (t in threshs){
      p_sigs=p_vals_adj<t
      signifs=which(p_sigs==T)
      signs=signs_m[which(p_sigs==T)]
      cols=colos[which(p_sigs==T)]
      for (d in done){
        cols=cols[signifs!=d]
        signs=signs[signifs!=d]
        signifs=signifs[signifs!=d]
      }
      if (func!='all'){
        if (t==0.05 & length(signifs)>0){
          text(signifs, y=0.85*y_max, signs_3[signs] , cex=3, col=cols)
        } else if (t==0.01 & length(signifs)>0){
          text(signifs, y=0.85*y_max, signs_2[signs] , cex=3, col=cols)
        } else if (t==0.001 & length(signifs)>0){
          text(signifs, y=0.85*y_max, signs_1[signs] , cex=3, col=cols)
        }
      }
      
      for (sig in signifs){
        done=c(done, sig)
      }
      co=co+1
    }
    
    abline(h=0, lwd=2)
    p_vals=NULL
    signs_m=NULL
    for (dt in to_plot){
      if (length(dt)>2){
        p_val=wilcox.test(dt, mu = 0)
        p_vals=c(p_vals,p_val$p.value)
      } else{
        p_vals=c(p_vals,NA)
      }
      signs_m=c(signs_m, as.integer(mean(dt)>0)+1)
    }
    p_vals_adj=p.adjust(p_vals, method = 'bonferroni')
    threshs=c(0.001)#, 0.01, 0.05)
    signs=c('*')# '**', '*')
    co=1
    done=NULL
    if (v=='mean1'){
      colos=c('blue', 'red')
    } else {
      colos=c('darkorange', 'darkblue')
    } 
    for (t in threshs){
      p_sigs=p_vals_adj<t
      signifs=which(p_sigs==T)
      cols=signs_m[which(p_sigs==T)]
      cols=colos[cols]
      for (d in done){
        cols=cols[signifs!=d]
        signifs=signifs[signifs!=d]
      }
      text(signifs, y=y_max, signs[co] , cex=4, col=cols)
      for (sig in signifs){
        done=c(done, sig)
      }
      co=co+1
    }
    
    
    
  }
  # plot phate
  to_plot1=rep(list(NULL), length(taxos))
  to_plot2=rep(list(NULL), length(taxos))
  names(to_plot)=taxos
  for (tx in taxos){
    dt1=data_expression[['phate-1']][metabo==func & taxo_vec==tx]
    dt1=dt1[!is.na(dt1)]
    dt2=data_expression[['phate-2']][metabo==func & taxo_vec==tx]
    dt2=dt2[!is.na(dt2)]
    if (length(dt)==0){
      to_plot1[[tx]]=0
      to_plot2[[tx]]=0
    } else{
      to_plot1[[tx]]=dt1
      to_plot2[[tx]]=dt2
    }
    
  }
  
  for (i in 1:length(taxos)){
    tx=taxos[i]
    if (i==1){
      plot(to_plot1[[tx]], to_plot2[[tx]], xlim=c(min(data_expression[['phate-1']]), max(data_expression[['phate-1']])), ylim=c(min(data_expression[['phate-2']]), max(data_expression[['phate-2']])), col=col_table$col[col_table$taxon==tx],
           xlab='PHATE 1', ylab='PHATE 2', pch=19)
    } else{
      points(to_plot1[[tx]], to_plot2[[tx]], col=col_table$col[col_table$taxon==tx], pch=19)
    }
    
  }
}
dev.off()

pdf(paste('violin_plots_metabolisms_taxon_and_community_medians_bis_',taxo,'.pdf'), height = 6)
for (func in names_pfams_list1){
  for (v in vals_to_plot){
    to_plot=rep(list(NULL), length(taxos))
    names(to_plot)=taxos
    for (tx in taxos){
      if (func !='all'){
        if (tx=="79_Bathycoccaceae"){
          dt=data_expression[[v]][metabo_bis==func & taxo_vec==tx | taxo_vec=="439_Bathycoccaceae"]
        } else{
          dt=data_expression[[v]][metabo_bis==func & taxo_vec==tx ]
        }
        
      } else{
        if (tx=="79_Bathycoccaceae"){
          dt=data_expression[[v]][taxo_vec==tx | taxo_vec=="439_Bathycoccaceae"]
        } else{
          dt=data_expression[[v]][taxo_vec==tx]
        }
        
      }
      dt=dt[!is.na(dt)]
      if (length(dt)==0){
        to_plot[[tx]]=0
      } else{
        to_plot[[tx]]=dt
      }
      
    }
    vec_to_plot=unlist(to_plot)
    g_to_plot=NULL
    lowercase_alphabet <- letters
    for (tx in taxos){
      if (length(to_plot[[tx]])>0){
        g_to_plot=c(g_to_plot, rep(tx, length(to_plot[[tx]])))
      }
      
    }
    
    t_test=pairwise.wilcox.test(x = vec_to_plot,g=g_to_plot, p.adjust.method = 'bonferroni')
    if (v=='atlanticity'){
      y_max=1.5
    } else if (v=='mean1'){
      y_max=max(data_expression[[v]], na.rm=T)
    } else if (v=='delta1'){
      y_max=0.7*max(data_expression[[v]], na.rm=T)
    }
    #y_max=1.15*max(
    
    vioplot::vioplot(to_plot, las=2,  main=paste(func, v, sep=' '), cex.axis=2, ylim=c(-y_max, y_max), axes=F)
    stripchart(to_plot,cex=0.1, vertical = TRUE, method = "jitter",
               pch = 19, add = TRUE, col = rep('black', length(taxos)))
    
    
    
    if (length(t_test$p.value)!=0){
      lowercase_alphabet <- letters
      sigs_x=NULL
      sigs_text=NULL
      c=1
      cs=rep(1,length(taxos))
      
      for (i in 1:dim(t_test$p.value)[1]){
        for (j in 1:dim(t_test$p.value)[1]){
          if (j<=i){
            i_t=match(rownames(t_test$p.value)[i],taxos )
            j_t=match(colnames(t_test$p.value)[j],taxos )
            ci=cs[i_t]
            cj=cs[j_t]
            
            if (!is.na(t_test$p.value[i,j]) & t_test$p.value[i,j]<0.01){
              sigs_x=c(sigs_x, i_t+sep_vals[ci])
              sigs_x=c(sigs_x, j_t+sep_vals[cj])
              sigs_text=c(sigs_text, letters[c])
              sigs_text=c(sigs_text, letters[c])
              c=c+1
              cs[i_t]=cs[i_t]+1
              cs[j_t]=cs[j_t]+1
            }
          }
        }
        
      }
      if (length(sigs_x)>0){
        text(x=sigs_x, y=1.12*min(vec_to_plot), sigs_text, cex=cex_taxo)
      }
      
      
      j=1
      for (tx in taxos){
        segments(x0 = j-0.5, x1 = j+0.5, y0 = median(data_expression[[v]][taxo_vec==tx], na.rm=T),y1= median(data_expression[[v]][taxo_vec==tx], na.rm=T), lwd=4, col=col_table$col[col_table$taxon==tx])
        j=j+1
      }
      
      
      p_vals=NULL
      signs_m=NULL
      j=1
      for (dt in to_plot){
        tx=taxos[j]
        med_val=median(data_expression[[v]][taxo_vec==tx], na.rm=T)
        if (length(dt)>2){
          p_val=wilcox.test(dt, mu = med_val)
          p_vals=c(p_vals,p_val$p.value)
        } else{
          p_vals=c(p_vals,NA)
        }
        signs_m=c(signs_m, as.integer(median(dt)>med_val)+1)
        j=j+1
      }
      
      p_vals_adj=p.adjust(p_vals, method = 'bonferroni')
      threshs=c(0.001)
      signs_1=c('-', '+')
      signs_2=c('-','+')
      signs_3=c('-', '+')
      co=1
      done=NULL
      colos=NULL
      for (tx in taxos){
        colos=c(colos, col_table$col[col_table$taxon==tx])
      }
      
      for (t in threshs){
        p_sigs=p_vals_adj<t
        signifs=which(p_sigs==T)
        #print(done)
        signs=signs_m[which(p_sigs==T)]
        #print(signs)
        cols=colos[which(p_sigs==T)]
        for (d in done){
          cols=cols[signifs!=d]
          signs=signs[signifs!=d]
          signifs=signifs[signifs!=d]
          #print(signs)
        }
        if (t==0.05 & length(signifs)>0){
          text(signifs, y=0.85*y_max, signs_3[signs] , cex=3, col=cols)
        } else if (t==0.01 & length(signifs)>0){
          text(signifs, y=0.85*y_max, signs_2[signs] , cex=3, col=cols)
        } else if (t==0.001 & length(signifs)>0){
          text(signifs, y=0.85*y_max, signs_1[signs] , cex=3, col=cols)
        }
        
        for (sig in signifs){
          done=c(done, sig)
        }
        co=co+1
      }
      
      abline(h=0, lwd=2)
      p_vals=NULL
      signs_m=NULL
      for (dt in to_plot){
        if (length(dt)>2){
          p_val=wilcox.test(dt, mu = 0)
          p_vals=c(p_vals,p_val$p.value)
        } else{
          p_vals=c(p_vals,NA)
        }
        signs_m=c(signs_m, as.integer(mean(dt)>0)+1)
      }
      p_vals_adj=p.adjust(p_vals, method = 'bonferroni')
      threshs=c(0.001)#, 0.01, 0.05)
      signs=c('*')#, '**', '*')
      co=1
      done=NULL
      if (v=='mean1'){
        colos=c('blue', 'red')
      } else {
        colos=c('darkorange', 'darkblue')
      } 
      for (t in threshs){
        p_sigs=p_vals_adj<t
        signifs=which(p_sigs==T)
        cols=signs_m[which(p_sigs==T)]
        cols=colos[cols]
        for (d in done){
          cols=cols[signifs!=d]
          signifs=signifs[signifs!=d]
        }
        #print(cols)
        #print(signifs)
        text(signifs, y=y_max, signs[co] , cex=4, col=cols)
        for (sig in signifs){
          done=c(done, sig)
        }
        co=co+1
      }
    }
  }
  
  # plot phate
  to_plot1=rep(list(NULL), length(taxos))
  to_plot2=rep(list(NULL), length(taxos))
  names(to_plot)=taxos
  for (tx in taxos){
    dt1=data_expression[['phate-1']][metabo_bis==func & taxo_vec==tx]
    dt1=dt1[!is.na(dt1)]
    dt2=data_expression[['phate-2']][metabo_bis==func & taxo_vec==tx]
    dt2=dt2[!is.na(dt2)]
    if (length(dt)==0){
      to_plot1[[tx]]=0
      to_plot2[[tx]]=0
    } else{
      to_plot1[[tx]]=dt1
      to_plot2[[tx]]=dt2
    }
    
  }
  
  for (i in 1:length(taxos)){
    tx=taxos[i]
    if (i==1){
      plot(to_plot1[[tx]], to_plot2[[tx]], xlim=c(min(data_expression[['phate-1']]), max(data_expression[['phate-1']])), ylim=c(min(data_expression[['phate-2']]), max(data_expression[['phate-2']])), col=col_table$col[col_table$taxon==tx],
           xlab='PHATE 1', ylab='PHATE 2', pch=19)
    } else{
      points(to_plot1[[tx]], to_plot2[[tx]], col=col_table$col[col_table$taxon==tx], pch=19)
    }
    
  }
  
}
dev.off()


y_max=max(abs(raw_expr[taxo_vec %in% taxos,]))
stats=c('144', '145', '147', '148', '150', '151', '152', '155', '158', '163', '168', '173', '175', '178', '180', '188', '189',
        '193', '194', '196')
pdf(paste('expression_vioplot_stations_',taxo,'.pdf', sep=''), height = 6)
for (func in names_pfams_list){
  for (tx in taxos){
    if (func!='all'){
      raw_expr_phot=raw_expr[taxo_vec==tx & metabo==func,]
    } else{
      raw_expr_phot=raw_expr[taxo_vec==tx,]
    }
    
    min_x=min(raw_expr)/2
    raw_expr_phot[raw_expr_phot==0]=NA
    dat_to_plot=rep(list(NULL), dim(raw_expr_phot)[2])
    for (i in 1:dim(raw_expr_phot)[2]){
      dat_to_plot[[i]]=raw_expr_phot[[i]]
      if (sum(raw_expr_phot[[i]], na.rm=T)==0){
        dat_to_plot[[i]]=0
      }
    }
    mean_arc=median(unlist(raw_expr_phot[,9:20]), na.rm=T)
    mean_atl=median(unlist(raw_expr_phot[,1:8]), na.rm=T)
    vioplot(dat_to_plot, las=2,  main=paste(func, tx, sep=' '), cex.axis=1, ylim=c(-y_max, y_max), names=stats,
            col=c(rep('darkorange', 8), rep('navyblue', 12)))
    #box(lwd = 3)
    #if (func!='all'){
    stripchart(dat_to_plot,cex=0.05, vertical = TRUE, method = "jitter",
                 pch = 19, add = TRUE, col = rep('black', length(taxos)))
    #}
    
    abline(h=0, cex=2)
    segments(x0 = 0, x1=8.5, y0 = mean_atl,y1=mean_atl, lwd=3, col='darkorange')
    segments(x0=8.5, x1=21,y0=mean_arc, y1=mean_arc, lwd=3, col='navyblue')
    
    
  }
  
}
dev.off()


pdf(paste('expression_vioplot_stations_bis_',taxo,'.pdf', sep=''), height = 6)
for (func in names_pfams_list1){
  for (tx in taxos){
    if (func!='all'){
      raw_expr_phot=raw_expr[taxo_vec==tx & metabo_bis==func,]
    } else{
      raw_expr_phot=raw_expr[taxo_vec==tx,]
    }
    
    min_x=min(raw_expr)/2
    raw_expr_phot[raw_expr_phot==0]=NA
    dat_to_plot=rep(list(NULL), dim(raw_expr_phot)[2])
    for (i in 1:dim(raw_expr_phot)[2]){
      dat_to_plot[[i]]=raw_expr_phot[[i]]
      if (sum(raw_expr_phot[[i]], na.rm=T)==0){
        dat_to_plot[[i]]=0
      }
    }
    mean_arc=median(unlist(raw_expr_phot[,9:20]), na.rm=T)
    mean_atl=median(unlist(raw_expr_phot[,1:8]), na.rm=T)
    vioplot(dat_to_plot, las=2,  main=paste(func, tx, sep=' '), cex.axis=1, ylim=c(-y_max, y_max), names=stats,
            col=c(rep('darkorange', 8), rep('navyblue', 12)))
    #box(lwd = 3)
    stripchart(dat_to_plot,cex=0.05, vertical = TRUE, method = "jitter",
               pch = 19, add = TRUE, col = rep('black', length(taxos)))
    abline(h=0, cex=2)
    segments(x0 = 0, x1=8.5, y0 = mean_atl,y1=mean_atl, lwd=3, col='darkorange')
    segments(x0=8.5, x1=21,y0=mean_arc, y1=mean_arc, lwd=3, col='navyblue')
    
    
  }
  
}
dev.off()



