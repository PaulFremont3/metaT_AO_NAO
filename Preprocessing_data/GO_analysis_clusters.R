#!/bin/env/usr/env Rscript
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/BioAdvection_II/MetaT_4")

GO_table <- readRDS('GO_table.rds')

fractions <- c('GGZZ')
n_clusts <- 6
taxo <- commandArgs(trailingOnly = T)[1]
bis <- commandArgs(trailingOnly = T)[2]

types <- c('T', 'G')

#tests=c('hommel', 'fdr')

if (taxo=='taxo_groups'){
  groups <-c(NULL,'Hexanauplia', 'Bacillariophyta', 'Coccosphaerales', 'Ciliophora', 'Tunicata',
             'Mamiellales', 'Dinophyceae', 'unknown', 'other')
} else if (taxo=='taxo_groups2'){
  groups <- c(NULL, 'Hexanauplia', 'Bacillariophyta','Pelagophyceae','Crysophyceae' ,'Coccosphaerales','Phaeocystales' ,
           'Ciliophora', 'Tunicata','Mamiellales', 'Dinophyceae', 'Bacteria', 'Streptophyta', 'Cnidaria','Rhizaria',
           'Insecta', 'Cnidaria', 'Cryptophyta' , 'Euglenozoa', 'Amoebozoa')
} else if (taxo=='taxo_groups3'){
  groups <- c(NULL,'Hexanauplia', 'Bacillariophyta','Pelagophyceae','Crysophyceae' ,'Coccosphaerales','Phaeocystales' ,
           'Ciliophora', 'Tunicata','Mamiellales', 'Dinophyceae', 'Bacteria', 'Streptophyta', 'Cnidaria','Rhizaria',
           'Insecta', 'Cryptophyta' , 'Euglenozoa', 'Amoebozoa','Craniata','unknown','other',
           'Fungi', 'Viruses','Archaea' ,'Eukaryota (unclassified)','root' ,paste('other ' , c('Haptophyta', 'Opisthokonta', 'Stramenopiles', 'Alveolata', 'Viridiplantae') , sep='' ) )
}else{
  groups <-NULL
}

for (fraction in fractions){
  for (type in types){

    GO_data <- readRDS(paste('GO_station_table_',fraction,'_',n_clusts,'_',
                                  type,'_',taxo,bis,'.rds', sep=''))
    
    
    if (fraction=='SSUU'){
      arctic_stations <- c("158SUR"  ,"173SUR", "175SUR", "178SUR","180SUR", "188SUR", "189SUR",
                           "193SUR", "194SUR", "196SUR")
      atlantic_stations <- c("143SUR", "144SUR", "145SUR", "146SUR", "147SUR", "148SUR", "149SUR", "150SUR" ,
                            "151SUR", "152SUR")
      outlier_stations <-c('155SUR','163SUR', '168SUR')
    } else if (fraction=='QQSS'){
      arctic_stations <- c( "158SUR"   ,"168SUR" ,"173SUR", "175SUR", "178SUR","180SUR", "188SUR", "189SUR", "191SUR",
                           "193SUR", "194SUR", "196SUR")
      atlantic_stations <- c("142SUR", "143SUR", "144SUR","145SUR","146SUR", "147SUR", "148SUR", "149SUR", "150SUR" ,
                             "151SUR" ,"152SUR" ,"155SUR")
      outlier_stations <-c('')
    } else if (fraction=='GGZZ'){
      arctic_stations <- c( "168SUR" ,"173SUR", "175SUR", "178SUR","180SUR", "188SUR", "189SUR",
                           "193SUR", "194SUR", "196SUR")
      atlantic_stations <- c("142SUR", "144SUR", "145SUR",  "147SUR", "148SUR" ,
                             "151SUR", "152SUR")
      outlier_stations <-c('155SUR',"158SUR",'168SUR', "150SUR")
    } else if (fraction=='KKQQ'){
      arctic_stations <- c("155SUR" ,"158SUR","168SUR" ,"173SUR", "175SUR", "178SUR","180SUR", "188SUR", "189SUR",
                           "193SUR", "194SUR")
      atlantic_stations <- c("142SUR","143SUR" ,"144SUR", "145SUR","146SUR" , "147SUR", "148SUR","149SUR" ,"150SUR" ,
                             "151SUR", "152SUR")
      outlier_stations <-c('196SUR')
    }
    arctic_stations <- c("158SUR","168SUR" ,"163SUR","173SUR", "175SUR", "178SUR","180SUR", "188SUR", "189SUR","191SUR",
                           "193SUR", "194SUR", "196SUR")
    atlantic_stations <- c("142SUR","143SUR" ,"144SUR", "145SUR","146SUR" , "147SUR", "148SUR","149SUR" ,"150SUR" ,
                             "151SUR", "152SUR", '155SUR')
 
    
    outliers_functions <- function(GO_data, gr){
      p_vals <- NULL
      GOs <- NULL
      loc <- NULL
      GOs_all <- NULL
      for (nam in names(GO_data)){
        t0 <- GO_data[[nam]]$all
        if (!is.null(gr)){
          t0 <- t0[rownames(t0)==gr,]
        }
        if (!is.null(t0)){
          d_t0 <- dim(t0)[1]
          t0 <- t0[,!(colnames(t0) %in% c('142SUR','201SUR','205SUR',
                                          '206SUR', '208SUR', '209SUR', '210SUR'))]
          if (d_t0==1){
            t1 <- t(t0)
            sum_metaT <- t1
          } else{
            sum_metaT <- apply(t0, 2, sum)
          }
          
          
          if (sum(sum_metaT==0)==0){
            cond_atl <- colnames(t0) %in% atlantic_stations
            cond_arc <- colnames(t0) %in% arctic_stations
            w_test <- t.test(log(sum_metaT[cond_atl]), log(sum_metaT[cond_arc]))$p.value
            if (mean(sum_metaT[cond_atl]) > mean(sum_metaT[cond_arc])){
              loca <- 'atlantic'
            } else{
              loca <- 'arctic'
            }
            p_vals <- append(p_vals, w_test)
            GOs <- append(GOs, nam)
            loc <- append(loc, loca)
            GOs_all <- append(GOs_all, nam)
          } else{
            GOs_all <- append(GOs_all, nam)
          }
        }
      }
      saveRDS(GOs_all, paste('baseline_GO_all_',fraction,bis,'.rds', sep=''))
      test_holm <- p.adjust(p_vals, method = 'hommel')
      saveRDS(test_holm, paste('baseline_GO_all_',fraction,bis,'_p_values.rds', sep='')) 
     
      arc <- GO_table[GO_table[,1] %in% GOs[test_holm<0.05 & loc=='arctic'] ,1]
      atl <- GO_table[GO_table[,1] %in% GOs[test_holm<0.05 & loc=='atlantic'] ,1]
      
      arc_gonames <- GO_table[GO_table[,1] %in% GOs[test_holm<0.05 & loc=='arctic'] ,2]
      atl_gonames <- GO_table[GO_table[,1] %in% GOs[test_holm<0.05 & loc=='atlantic'] ,2]
      
      if (is.null(gr)){
	      print(gr)
        saveRDS(arc, paste('GO_representative_',type,'_arctic_', fraction,'_' ,taxo,bis, '.rds', sep=''))
        saveRDS(atl, paste('GO_representative_',type,'_atlantic_', fraction, '_' ,taxo, bis,'.rds', sep=''))
      } else{
        saveRDS(arc, paste('GO_representative_',type,'_arctic_', fraction,'_' ,taxo,'_', gr,bis,'.rds', sep=''))
        saveRDS(atl, paste('GO_representative_',type,'_atlantic_', fraction, '_' ,taxo,'_', gr,bis, '.rds', sep=''))
      }

      test_fdr <- p.adjust(p_vals, method = 'fdr')
      #saveRDS(test_holm, paste('baseline_GO_all_',fraction,bis,'_p_values.rds', sep=''))

      arc <- GO_table[GO_table[,1] %in% GOs[test_fdr<0.05 & loc=='arctic'] ,1]
      atl <- GO_table[GO_table[,1] %in% GOs[test_fdr<0.05 & loc=='atlantic'] ,1]

      arc_gonames <- GO_table[GO_table[,1] %in% GOs[test_fdr<0.05 & loc=='arctic'] ,2]
      atl_gonames <- GO_table[GO_table[,1] %in% GOs[test_fdr<0.05 & loc=='atlantic'] ,2]

      if (is.null(gr)){
              print(gr)
        saveRDS(arc, paste('GO_representative_',type,'_arctic_', fraction,'_' ,taxo,bis, '_fdr.rds', sep=''))
        saveRDS(atl, paste('GO_representative_',type,'_atlantic_', fraction, '_' ,taxo, bis,'_fdr.rds', sep=''))
      }
 
    }
    outliers_functions(GO_data = GO_data,gr = NULL)
    for (gr in groups){
      outliers_functions(GO_data = GO_data,gr = gr)
      print(gr)
    }
  }
}




