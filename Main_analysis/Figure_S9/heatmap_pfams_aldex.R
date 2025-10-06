library('pheatmap')
library('ALDEx2')
library('gplots')
library(gridExtra)

frac <-commandArgs(trailingOnly = T)[1]
aldex_files<- list.files(pattern=paste('../data/mTG_',frac,'.rds', sep=''))
aldex_files <- aldex_files[!grepl('unigenes', aldex_files)]
aldex_files <- aldex_files[c(6, 5, 8, 7, 1, 4, 3,9)]
taxos <- c('All', 'Mamiellales', 'Phaeocystales', 'Pelagophyceae', 'Bacillariophyta', 'Dinophyceae', 'Ciliophora', 'unknown')

photosynthesis_pfams <- c('PF11623','PF02276' , 'PF01716','PF00504','PF00556' ,'PF00124', 'PF03967', 'PF00223', 'PF02605', 'PF07465', 'PF05479', 'PF00737', 'PF02532', 'PF02533', 'PF05151'
                          , 'PF02468', 'PF04725', 'PF01405', 'PF06514', 'PF07123', 'PF06596', 'PF06298' , 'PF00421', 'PF05969', 'PF13326', 'PF00796',
                          'PF02427', 'PF02507', 'PF03244', 'PF01701', 'PF01241', 'PF10657', 'PF11947')
names_photosynthesis_pfams <- c('NdhS','CytoC_RC', 'MSP', 'Chla_b-bind','LHC','Photo_RC', 'PRCH', 'PsaA_PsaB', 'PsaL', 'PsaM', 'PsaN', 'PsbH', 'PsbI', 'PsbK', 'PsbM',
                                'PsbN', 'PsbR', 'PsbT', 'PsbU', 'PsbW', 'PsbX', 'PsbY', 'PSII', 'PSII_Ycf12', 'PSII_Pbs27', 'PSI_8' ,
                                'PSI_PsaE', 'PSI_PsaF', 'PSI_PsaH', 'PSI_PsaJ', 'PSI_PSAK', 'RC-P840_PscD', 'DUF3464')
locations <- c('absent', 'absent', 'absent', 'absent', 'absent', 'chromosome', 'chloroplast', 'chromosome', 'chloroplast', 'chloroplast', 'chloroplast', 'chromosome',
               'chloroplast', 'chromosome', 'absent', 'absent', 'chromosome', 'chromosome', 'chromosome', 'chloroplast', 'chloroplast', 'absent', 'chloroplast',
               'chromosome', 'chromosome', 'chromosome', 'chloroplast', 'chloroplast', 'absent', 'absent')
colos <- c('darkgreen', 'blue', 'black')
nas <- c('chloroplast', 'chromosome', 'absent')

nitrogen_pfams <-c('PF00142','PF04891', 'PF03206', 'PF04319', 'PF06988', 'PF07732', 'PF13473', 'PF00115', 'PF13442',
                   'PF04879', 'PF02335', 'PF01077','PF03460', 'PF00384',
                   'PF01568', 'PF04324', 'PF00174', 'PF13435', 'PF07992', 'PF00355')

names_nitrogen_pfams <- c('Fer4_NifH', 'NifQ','NifW','NifZ' , 'NifT', 'Cu-oxidase_3 (nirK)', 'Cupredoxin_1 (nosZ)', 'COX1 (norB)', 'Cytochrome_CBB3 (nirS)',
                          'Molybdop_Fe4S4 (napA/nasA)', 'Cytochrom_C552 (nrfA)', 'NIR_SIR (nirA/NIT-6)', 'NIR_SIR_ferr (nirA/NIT-6)', 'Molybdopterin (napA/nasA)',
                          'Molydop_binding (napA/nasA)', 'Fer2_BFD (nasA/NIT-6)', 'Oxydored_molyb (NR)', 'Cytochrome_C554 (hao)', 'Pyr_redox_2 (NIT-6)',
                          'Rieske (NIT-6)')

carbon_fixation_pfams <- c('PF13452', 'PF00485', 'PF00016', 'PF00101', 'PF06240', 'PF07690' , 'PF00194', 'PF00126')
names_carbon_fixation_pfams <- c('MaoC_dehydrat_N (meh)', 'PRK', 'RuBisCO_large', 'RuBisCO_small','COXG (coxSML)', 'MFS_1 (mct)', 'Carb_anhydrase',
                                 'RuBisCO cbbR')

iron_functions <- c( 'PF00111', 'PF01799', 'PF13085', 'PF00258', 'PF02525')
names_iron_function <- c( 'Fer2', 'Fer2_3', 'Fer2_3', 'Flavodoxin_1' , 'Flavodoxin_2')

other_functions <- c('PF03842', 'PF16867', 'PF00127', 'PF16983','PF11124',  'PF01384', 'PF16974')
names_other_function <- c('Silic_transp', 'DMSP_lyase', 'Copper-bind', ' Molybdate transporter MFS_MOT1', 'Pho86','Phosphate transporter' ,
                          'Nitrate transporter (NAR2)')

temperature_pfams <- c('PF00360' ,'PF03952','PF01786','PF06415','PF05971', 'PF01180', 'PF05351', 'PF05035', 'PF01346', 'PF03104', 'PF11999',
                       'PF00249', 'PF03259', 'PF02209', 'PF00241', 'PF00626', 'PF01119','PF10551', 'PF01979',
                       'PF01964', 'PF02358','PF14249', 'PF02649',
                       'PF00313', 'PF05562', 'PF14169')
names_temperature_pfams <- c('PHY' ,'Enolase','AOX','PGAM','Methyltransf_10', 'DHO_dh', 'GMP_PDE_delta', 'DGOK', 'FKBP_N', 'DNA_pol_B_exo1', 'Ice_binding',
                             'Myb_DNA-binding (CCA1)', 'Robl_LC7', 'VHP', 'Cofilin_ADF', 'Gelsolin', 'DNA_mis_repair', 'MULE', 'Amidohydro_1 (allantoinase)',
                             'Thiamine biosynthesis','Trehalose_PPase', 'Tocopherol_cycl', 'GCHY-1 (tetrahydrofolate biosynthesis)',
                             'CSD', 'WCOR413', 'Cold-inducible protein YdjO')
cell_cycle_pfams <- c('PF00134','PF08613', 'CL0065','PF09080','PF16899','PF02234','PF12214')
names_cell_cycles_pfams <- c('Cyclin_N', 'Cyclin','Cyclin like','K-cyclin_vir_C','Cyclin_C_2','CDI','TPX2_importin')

clock_genes_pfams <- c('PF06203',  'PF07536')
names_clock_genes <- c('CCT (TOC1/CO-like)',  'HWE_HK')

all_pfams <- c(photosynthesis_pfams, nitrogen_pfams, carbon_fixation_pfams, iron_functions,
               other_functions, temperature_pfams, cell_cycle_pfams, clock_genes_pfams)
names_all_pfams <- c(names_photosynthesis_pfams, names_nitrogen_pfams, names_carbon_fixation_pfams, names_iron_function,
                     names_other_function, names_temperature_pfams, names_cell_cycles_pfams, names_clock_genes)

mg_data <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
pvals_g <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
rownames(mg_data) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(mg_data) <- taxos
rownames(pvals_g) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(pvals_g) <- taxos
mt_data <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
pvals_t <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
rownames(mt_data) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(mt_data) <- taxos
rownames(pvals_t) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(pvals_t) <- taxos
me_data <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
pvals_e <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
rownames(me_data) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(me_data) <- taxos
rownames(pvals_e) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(pvals_e) <- taxos

ut_data <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
pvals_ut <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
rownames(ut_data) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(ut_data) <- taxos
rownames(pvals_ut) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(pvals_ut) <- taxos

ug_data <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
pvals_ug <- matrix(NA, nrow = length(all_pfams), ncol=length(aldex_files))
rownames(ug_data) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(ug_data) <- taxos
rownames(pvals_ug) <- paste(all_pfams, ':', names_all_pfams, sep='')
colnames(pvals_ug) <- taxos

c=1
for (f in aldex_files){
  
  v <- readRDS(f)
  tx <- taxos[c]
  if (tx=='All'){
    tx <-'groups3'
  }
  vu <- readRDS(paste('../data/unigenes_',tx,'_mTG_',frac,'.rds', sep=''))
  
  mg_data[, c] <- v[[4]]$diff.btw[match(all_pfams, rownames(v[[4]]))]
  mt_data[, c] <- v[[2]]$diff.btw[match(all_pfams, rownames(v[[2]]))]
  me_data[, c] <- v[[6]]$diff.btw[match(all_pfams, rownames(v[[6]]))]
  pvals_g[,c] <- v[[4]]$we.eBH[match(all_pfams, rownames(v[[4]]))]<0.05
  pvals_t[,c] <- v[[2]]$we.eBH[match(all_pfams, rownames(v[[2]]))]<0.05
  pvals_e[,c] <- v[[6]]$we.eBH[match(all_pfams, rownames(v[[6]]))]<0.05
  
  vls <- vu[match(all_pfams, rownames(vu)) , 8]
  signs <- sign(vls)
  lg_vls <- log2(abs(vls))
  dvals <- lg_vls*signs
  ug_data[, c] <- dvals
  vls <- vu[match(all_pfams, rownames(vu)) , 4]
  signs <- sign(vls)
  lg_vls <- log2(abs(vls))
  dvals <- lg_vls*signs
  ut_data[, c] <- dvals
  pvals_ug[,c] <- vu[match(all_pfams, rownames(vu)), 7]<0.05
  pvals_ut[,c] <- vu[match(all_pfams, rownames(vu)), 3]<0.05
  
  c=c+1
}
sel0 <- !is.na(mg_data[,1]) & !is.na(mt_data[,1]) 
pvals_g<- pvals_g[sel0,]
pvals_t <- pvals_t[sel0,]
pvals_e <- pvals_e[sel0,]
pvals_ug<- pvals_ug[sel0,]
pvals_ut <- pvals_ut[sel0,]

photosynthesis_pfams1 <- photosynthesis_pfams[sel0[1:length(photosynthesis_pfams)]]
nit <- all_pfams %in% nitrogen_pfams
nitrogen_pfams1 <- nitrogen_pfams[sel0[nit]]
cab <- all_pfams %in% carbon_fixation_pfams
carbon_fixation_pfams1<- carbon_fixation_pfams[sel0[cab]]
ir <- all_pfams %in% iron_functions
iron_functions1<- iron_functions[sel0[ir]]
ot <- all_pfams %in% other_functions
other_functions1<- other_functions[sel0[ot]]
temp <- all_pfams %in% temperature_pfams
temperature_pfams1<- temperature_pfams[sel0[temp]]
cc <- all_pfams %in% cell_cycle_pfams
cell_cycle_pfams1<- cell_cycle_pfams[sel0[cc]]
cr <- all_pfams %in% clock_genes_pfams
clock_genes_pfams1<- clock_genes_pfams[sel0[cr]]
all_pfams1 <- all_pfams[sel0]
mg_data <- mg_data[sel0,]
mt_data <- mt_data[sel0,]
me_data <- me_data[sel0,]
ug_data <- ug_data[sel0,]
ut_data <- ut_data[sel0,]

set <- colorRampPalette(colors = c("darkblue", "white", 'darkorange'))(50)
set_2 <- colorRampPalette(colors = c("chartreuse", "white", 'darkred'))(50)



Functionsa <- c(rep('photosynthesis', length(photosynthesis_pfams1)),rep('Nitrogen metabolism', length(nitrogen_pfams1)) ,
           rep('Carbon fixation', length(carbon_fixation_pfams1)),rep('Iron metabolism', length(iron_functions1)),
           rep('Nutrient transporters', length(other_functions1)), rep('Cold acclimation', length(temperature_pfams1)),
           rep('Cell cycle', length(cell_cycle_pfams1)), rep('Circadian rythms', length(clock_genes_pfams1)))
Functionsa <- as.data.frame(Functionsa)
rownames(Functionsa)<- rownames(mg_data)
annot_g<- matrix(ifelse(pvals_g==T, "*", ""), nrow=dim(mg_data)[1])
annot_g[is.na(annot_g) | is.na(mg_data)] <- ""
annot_t<- matrix(ifelse(pvals_t==T, "*", ""), nrow=dim(mg_data)[1])
annot_t[is.na(annot_t) | is.na(mt_data)] <- ""
annot_e<- matrix(ifelse(pvals_e==T, "*", ""), nrow=dim(mg_data)[1])
annot_e[is.na(annot_e) | is.na(me_data)] <- ""

annot_ug<- matrix(ifelse(pvals_ug==T, "*", ""), nrow=dim(mg_data)[1])
annot_ug[is.na(annot_ug) | is.na(ug_data)] <- ""
annot_ut<- matrix(ifelse(pvals_ut==T, "*", ""), nrow=dim(mg_data)[1])
annot_ut[is.na(annot_ut) | is.na(ut_data)] <- ""

cols <- list(Functions=scales::hue_pal()(length( unique(Functionsa$Functions))))
names(cols$Functions)<- unique(Functionsa$Functions)

i=1
plot_list <- rep(list(NULL), length(unique(Functionsa$Functionsa))*3)
for (f in unique(Functionsa$Functionsa)){
  sel <- Functionsa==f
  Functions <- Functionsa[sel,]
  Functions <- as.data.frame(Functions)
  rownames(Functions)<- rownames(mg_data)[sel]
  Functions$Functions<- as.factor(Functions$Functions)
  a<- pheatmap(mg_data[sel,],cluster_cols = F,cluster_rows = F, color=set,
           breaks = seq(-max(abs(c(mg_data[sel,], mt_data[sel,],me_data[sel,] )), na.rm=T), 
                        max(abs(c(mg_data[sel,], mt_data[sel,],me_data[sel,] )), na.rm=T),length.out =  51),
           display_numbers =annot_g[sel, ], fontsize_number = 17, annotation_names_row = F,  
           show_rownames = F, legend = F)
  b<- pheatmap(mt_data[sel,],cluster_cols = F,cluster_rows = F, color=set,
           breaks = seq(-max(abs(c(mg_data[sel,], mt_data[sel,],me_data[sel,] )), na.rm=T), 
                        max(abs(c(mg_data[sel,], mt_data[sel,],me_data[sel,] )), na.rm=T),length.out =  51),
           display_numbers =annot_t[sel, ], fontsize_number = 17, annotation_names_row = F, 
           show_rownames = F, legend = F)
  c <- pheatmap(me_data[sel,],cluster_cols = F,cluster_rows = F, color=set,
           breaks = seq(-max(abs(c(mg_data[sel,], mt_data[sel,],me_data[sel,] )), na.rm=T), 
                        max(abs(c(mg_data[sel,], mt_data[sel,],me_data[sel,] )), na.rm=T),length.out =  51),
           display_numbers =annot_e[sel, ], fontsize_number = 17, annotation_names_row = F, legend = F, show_rownames = F)
  dt <- mt_data[sel,1:2]
  dt[,1:2]<-0
  d <- pheatmap(dt,cluster_cols = F,cluster_rows = F,annotation_row = Functions, color=set,
                breaks = seq(-max(abs(c(mg_data[sel,], mt_data[sel,],me_data[sel,] )), na.rm=T), 
                             max(abs(c(mg_data[sel,], mt_data[sel,],me_data[sel,] )), na.rm=T),length.out =  51),
                annotation_colors =cols,show_rownames = T, show_colnames = T,
                 fontsize_number = 17, annotation_names_row = F)
  
  e <- pheatmap(ug_data[sel,],cluster_cols = F,cluster_rows = F, color=set_2,
                breaks = seq(-max(abs(c(ug_data[sel,], ut_data[sel,] )), na.rm=T), 
                             max(abs(c(ug_data[sel,], ut_data[sel,] )), na.rm=T),length.out =  51),
                display_numbers =annot_ug[sel, ], fontsize_number = 17, annotation_names_row = F,  
                show_rownames = F, legend = F)
  f <- pheatmap(ut_data[sel,],cluster_cols = F,cluster_rows = F, color=set_2,
                breaks = seq(-max(abs(c(ug_data[sel,], ut_data[sel,] )), na.rm=T), 
                             max(abs(c(ug_data[sel,], ut_data[sel,] )), na.rm=T),length.out =  51),
                display_numbers =annot_ut[sel, ], fontsize_number = 17, annotation_names_row = F,  
                show_rownames = F, legend = F)
  dt <- ut_data[sel,1:2]
  dt[,1:2]<-0
  g <- pheatmap(dt,cluster_cols = F,cluster_rows = F,annotation_row = Functions, color=set_2,
                breaks = seq(-max(abs(c(ug_data[sel,], ut_data[sel,] )), na.rm=T), 
                             max(abs(c(ug_data[sel,], ut_data[sel,] )), na.rm=T),length.out =  51),
                annotation_colors =cols,show_rownames = T, show_colnames = T,
                fontsize_number = 17, annotation_names_row = F)
  
  plot_list[[i]] <- a[[4]]
  plot_list[[i+1]] <- b[[4]]
  plot_list[[i+2]] <- c[[4]]
  plot_list[[i+3]] <- d[[4]]
  plot_list[[i+4]] <- e[[4]]
  plot_list[[i+5]] <- f[[4]]
  plot_list[[i+6]] <- g[[4]]
  i=i+7
}
pdf(paste('pheatmap_pfams_aldex_', frac,'.pdf',sep=''), width = 25)
g=1
for (i in 1:length(unique(Functionsa$Functionsa))){
  grid.arrange(arrangeGrob(grobs= list(plot_list[[g]], plot_list[[g+1]],plot_list[[g+2]], plot_list[[g+3]]),ncol=4))
  g=g+7
}
dev.off()

pdf(paste('pheatmap_unigenes_pfams_aldex_', frac,'.pdf',sep=''), width = 15)
g=1
for (i in 1:length(unique(Functionsa$Functionsa))){
  grid.arrange(arrangeGrob(grobs= list(plot_list[[g+4]], plot_list[[g+5]], plot_list[[g+6]]),ncol=3))
  g=g+7
}
dev.off()


