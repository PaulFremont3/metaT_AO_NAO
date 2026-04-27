This folder contains all codes to preprocess data necessary to generate the figures from the `Main_analysis/` folder. Note that some scripts run using the parallel library of R, you need to set the number of core used.

## 1. Download raw datasets of MATOU v1.5 and MGT v1.5 at: https://www.genoscope.cns.fr/tara/

## 2. Parse the file to select stations of interest
- run `Rscript extract_data_seq.R GGZZ 0`  
This will create 792 files GGZZmetaTnMetaG_$i.rds each containing a subset of the dataset (158M unigenes)  
- create a directory to store them: `mkdir subset_metat_GGZZ/`
- `mv GGZZmetaTnMetaG_*rds subset_metat_GGZZ/`

## 3. Save taxonomy as rds
- `Rscript save_taxonomy_rds.R` output: `taxID_uni.rds`

## 4. Calculate the clr
- calculate the geometric mean for each station: `Rscript geometric_mean_metaTG.R GGZZ`
- calculate the clr: `Rscript save_geometric_mean_metaTG.R GGZZ`
- `mv GGZZmetaTnMetaG_*clr.rds subset_metat_GGZZ/`

## 5. Save expressed unigenes
- run `Rscript pre_MetaT_analysis_0.R GGZZ $i` i from 1 to 792. can be done using pegasus on an HPC: `./pre_submit_metat_analysis_0.sh GGZZ`
- `mv expressed_unilist*GGZZ*rds subset_metat_GGZZ/`

## 6. Save taxonomy at the desired level
- run `Rscript save_data_unis_group.R groups3`, output: `data_uni_groups3.rds`
- save for the subset: run `Rscript save_rds_taxID_groups3.R GGZZ`
- `mv TaxID_groups3_*GGZZ*.rds subset_metat_GGZZ/`

## 7. Unigenes expressed in more than 4 stations
- run `Rscript pre_MetaT_analysis_1.R GGZZ $i` i from 1 to 792. can be done using pegasus on an HPC: `./pre_submit_metat_analysis_1.sh GGZZ`
- `mv pre_unilist*GGZZ*rds subset_metat_GGZZ/`

## 8. Classify unigenes based on basin
- run `Rscript pre_MetaT_analysis_6.R GGZZ $i` i from 1 to 792. can be done using pegasus on an HPC: `./pre_submit_metat_analysis_6.sh GGZZ`
- `mv class_unilist_6_*GGZZ*rds subset_metat_GGZZ/`

## 9. MATOU annotation pfams
- run `Rscript matou-rds.R`
- run `Rscript subset_matou.R GGZZ i1 i2` i1 and i2 allow to treat the 792 files separately: treat files from i1 to i2
- `mv matou*GGZZ*rds subset_metat_GGZZ/`

## 10. Run the main metat analysis (correlation of gene expression with environmental parameters)
- run `./submit_metat_analysis.sh GGZZ` (pegasus), output: `metat_analysis_GGZZ.txt`
- run `Rscript save_metat_analysis_rds.R GGZZ`, output: `metat_analysis_GGZZ.rds`

## 11. Create pfam X taxo X stations tables  
- create matou baseline (list of unique pfams): run `Rscript matou_baseline.R GGZZ`
- create pfam station table no taxo (bis) for each 792 files: `Rscript create_pfam_table_allpf_bis.R GGZZ`
- create pfam station table for each 792 files: run `Rscript create_pfam_table_allpf_taxo.R GGZZ groups3`
- create pfam station table (bis) for each 792 files: run `Rscript create_pfam_table_allpf_taxo_bis.R GGZZ groups3`
- create pfam station table atl arc and com for each 792 files: run `Rscript create_pfam_table_allpf_atl_vs_arc_taxo.R GGZZ 6 taxo_groups3`
- create pfam station table atl arc and com (bis) for each 792 files: run `Rscript create_pfam_table_allpf_atl_vs_arc_taxo_bis.R GGZZ 6 taxo_groups3`
- concatenate pfam station table no taxo (bis) (calculate the sum >1 for bis*): run `Rscript save_pfam_table_allpf_bis.R GGZZ 0`
- concatenate pfam station table: run `Rscript save_pfam_table_allpf.R GGZZ taxo_groups3`
- concatenate pfam station table (bis): run `Rscript save_pfam_table_allpf_bis.R GGZZ taxo_groups3`
- concatenate pfam station table atl arc and com: run `Rscript save_pfam_table_allpf_atl_vs_arc.R GGZZ 6 taxo_groups3`
- concatenate pfam station table atl arc and com (bis): run `Rscript save_pfam_table_allpf_atl_vs_arc_bis.R GGZZ 6 taxo_groups3`  
*not bis: abundance of pfam is divided by number of pfam of each unigene  
bis: abundance of all pfams are counted => sum of abundance >1 (then renormalized)

## 12. Create GO abundances tables and run analyssis of differential abundance
- pre processing of GO table: run `Rscript GO_table_pre_process.R`
- create taxo X GO table X station table: run `Rscript create_GO_table_allpf_allGO_taxo.R GGZZ 6 taxo_groups3 _bis`
- differential abundance analysis: run `Rscript GO_analysis_clusters.R taxo_groups3 _bis`
- save GO names: run `Rscript save_GO_representative_names.R`

## 13. Run the PHATE analysis
- prepare traing dataset: run `Rscript prepare_training_phate_expressions.R GGZZ 0 0 0 0`
- run the phate analysis: run `python phate_expressions.py $1 GGZZ 0 0 0 0 500` $1 is the number of cores used

## 14. Run the ALDEx analysis
- preprocess:  run `Rscript preprocess_ALDEx_taxo.R`
- number of unigenes comparison: `Rscript number_unigenes_comparison.R`
- run the analysis: run `Rscript ALDEx.R GGZZ`

## 15. Calculate the Bray-curtis index based on unigenes abundances
- run the calculation: `Rscript bray-curtis_unigenes.R GGZZ`






