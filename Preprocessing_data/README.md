This folder contains all codes to preprocess data necessary to generate the figures from the `Main_analysis/` folder:

## 1. Download raw datasets of MATOU v1.5 and MGT v1.5 at: https://www.genoscope.cns.fr/tara/

## 2. Parse the file to select stations of interest
- run `Rscript extract_data_seq.R GGZZ 0`  
This will create 792 files GGZZmetaTnMetaG_$i.rds each containing a subset of the dataset (158M unigenes)  
- create a directory to store them: `mkdir subset_metat_GGZZ/`
- `mv GGZZmetaTnMetaG_*rds subset_metat_GGZZ/`

## 3. save taxonomy as rds
- `Rscript save_taxonomy_rds.R` output: `taxID_uni.rds`

## 4. Calculate the clr
- calculate the geometric mean for each station: `Rscript geometric_mean_metaTG.R GGZZ`
- calculate the clr: `Rscript save_geometric_mean_metaTG.R GGZZ`
- `mv GGZZmetaTnMetaG_*clr.rds subset_metat_GGZZ/`

## 5. save expressed unigenes
- run `Rscript pre_MetaT_analysis_0.R GGZZ $i` i from 1 to 792. can be done using pegasus on an HPC: `./pre_submit_metat_analysis_0.sh GGZZ`
- `mv expressed_unilist*GGZZ*rds subset_metat_GGZZ/`

## 6. save taxonomy at the desired level
- run `Rscript save_data_unis_group.R groups3`, output: `data_uni_groups3.rds`
- 






