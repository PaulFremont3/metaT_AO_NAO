This folder contains all codes to preprocess data necessary to generate the figures from the `Main_analysis/` folder:

# 1. Download raw datasets of MATOU v1.5 and MGT v1.5 at: https://www.genoscope.cns.fr/tara/

# 2. Parse the file to select stations of interest
run `Rscript extract_data_seq.R GGZZ 0`  
This will create 792 files GGZZmetaTnMetaG_$i.rds each containing a subset of the dataset (158M unigenes)  
create a directory to store them: `mkdir subset_metat_GGZZ/`
`mv GGZZmetaTnMetaG_*rds subset_metat_GGZZ/` 




