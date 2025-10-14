This repository contains code and data necessary to replicate the study *Changes in gene expression in eukaryotic phytoplankton at the Atlantic-Arctic polar front*  
  
Preprint: https://www.biorxiv.org/content/10.1101/2022.11.01.514737v2  
  
The catalog of unigenes (fasta and abundance files) of MATOU v1.5 and the MGT v1.5 are available at https://www.genoscope.cns.fr/tara/

Large .rds files necessary to generate the figures are available on Zenodo at: https://doi.org/10.5281/zenodo.17316416. They must be placed in the folder `Main_analysis/data` to be able to regenerate the figures.

The repository contains two main folders:
- `Main_analysis/`: code and data to reproduce the figures
- `Preprocessing_data/`: code to preprocess data: generate data necessary to reproduce the figures (no need to be run to reproduce the figures as data are available in the `Main_analysis/data/` folder

All preprocessing and analysis were run on the Inti machine at Genoscope

- R dependencies:  
  
treemap  
gplots  
FactoMineR  
stringr  
RColorBrewer  
parallel  
stringr  
mapproj  
mapplots  
maptools  
SDMTools  
ncdf4  
scales  
bestglm  
reshape2  
gbm  
dismo  
dplyr  
parallel  
viridis  
viridisLite  
cluster  
ClusterR  
ggplot2  
tidyverse  
ggpointdensity  
paletteer  
pals  
MASS  
plyr  
shipunov  
data.table  
sm    
ALDEx2  
tidyr  
matlab  
gridExtra  
ggrepel  
raster  
pals  
viridis  

- Python dependencies  

matplotlib  
phate  
pandas  
numpy  
scprep  
sys  
seaborn  



