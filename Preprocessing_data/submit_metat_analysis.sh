#!/bin/bash
#MSUB -r metat_analysis
#MSUB -q normal
#MSUB -n 1
#MSUB -x
#MSUB -c 35
#MSUB -T 86400  # secondes
#MSUB -o pegasus_job_%j.out
#MSUB -e pegasus_job_%j.err

# gestion des erreurs
set -e;
set -o pipefail;
set -u;

# calcul le nombre de cœur de la réservation
TOTAL_CORE=$(sinfo -eo '%c %n' -n $SLURM_JOB_NODELIST | awk '{ SUM+=$1 } END {print SUM+1}') 

# chargement pegasus et configuration de l'env
module load pegasus/4.8.0.a

module load r
# lancement du pegasus
mpirun -oversubscribe -n $TOTAL_CORE pegasus-mpi-cluster ./comds_metat_${1}.pegasus
