#!/bin/bash
#
#SBATCH --job-name=DM.R
#SBATCH --workdir /share/lasallelab/Ben/
#SBATCH --ntasks=2 # Number of cores/threads
#SBATCH --mem=64000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --time=5-00:00:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

################
# Load Modules #
################

module load R

######################
# Set Up Environment #
######################

aklog

########
# DM.R #
########

call="Rscript \
--vanilla \
/share/lasallelab/programs/DMRichR/DM.R \
--genome hg38 \
--coverage 1 \
--testCovariate Diagnosis \
--adjustCovariate Age \
--matchCovariate Sex \
--cores 2"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
