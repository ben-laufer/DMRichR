#!/bin/bash
#
#SBATCH --job-name=DM.R
#SBATCH --chdir /share/lasallelab/
#SBATCH --ntasks=20 # Number of cores/threads
#SBATCH --mem=192000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --time=4-00:00:00

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

module load R/3.6.1 

########
# DM.R #
########

call="Rscript \
--vanilla \
/share/lasallelab/programs/DMRichR/DM.R \
--genome hg38 \
--coverage 1 \
--perGroup '1' \
--minCpGs 5 \
--maxPerms 10 \
--cutoff '0.05' \
--testCovariate Diagnosis \
--adjustCovariate 'Sex;Age' \
--cores 20"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
