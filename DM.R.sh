#!/bin/bash
#
#SBATCH --job-name=DM.R_Rett
#SBATCH --workdir /share/lasallelab/Ben/
#SBATCH --ntasks=2 # Number of cores/threads
#SBATCH --mem=64000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --output=DM.R_%A.out # File to which STDOUT will be written
#SBATCH --error=DM.R_%A.err # File to which STDERR will be written
#SBATCH --time=4-00:00:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
# Last Update Date: 09-20-2018
# Version: 0.98
#
# DMR inference and data visualization for CpG_Me output and bismark cytosine reports
#
# If you use this, please cite:
##########################################################################################

##############
# Initialize #
##############

# 

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
/share/lasallelab/programs/CpG_Me/DM.R"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
