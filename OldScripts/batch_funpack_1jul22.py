#!/bin/bash
######## Job Name: May24_20k ########
#SBATCH -J Jun6_funpack
######## Job Output File: May19_out.oJOBID ########
#SBATCH -o Jun6_funpack_out.o%j
######## Job Error File: May19_error.eJOBID ########
#SBATCH -e Jun6_funpack_error.e%j
######## Number of nodes: 8 ########
#SBATCH -N 1
######## Number of tasks: 32 ########
#SBATCH -n 1
######## Memory per node: 20 GB ########
#SBATCH --mem 100G
######## Walltime: 100 hours ########
#SBATCH -t 5:00:00
 
######## Load module environment required for the job ########
module load python
source activate kayla_env

######## Run the job ########
python /scratch/khannon/funpack_extract_Jul1
