#!/bin/bash
#SBATCH -p batch
#SBATCH --mem=10G
#SBATCH -t 2:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=1-120

##cluster specific python 3 loading
##python environment requirements: pandas, numpy
module load miniforge3/25.3.0-3
source ${MAMBA_ROOT_PREFIX}/etc/profile.d/conda.sh

python window_H.py ${SLURM_ARRAY_TASK_ID} 









