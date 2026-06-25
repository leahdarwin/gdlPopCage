#!/bin/bash
#SBATCH -p batch
#SBATCH --mem=10G
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1

##output directory
glm_dir="output/glm/"
output_file="${glm_dir}${model}.glm"

##directory containing sync file (not included in data upload too large, can be regenerated from raw reads)
sync_dir="${base_dir}aligned_reads_6.32/" 
input_file="${sync_dir}joined.sync"

model="treatment_time_repl"
script="poolFreqDiff_${model}.py"

##orginal (NOT MINE) and modified code from poolFreqDiff
##source: https://github.com/RAWWiberg/poolFreqDiff/ 
code_dir="tools/poolFreqDiff/"

##number of pool sequencing samples from selection experiment
nsamps=120

##more compute cluster specific r and python loads
##yaml for python environment is given in the github
eval "$(conda shell.bash hook)"
module load r
module load miniforge3/25.3.0-3
source ${MAMBA_ROOT_PREFIX}/etc/profile.d/conda.sh
conda activate py27

python "${code_dir}${script}" -filename "$input_file" -npops $nsamps -nlevels 1 -n 200 -mincnt 300 -minc 200 -maxc 100000 -rescale nr -zeroes 1 > "${output_file}.rin"

Rscript "${output_file}.rin" > "$output_file"