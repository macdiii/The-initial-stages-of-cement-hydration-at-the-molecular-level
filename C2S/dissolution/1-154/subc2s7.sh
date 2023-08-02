#!/bin/bash
#SBATCH -p amd_512
#SBATCH -n 1
#SBATCH -c 4
source /public3/soft/modules/module.sh
module load anaconda/3-Python3.7.4-fenggl
source activate py38
python -u 7.py

