#!/bin/bash

#SBATCH --job-name="tlezuo_extract_hf"
#SBATCH --nodes=2
#SBATCH --output="job.out"
#SBATCH --time=24:00:00
#SBATCH --partition=postproc
#SBATCH --cpus-per-task=12

python calc_htd_icon_std.py
