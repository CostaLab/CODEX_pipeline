#!/bin/bash

#Request cores. In this case - 2 cores will be requested
#SBATCH --cpus-per-task=1

#Request memory per thread. Here it's 1 GB per requested core
#SBATCH --mem-per-cpu=300G

# Set a job name
#SBATCH --job-name=CODEXvis

#Specify time before the job is killed by scheduler (in case is hangs). In this case - 24 hours
#SBATCH --time=24:00:00 

#Declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.out

#Email
#SBATCH --mail-type=END
#SBATCH --mail-user=aryam950812@gmail.com

#Run something
#export PATH="/home/hu367653/miniconda3/bin:$PATH"

#conda activate CODEX_env
python ./codex_vis.py