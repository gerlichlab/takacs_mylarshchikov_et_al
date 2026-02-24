#!/usr/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mincpus=8
#SBATCH --time=0-08:00:00
#SBATCH --mem=70gb
#SBATCH --output=logs/pups_5mb.out
#SBATCH --error=logs/pups_5mb.err

export CONTAINER="/groups/gerlich/labinfo/Methods/singularity_images/gerlichlab-replchromconf-jupyter-v0.2_master.sif"

singularity run $CONTAINER /bin/bash get_pups_5Mb.sh
