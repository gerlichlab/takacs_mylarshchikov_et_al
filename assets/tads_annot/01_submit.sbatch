#!/usr/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mincpus=1
#SBATCH --time=0-01:00:00
#SBATCH --mem=40gb
#SBATCH --output=logs/tads.out
#SBATCH --error=logs/tads.err

export CONTAINER="/groups/gerlich/labinfo/Methods/singularity_images/gerlichlab-replchromconf-jupyter-v0.1_master_v4.5.0.sif"

singularity run $CONTAINER /bin/bash call_tads.sh
