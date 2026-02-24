#!/usr/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mincpus=1
#SBATCH --time=0-01:00:00
#SBATCH --mem=20gb
#SBATCH --output=logs/tads_%x_%a.out
#SBATCH --error=logs/tads_%x_%a.err
#SBATCH --array=0-5

export CONTAINER="/groups/gerlich/labinfo/Methods/singularity_images/gerlichlab-replchromconf-jupyter-v0.1_master_v4.5.0.sif"

singularity run $CONTAINER /bin/bash call_tads.sh $SLURM_ARRAY_TASK_ID
