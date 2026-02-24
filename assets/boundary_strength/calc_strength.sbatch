#!/usr/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mincpus=8
#SBATCH --time=0-01:00:00
#SBATCH --mem=10gb
#SBATCH --qos=rapid
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err

export CONTAINER="/groups/gerlich/labinfo/Methods/singularity_images/gerlichlab-replchromconf-jupyter-v0.2_master.sif"
export COOLER=$1
export FLANK=$2
export IGNORE_DIAGS=$3
export OUTPUT_PATH=$4
singularity run $CONTAINER /bin/bash calc_strength.sh $COOLER $FLANK $IGNORE_DIAGS $OUTPUT_PATH
