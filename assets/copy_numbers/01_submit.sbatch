#!/usr/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mincpus=8
#SBATCH --time=0-01:00:00
#SBATCH --mem=20gb
#SBATCH --output=logs/cnvs.out
#SBATCH --error=logs/cnvs.err

export CONTAINER="/groups/gerlich/labinfo/Methods/singularity_images/gerlichlab-replchromconf-jupyter-v0.2_master.sif"
# Create cache folder and download required files b/c neoloopfinder can't do that inside singularity
mkdir -p .cache
wget -O .cache/hg19_mappability_100mer.1kb.bw -L http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19_mappability_100mer.1kb.bw
wget -O .cache/hg19_1kb_GC.bw -L https://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19_1kb_GC.bw
wget -O .cache/hg19.DpnII.npz -L http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.MboI.npz

singularity run $CONTAINER /bin/bash call_cnvs.sh
