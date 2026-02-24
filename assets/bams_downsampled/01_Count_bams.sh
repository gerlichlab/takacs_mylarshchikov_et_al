#! /usr/bin/env bash
#SBATCH --job-name="count_bams"
#SBATCH -c 2
#SBATCH --mem 16G
#SBATCH --array=1-6
#SBATCH --qos=short
#SBATCH --time=0-00:30:00

module load build-env/f2022
module load samtools/1.15-gcc-11.2.0

# set paths
seqDir="/groups/gerlich/experiments/Experiments_006500/006575/bam_repo/RAD21"
outDir="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/bams_downsampled"

case ${SLURM_ARRAY_TASK_ID} in
	1) path="WT_G1/bio_rep1" condition="RAD21_WT_G1_BR1";;
	2) path="WT_G1/bio_rep2" condition="RAD21_WT_G1_BR2";;
	3) path="WT_G2/bio_rep1" condition="RAD21_WT_G2_BR1";;
	4) path="WT_G2/bio_rep2" condition="RAD21_WT_G2_BR2";;
	5) path="WT_G1" condition="RAD21_WT_G1";;
	6) path="WT_G2" condition="RAD21_WT_G2";;
esac

# Count reads and append to consolidated file with sample name
count=$(samtools view -@ 2 -c ${seqDir}/${path}/${condition}.bam)
echo -e "${condition}\t${count}" >> ${outDir}/bam_counts.txt