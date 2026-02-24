#! /usr/bin/env bash
#SBATCH --job-name="Downsample_bams"
#SBATCH -c 4
#SBATCH --mem 20G
#SBATCH --array=1-2
#SBATCH --qos=short
#SBATCH --time=0-03:00:00
#SBATCH --output=/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/bams_downsampled/logs/downsample_%a.out
#SBATCH --error=/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/bams_downsampled/logs/downsample_%a.err

module load build-env/f2022
module load samtools/1.15-gcc-11.2.0

# Enable exit on error
set -euo pipefail

# set paths
seqDir="/groups/gerlich/experiments/Experiments_006500/006575/bam_repo/RAD21"
outDir="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/bams_downsampled/RAD21"

#Define target read count
target_reads=160000000

#Set seed for downsampling for reproducibility
seed=42

case ${SLURM_ARRAY_TASK_ID} in
	1) path="WT_G1" name="RAD21_WT_G1";;
	2) path="WT_G2" name="RAD21_WT_G2";;
	*) echo "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"; exit 1;;
esac

#Create output dir if it doesn't exist
mkdir -p ${outDir}/${path}

# Check if input BAM exists
if [[ ! -f ${seqDir}/${path}/${name}.bam ]]; then
    echo "Error: Input BAM file not found: ${seqDir}/${path}/${name}.bam"
    exit 1
fi

echo "=== Starting processing for ${name} ==="

# Calculate total reads
total_reads=$(samtools idxstats ${seqDir}/${path}/${name}.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}')
echo "Total reads in input: ${total_reads}"

# Check if target exceeds total
if (( target_reads > total_reads )); then
    echo "Error: Target reads (${target_reads}) exceed total reads (${total_reads}) in ${name}"
    exit 1
fi

# Calculate scaling factor
fraction=$(echo "scale=10; ${target_reads} / ${total_reads}" | bc)
echo "Calculated fraction: ${fraction}"

# Check if fraction is valid (should be <= 1)
if (( $(echo "${fraction} > 1" | bc -l) )); then
    echo "Error: Calculated fraction (${fraction}) is greater than 1"
    exit 1
fi

# Convert fraction to samtools format (remove leading 0.)
# samtools expects format like: seed.fraction (e.g., 42.5 for 50%)
trimmed_fraction=${fraction#0.}
trimmed_fraction=${trimmed_fraction#.}
echo "Trimmed fraction: ${trimmed_fraction}"

echo "Processing ${name}: target=${target_reads}, fraction=${fraction}, samtools_fraction=${trimmed_fraction}"

# Downsample to temp file
echo "Starting downsampling..."
if ! samtools view -@ 2 -bs ${seed}.${trimmed_fraction} ${seqDir}/${path}/${name}.bam > ${outDir}/${path}/${name}.temp.bam; then
    echo "Error: Downsampling failed for ${name}"
    exit 1
fi
echo "Completed downsampling: ${name}"

# Sort
echo "Starting sorting..."
if ! samtools sort -o ${outDir}/${path}/${name}.bam -@ 2 ${outDir}/${path}/${name}.temp.bam; then
    echo "Error: Sorting failed for ${name}"
    rm -f ${outDir}/${path}/${name}.temp.bam
    exit 1
fi
echo "Completed sorting: ${name}"

# Clean up temp file
rm ${outDir}/${path}/${name}.temp.bam

# Count reads
echo "Counting reads in final BAM..."
read_count=$(samtools view -c -@ 2 ${outDir}/${path}/${name}.bam)
echo "Final BAM contains: ${read_count} reads (target was ${target_reads})"

# Calculate difference
diff=$((target_reads - read_count))
percent_diff=$(echo "scale=2; (${diff} / ${target_reads}) * 100" | bc)
echo "Difference from target: ${diff} reads (${percent_diff}%)"

# Index
echo "Starting indexing..."
if ! samtools index -@ 2 ${outDir}/${path}/${name}.bam; then
    echo "Error: Indexing failed for ${name}"
    exit 1
fi
echo "Completed indexing: ${name}"

echo "=== ${name} processing complete ==="