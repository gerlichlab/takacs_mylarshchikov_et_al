# active conda env in the container
source /srv/conda/etc/profile.d/conda.sh
conda activate notebook

# add OnTAD to PATH
export PATH="/users/dmitry.mylarshchikov/software/OnTAD/src:${PATH}"
# define input cooler path
export INPUT_NAME="/groups/gerlich/experiments/Experiments_006500/006575/coolers_repo/WT_G2/WT_G2.all.mcool"
# define current working folder (we need absolute paths)
export OUTFOLDER="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/tads_annot"
# define output path
export OUTPUT_NAME="${OUTFOLDER}/nested_tads.bedpe"
# run ontad
cooler_ontad --short_name --binsize 10000 --maxsz 400 --minsz 5 --penalty 0.1 --o $OUTPUT_NAME $INPUT_NAME
