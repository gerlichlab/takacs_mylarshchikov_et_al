# active conda env in the container
source /srv/conda/etc/profile.d/conda.sh
conda activate notebook

MTX="/groups/gerlich/experiments/Experiments_006200/006288/coolers/WT/WT_G2.all.1000.mcool::/resolutions/100000"

# OUTFOLDER="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/copy_numbers"
OUTFOLDER="."
CACHEFOLDER="${OUTFOLDER}/.cache"
PROF_FILE="${OUTFOLDER}/cnv_profile.bedGraph"
PROF_LOG="${OUTFOLDER}/cnv_profile.log"
SEG_FILE="${OUTFOLDER}/cnv_seg.bedGraph"
SEG_LOG="${OUTFOLDER}/cnv_seg.log"
FIG_FILE="${OUTFOLDER}/cnv_plot.png"

calculate-cnv -H $MTX -g hg19 -e DpnII --output $PROF_FILE --logFile $PROF_LOG --cachefolder $CACHEFOLDER
segment-cnv --cnv-file $PROF_FILE --binsize 100000 --ploidy 3 --output $SEG_FILE --nproc 8 --logFile $SEG_LOG
plot-cnv --cnv-profile $PROF_FILE --cnv-segment $SEG_FILE --output-figure-name $FIG_FILE
