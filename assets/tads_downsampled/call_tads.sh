# active conda env in the container
source /srv/conda/etc/profile.d/conda.sh
conda activate notebook

# add OnTAD to PATH
export PATH="/users/dmitry.mylarshchikov/software/OnTAD/src:${PATH}"

# define current working folder (we need absolute paths)
export OUTFOLDER="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/tads_downsampled"

case $1 in
    0)
        export INPUT_NAME="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/coolers_downsampled/WT_G2.mcool"
        export OUTPUT_NAME="${OUTFOLDER}/nested_tads.WT_G2.bedpe"
        ;;
    1)
        export INPUT_NAME="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/coolers_downsampled/dSororin_G2.mcool"
        export OUTPUT_NAME="${OUTFOLDER}/nested_tads.dSororin_G2.bedpe"
        ;;
    2)
        export INPUT_NAME="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/coolers_downsampled/WT_G1.mcool"
        export OUTPUT_NAME="${OUTFOLDER}/nested_tads.WT_G1.bedpe"
        ;;
    3)
        export INPUT_NAME="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/coolers_downsampled/Prometa.mcool"
        export OUTPUT_NAME="${OUTFOLDER}/nested_tads.Prometa.bedpe"
        ;;
    4)
        export INPUT_NAME="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/coolers_downsampled/WT_G2_BR1n2.mcool"
        export OUTPUT_NAME="${OUTFOLDER}/nested_tads.WT_G2_BR1n2.bedpe"
        ;;
    5)
        export INPUT_NAME="/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/coolers_downsampled/WT_G2_BR3_ds_1n2.mcool"
        export OUTPUT_NAME="${OUTFOLDER}/nested_tads.WT_G2_BR3.bedpe"
        ;;
esac

# run ontad
cooler_ontad --short_name --binsize 20000 --maxsz 200 --minsz 2 --penalty 0.1 --o $OUTPUT_NAME $INPUT_NAME
