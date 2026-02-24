source /srv/conda/etc/profile.d/conda.sh
conda activate notebook
export COOLER=$1
export FLANK=$2
export IGNORE_DIAGS=$3
export OUTPUT_PATH=$4
python calc_strength.py $COOLER $FLANK $IGNORE_DIAGS $OUTPUT_PATH