source /srv/conda/etc/profile.d/conda.sh
conda activate notebook
export INPUT_PATH=$1
export OUTPUT_PATH=$2
python calc_scaling.py $INPUT_PATH $OUTPUT_PATH
