source /srv/conda/etc/profile.d/conda.sh
conda activate notebook

python get_pileups.py 10_000 1_000_000 ./flank_1Mb_res_10Kb/
