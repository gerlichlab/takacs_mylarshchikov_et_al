source /srv/conda/etc/profile.d/conda.sh
conda activate notebook

python get_pileups.py 20_000 5_000_000 ./flank_5Mb_res_20Kb/
