from subprocess import run
from pathlib import Path


clr_paths = {condition: f"../coolers_downsampled/{condition}.mcool"
              for condition in ('WT_G2', 'WT_G1', 'dSororin_G2', 'Prometa', 'WT_G2_BR1n2', 'WT_G2_BR3_ds_1n2')}
RESOLUTION = 20_000
for cond, path in clr_paths.items():
    clr_uri = path + f"::/resolutions/{RESOLUTION}"
    ignore_diags = 2
    flank = 500_000
    out_path = Path(f"downsampled/bs.{cond}.bw")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    run(f"sbatch -J {cond} calc_strength.sbatch {clr_uri} {flank} {ignore_diags} {out_path}", shell=True)
