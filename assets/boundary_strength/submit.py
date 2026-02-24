from subprocess import run
from pathlib import Path


clr_paths = {condition: {sis: f"/groups/gerlich/experiments/Experiments_006500/006575/coolers_repo/{condition}/{condition}.{sis}.mcool"
                         for sis in ('cis', 'trans')}
              for condition in ('WT_G2', 'dNIPBL_G2', 'dWAPL_G2', 'dCTCF_G2')}
RESOLUTION = 10_000
for cond, sis_data in clr_paths.items():
    if cond != "dCTCF_G2":
        continue
    for sis, path in sis_data.items():
        clr_uri = path + f"::/resolutions/{RESOLUTION}"
        if sis == 'trans':
            ignore_diags = 15
        else:
            ignore_diags = 2
        flank = 500_000
        out_path = Path(f"{cond}/bs.{cond}.{sis}.bw")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        run(f"sbatch -J {cond}_{sis} calc_strength.sbatch {clr_uri} {flank} {ignore_diags} {out_path}", shell=True)
