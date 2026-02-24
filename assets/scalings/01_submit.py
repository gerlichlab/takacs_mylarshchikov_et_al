from pathlib import Path
from subprocess import run


repo = Path('/groups/gerlich/experiments/Experiments_006500/006575/pairs_repo/')
conditions = ('WT_G2', 'dCTCF_G2', 'dNIPBL_G2', 'dWAPL_G2')

for cond in conditions:
    for sis in ('cis', 'trans'):
        input_path = repo / cond / f"{cond}.{sis}.pairs.gz"
        output_path = Path(f'./{cond}/{cond}.{sis}.scaling.parquet')
        output_path.parent.mkdir(parents=True, exist_ok=True)
        jobname = f"{cond}_{sis}_scaling"
        run(f"sbatch -J {jobname} calc_scaling.sbatch {input_path} {output_path}", shell=True)