import multiprocessing
from pathlib import Path
from sys import argv

import bioframe as bf
import cooler
import cooltools
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

# Parameters
CORES = 8
RES = int(argv[1])
FLANK = int(argv[2])
OUTDIR = Path(argv[3])

# Load boundaries
boundaries = pd.read_parquet('/groups/gerlich/labinfo/Papers/2025_ReplChromConf/analysis/assets/tads_annot/boundaries.parquet')

# Get coolers
clr_paths = {condition: {sis: f"/groups/gerlich/experiments/Experiments_006500/006575/coolers_repo/{condition}/{condition}.{sis}.mcool"
                         for sis in ('cis', 'trans', 'all')}
              for condition in ('WT_G2', 'dNIPBL_G2', 'dWAPL_G2', 'dCTCF_G2')}

# Get pileups
def get_pileups_stack(clr, regions_df, view_df, flank, nproc, expected=None):
    with multiprocessing.Pool(nproc) as p:
        stack = cooltools.pileup.__wrapped__(clr,
                                             regions_df,
                                             view_df=view_df,
                                             expected_df=expected,
                                             flank=flank,
                                             map_functor=p.map)
    return stack


hg19_cens = bf.fetch_centromeres('hg19')

for cond, sis_data in tqdm(clr_paths.items(), desc='getting pups'):
    if cond != "dCTCF_G2":
        continue
    for sis, path in tqdm(sis_data.items(), desc=cond):
        if sis != "all":
            continue
        clr = cooler.Cooler(path + f"::/resolutions/{RES}")
        chromsizes = clr.chroms()[:].rename(columns={'name': 'chrom'})
        arms = bf.make_chromarms(chromsizes, hg19_cens)
        for kind in tqdm(('obs', 'obsexp'), desc=sis):
            if kind == 'obsexp':
                expected_data = pd.read_parquet(f'./expected/expected.{cond}.{sis}.{RES}.parquet')
            else:
                expected_data = None
            for level, level_df in tqdm(boundaries.groupby('levelGroup'), desc=kind):
                    pileup_stack = get_pileups_stack(clr,
                                                    level_df,
                                                    arms,
                                                    FLANK,
                                                    CORES,
                                                    expected_data)
                    pileup = np.nanmean(pileup_stack, axis=0)
                    np.save(OUTDIR / f"{cond}.{sis}.{level}.{kind}.pileup.npy", pileup, allow_pickle=False)
