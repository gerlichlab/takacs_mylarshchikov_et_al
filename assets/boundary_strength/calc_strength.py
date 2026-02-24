from sys import argv

import bioframe as bf
import cooler
import cooltools
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

CLR_URI = argv[1]
FLANK = int(argv[2])
IGNORE_DIAGS = int(argv[3])
OUT_PATH = argv[4]
CORES = 8

clr = cooler.Cooler(CLR_URI)
WIDTH_BINS = (2 * FLANK) // clr.binsize + 1


def agg_stack(stack):
    return np.nansum(np.nansum(np.triu(stack), axis=2), axis=1)


def get_gw_pup(clr):
    bins = clr.bins()[:][['chrom', 'start', 'end']]
    chrom_data = {chrom: agg_stack(cooltools.pileup(clr,
                                                    chrom_df,
                                                    flank=FLANK,
                                                    min_diag=IGNORE_DIAGS,
                                                    nproc=CORES))
                  for chrom, chrom_df in tqdm(bins.groupby('chrom'))
                  if chrom_df.shape[0] >= WIDTH_BINS}

    agg_holder = list()
    for chrom, chrom_df in bins.groupby('chrom'):
        arr = chrom_data.get(chrom)
        if arr is None:
            continue
        else:
            arr_df = chrom_df.assign(value=arr)
            agg_holder.append(arr_df)
    agg_df = pd.concat(agg_holder, ignore_index=True)
    return agg_df


def get_insulation(clr):
    ins_df = cooltools.insulation(clr,
                                  [FLANK],
                                  nproc=CORES,
                                  ignore_diags=IGNORE_DIAGS,
                                  append_raw_scores=True)
    return ins_df


hg19_chromsizes = bf.fetch_chromsizes('hg19')

pup_df = get_gw_pup(clr)

diamond_col = f"sum_balanced_{FLANK}"
ins_df = get_insulation(clr)[['chrom', 'start', 'end', diamond_col]]

table = pd.merge(ins_df, pup_df, on=['chrom', 'start', 'end'], how='inner')
table['pup'] = np.where(table['value'] == 0, np.nan, table['value'])
table['diamond'] = np.where(table[diamond_col] == 0, np.nan, table[diamond_col])

with np.errstate(divide='ignore', invalid='ignore'):
    table['strength'] = -np.log2(table['diamond'] / table['pup'])

clean_table = bf.trim(table, hg19_chromsizes).dropna().reset_index(drop=True)
bf.to_bigwig(clean_table, hg19_chromsizes, OUT_PATH, value_field='strength', engine='bigtools')
