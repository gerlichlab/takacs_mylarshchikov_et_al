from sys import argv
import pairtools.lib.scaling as scaling
import bioframe as bf

input_path = argv[1]
output_path = argv[2]

hg19_chromsizes = bf.fetch_chromsizes('hg19')
hg19_centromeres = bf.fetch_centromeres('hg19')
hg19_arms = bf.make_chromarms(hg19_chromsizes, hg19_centromeres)


def calc_scaling(path):
    return scaling.compute_scaling(path,
                                   regions=hg19_arms,
                                   chromsizes=hg19_chromsizes,
                                   dist_range=(1, 1_000_000_000), 
                                   n_dist_bins_decade=8,
                                   ignore_trans=True,
                                   nproc_in=8,
                                   chunksize=int(1e6))[0]


calc_scaling(input_path).to_parquet(output_path)
