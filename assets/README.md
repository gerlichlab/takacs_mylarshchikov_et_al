# Assets
This folder contains code to create supporting datasets used in data analysis from processed ChIP-seq and Hi-C data.

## ChIP-seq data

### [bams_downsampled](/assets/bams_downsampled)
Precalculation of ChIP-seq signal tracks in bigWig file format from downsampled BAM files for Supplementary Figure S1F.

### [reproducible_peaks](/assets/reproducible_peaks/)
Annotation of reproducible ChIP-seq peaks for Figure 1E-F.

### [copy_numbers](/assets/copy_numbers)
Estimation of copy numbers from Hi-C data for normalization of ChIP-seq tracks. Applied to average profiles in Figures 2D, 3B, 4B, 4C, 6B, 6C, Supplementary Figure S6C.

## Hi-C data
### [coolers_downsampled](/assets/coolers_downsampled)
Downsampled Hi-C contact maps in cooler format for Figures 1A, 1B, 1C.

### [boundary_strength/downsampled](/assets/boundary_strength/downsampled)
Precalculation of insulation score tracks from downsampled cooler files for Figure 1B.

### [tads_downsampled](/assets/tads_downsampled)
TAD calling on downsampled coolers for Figure 1C.

### [tads_annot](/assets/tads_annot)
Annotation of TADs in WT G2 with OnTAD and extraction of TAD boundaries. First used in Figure 2C and then throughout the paper.

### [pileups](/assets/pileups)
Precalculation of average Hi-C contact maps, or pileups. Used in Figures 2E, 3C, 3D, 5A, 5B, 5C, 5D, Supplementary Figure S2.

### [boundary_strength](/assets/boundary_strength/downsampled)
Precalculation of insulation score tracks for Figure 3E.

### [scalings](/assets/scalings)
Precalculation of scaling data from cis- and trans-sister maps for Supplementary Figure 3K.

## How to run
Code in assets can be run with a dedicated Singularity image [gerlichlab/replchromconf-jupyterlab](https://github.com/gerlichlab/replchromconf-jupyterlab). Currently, OnTAD for TAD calling needs to be installed separately.