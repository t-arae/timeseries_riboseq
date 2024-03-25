# `timeseries_riboseq`

This repository provides codes for reproducing results in the article titled "Impact of translational regulation on diel expression revealed by time-series ribosome profiling in *Arabidopsis*".

## Contents

```
timeseries_riboseq
  ├── README.md       This file
  ├── script_julia    Julia scripts for removing soft clipped bases from BAM file
  ├── script_pre      Shell scripts for data pre-processing (and two html files,
  │                   which include code and output of pre-process.)
  ├── script_r        R scripts for provide some settings
  └── reports         Analyses reports, which include R code to reproduce results
```

## Requirements

- Prepare an environment which installed programs listed bellow.
  * `fastp 0.23.2`
  * `featureCounts v2.0.3`
  * `cutadapt 4.1`
  * `bedtools v2.30.0`
  * `bowtie2 version 2.4.4`
  * `STAR 2.7.9a`
  * `UMI-tools version: 1.1.2`
  * `ribotricer, version 1.3.2`
  * `gtfToGenePred`, `genePredToBed` (can be get from [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/))
  * `ubu` (can be get from [here](https://github.com/mozack/ubu))
  * `empirical-JTK_CYCLE` (slightly modified version from original [one](https://github.com/alanlhutchison/empirical-JTK_CYCLE-with-asymmetry) can be get from [here](https://github.com/t-arae/empirical-JTK_CYCLE-with-asymmetry))
- Download the miscellaneous files from [Here](https://github.com/t-arae/timeseries_riboseq/releases/tag/misc_v1.0.0), decompress and put them to the appropriate directory.
- Prepare an Julia v1.9 environment for running the Julia script.
- Prepare an R v4.2.1 environment for running R scripts.
    * To easily restore the R analysis environment, this repository includes the `renv.lock` file.
    * You can install the most packages which I used in this project, by using this `renv.lock` file and `renv::restore()` (Some packages may need to be installed manually).
    * I also use [the self-made R package](https://github.com/t-arae/ngsmisc) in my R scripts. 

## Reference
Aoyama, H., Arae, T., Yamashita, Y., Toyoda, A., Naito, S., Sotta, N., Chiba, Y., 2024. Impact of translational regulation on diel expression revealed by time‐series ribosome profiling in Arabidopsis.  Plant J. https://doi.org/10.1111/tpj.16716
