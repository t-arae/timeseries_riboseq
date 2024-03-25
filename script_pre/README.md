# `timeseries_riboseq/script_pre`

Scripts for the pre-processing/analyzing data.

The shell scripts listed bellow (and example codes in HTML reports) are assumed to be run in the directory where data analysis will be performed.

The HTML reports show example scripts of pre-processing/data analysis for each step.
You can see the report by download the HTML file and opening it in your web browser.
To run the example scripts, you have to prepare the appropriate runtime environments of R or Julia language.

## Contents

- `00_download_fastq_files.sh` Download FASTQ files from [DDBJ DRA](https://ddbj.nig.ac.jp/resource/sra-submission/DRA016840).
- `01_quality_filtering.sh` Filter reads by quality scores etc..
- `02_adaptor_trimming.sh` Trim adaptor from reads.
- `03_umi_trimming.sh` Trim UMI from reads and add it as a suffix of read name.
- `04_remove_ncRNA.sh` Remove reads mapped to the ncRNA sequences.
- `05_remove_marker_contami.sh` Remove reads mapped to the size marker.
- `06_mapping_by_star.sh` Align reads to the reference genome.
- `07_remove_softclip.html` (**HTML Report**) Remove soft-clipped bases from mapped reads.
- `08_umi_deduplication.sh` De-duplicate reads.
- `09_create_transcript_bed_file.sh` Make a transcript BED file from GTF file.
- `10_transfer_read_onto_transcripts.sh` Change genome mapped reads to the transcript coordinates.
- `11_riboWaltz.html` (**HTML Report**) Run Ribo-seq QC by `riboWaltz`.
- `12_filter_by_read_length.sh` Filter reads by its read length.
- `13_create_psite_bam.sh` Trim down reads to leave only p-site position.
