# `timeseries_riboseq/script_pre`

Scripts for the pre-processing/analyzing data.

The shell scripts listed bellow (and example codes in PDF reports) are assumed to be run in the directory where data analysis will be performed.

The PDF reports show example scripts of pre-processing/data analysis for each step.
To run the example scripts, you have to prepare the appropriate runtime environments of R or Julia language.

## Contents

- `00_download_fastq_files.sh` Download FASTQ files from [DDBJ DRA](https://ddbj.nig.ac.jp/resource/sra-submission/DRA016840).
- `01_quality_filtering.sh` Filter reads by quality scores etc..
- `02_adaptor_trimming.sh` Trim adaptor from reads.
- `03_umi_trimming.sh` Trim UMI from reads and add it as a suffix of read name.
- `04_remove_ncRNA.sh` Remove reads mapped to the ncRNA sequences.
- `05_remove_marker_contami.sh` Remove reads mapped to the size marker.
- `06_mapping_by_star.sh` Align reads to the reference genome.
- `07_remove_softclip.pdf` (**PDF Report**) Remove soft-clipped bases from mapped reads.
- `08_umi_deduplication.sh` De-duplicate reads.
- `09_create_transcript_bed_file.sh` Make a transcript BED file from GTF file.
- `10_transfer_read_onto_transcripts.sh` Change genome mapped reads to the transcript coordinates.
- `11_riboWaltz.pdf` (**PDF Report**) Run Ribo-seq QC by `riboWaltz`.
- `12_filter_by_read_length.sh` Filter reads by its read length.
- `13_create_psite_bam.sh` Trim down reads to leave only p-site position.
- `14_1_uorf_detection_with_ribotricer.sh` Search uORF candidates.
- `14_2_make_ribotricer_uorf_annotation.pdf` (**PDF Report**) Make an uORF annotation file.
- `15_1_readcount_rna.sh` Read counting for RNA-seq reads.
- `15_2_readcount_morf.sh` Read counting for mORF mapped P-site reads.
- `15_3_readcount_each_uorf.sh` Read counting for each uORF mapped P-site reads.
- `15_4_merge_readcount.pdf` (**PDF Report**) Merge read-counts data into one.

- `16_1_deseq2_rna.pdf` (**PDF Report**) Differential mRNA abundance analysis with RNA-seq data.
- `16_2_deseq2_morf.pdf` (**PDF Report**) Differential translational activity analysis with Ribo-seq data (for mORF).

- `17_1_ribotricer_active_uorf_annotation.pdf` (**PDF Report**) Make an active uORF annotation file.
- `17_2_readcount_uorf.sh` Read counting of P-site reads mapped to the uORFs of each gene.
- `17_3_merge_readcount.pdf` (**PDF Report**) Merge read-counts data into one.

- `18_1_ejtk_w_filtering.pdf` (**PDF Report**) Prepare data for the analysis of `empirical-JTK_CYCLE`.
- `18_2_ejtk_rna.sh` `empirical-JTK_CYCLE` analysis for RNA-seq reads.
- `18_3_ejtk_morf.sh` `empirical-JTK_CYCLE` analysis for mORF mapped P-site reads.
- `18_4_ejtk_uorf.sh` `empirical-JTK_CYCLE` analysis for P-site reads mapped to the uORFs of each gene.

- `19_1_deseq2_uorf.pdf` (**PDF Report**) Differential translational activity analysis with Ribo-seq data (for uORF).
- `19_2_deseq2_te.pdf` (**PDF Report**) Differential translational efficiency analysis (for mORF).
- `19_3_deseq2_te_morf_uorf.pdf` (**PDF Report**) Differential translational efficiency analysis (for uORF).

- `20_create_psite_coverage.sh` Make P-site reads coverage files using `deeptools bamCoverage`.

- `99_create_summary_list.pdf` (**PDF Report**) Merge all data into one CSV file.
- `99_run_fastqc.sh` Run `FastQC` for pre-processed data.
