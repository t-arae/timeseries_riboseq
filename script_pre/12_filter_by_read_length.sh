#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/bam_filtered_len

samtools view -H data_preproc/umi_dedup/zt0_1_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt0_1_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt0_1_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt0_1_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt0_2_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt0_2_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt0_2_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt0_2_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt3_1_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt3_1_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt3_1_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt3_1_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt3_2_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt3_2_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt3_2_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt3_2_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt6_1_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt6_1_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt6_1_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt6_1_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt6_2_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt6_2_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt6_2_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt6_2_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt12_1_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt12_1_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt12_1_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt12_1_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt12_2_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt12_2_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt12_2_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt12_2_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt18_1_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt18_1_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt18_1_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt18_1_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt18_2_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt18_2_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt18_2_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt18_2_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt21_1_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt21_1_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt21_1_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt21_1_ribo.sort.bam



samtools view -H data_preproc/umi_dedup/zt21_2_ribo.sort.bam > temp_sam_header.txt
samtools view data_preproc/umi_dedup/zt21_2_ribo.sort.bam | \
  awk '32 <= length($10)' | \
  awk 'length($10) <= 35' | \
  cat temp_sam_header.txt - | \
  samtools view -bS - > data_preproc/bam_filtered_len/zt21_2_ribo.sort.bam
rm temp_sam_header.txt
samtools index data_preproc/bam_filtered_len/zt21_2_ribo.sort.bam



