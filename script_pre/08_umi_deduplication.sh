#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/umi_dedup

umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt0_1_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt0_1_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt0_1_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt0_1_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt0_2_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt0_2_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt0_2_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt0_2_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt3_1_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt3_1_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt3_1_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt3_1_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt3_2_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt3_2_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt3_2_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt3_2_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt6_1_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt6_1_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt6_1_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt6_1_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt6_2_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt6_2_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt6_2_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt6_2_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt12_1_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt12_1_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt12_1_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt12_1_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt12_2_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt12_2_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt12_2_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt12_2_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt18_1_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt18_1_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt18_1_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt18_1_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt18_2_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt18_2_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt18_2_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt18_2_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt21_1_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt21_1_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt21_1_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt21_1_ribo \
  --method=unique \
  --random-seed=777


umi_tools dedup \
  --stdin data_preproc/rm_softclip/zt21_2_ribo.sort.bam \
  --stdout data_preproc/umi_dedup/zt21_2_ribo.sort.bam \
  --log data_preproc/umi_dedup/zt21_2_ribo.log \
  --output-stats=data_preproc/umi_dedup/zt21_2_ribo \
  --method=unique \
  --random-seed=777


