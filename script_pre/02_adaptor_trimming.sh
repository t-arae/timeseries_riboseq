#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/fastq_cl_adaptor

echo '### Processing "zt0_1_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt0_1_ribo.fastq.gz \
  data_preproc/fastq_qf/zt0_1_ribo.fastq.gz


echo '### Processing "zt0_2_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt0_2_ribo.fastq.gz \
  data_preproc/fastq_qf/zt0_2_ribo.fastq.gz


echo '### Processing "zt3_1_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt3_1_ribo.fastq.gz \
  data_preproc/fastq_qf/zt3_1_ribo.fastq.gz


echo '### Processing "zt3_2_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt3_2_ribo.fastq.gz \
  data_preproc/fastq_qf/zt3_2_ribo.fastq.gz


echo '### Processing "zt6_1_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt6_1_ribo.fastq.gz \
  data_preproc/fastq_qf/zt6_1_ribo.fastq.gz


echo '### Processing "zt6_2_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt6_2_ribo.fastq.gz \
  data_preproc/fastq_qf/zt6_2_ribo.fastq.gz


echo '### Processing "zt12_1_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt12_1_ribo.fastq.gz \
  data_preproc/fastq_qf/zt12_1_ribo.fastq.gz


echo '### Processing "zt12_2_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt12_2_ribo.fastq.gz \
  data_preproc/fastq_qf/zt12_2_ribo.fastq.gz


echo '### Processing "zt18_1_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt18_1_ribo.fastq.gz \
  data_preproc/fastq_qf/zt18_1_ribo.fastq.gz


echo '### Processing "zt18_2_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt18_2_ribo.fastq.gz \
  data_preproc/fastq_qf/zt18_2_ribo.fastq.gz


echo '### Processing "zt21_1_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt21_1_ribo.fastq.gz \
  data_preproc/fastq_qf/zt21_1_ribo.fastq.gz


echo '### Processing "zt21_2_ribo" ###'
cutadapt \
  --trimmed-only \
  --minimum-length 15 \
  --cores 8 \
  -a AGCTAAGATCGGAAGAGCACACGTCTGAA \
  -o data_preproc/fastq_cl_adaptor/zt21_2_ribo.fastq.gz \
  data_preproc/fastq_qf/zt21_2_ribo.fastq.gz


