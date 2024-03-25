#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/fastq_cl_umi

echo '### Processing "zt0_1_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt0_1_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt0_1_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt0_1_ribo.fastq.gz


echo '### Processing "zt0_2_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt0_2_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt0_2_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt0_2_ribo.fastq.gz


echo '### Processing "zt3_1_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt3_1_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt3_1_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt3_1_ribo.fastq.gz


echo '### Processing "zt3_2_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt3_2_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt3_2_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt3_2_ribo.fastq.gz


echo '### Processing "zt6_1_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt6_1_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt6_1_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt6_1_ribo.fastq.gz


echo '### Processing "zt6_2_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt6_2_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt6_2_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt6_2_ribo.fastq.gz


echo '### Processing "zt12_1_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt12_1_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt12_1_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt12_1_ribo.fastq.gz


echo '### Processing "zt12_2_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt12_2_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt12_2_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt12_2_ribo.fastq.gz


echo '### Processing "zt18_1_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt18_1_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt18_1_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt18_1_ribo.fastq.gz


echo '### Processing "zt18_2_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt18_2_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt18_2_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt18_2_ribo.fastq.gz


echo '### Processing "zt21_1_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt21_1_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt21_1_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt21_1_ribo.fastq.gz


echo '### Processing "zt21_2_ribo" ###'
umi_tools extract \
  --extract-method=string \
  -p NNNNN \
  --3prime \
  -L data_preproc/fastq_cl_umi/zt21_2_ribo.log \
  -I data_preproc/fastq_cl_adaptor/zt21_2_ribo.fastq.gz \
  -S data_preproc/fastq_cl_umi/zt21_2_ribo.fastq.gz


