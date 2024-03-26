#!/bin/bash
set -eu
set -o pipefail

mkdir -p report_fastqc/fastq_dra


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt0_1_ribo.fastq.gz


echo '### Processing "zt0_1_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt0_1_rna.fastq.gz


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt0_2_ribo.fastq.gz


echo '### Processing "zt0_2_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt0_2_rna.fastq.gz


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt12_1_ribo.fastq.gz


echo '### Processing "zt12_1_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt12_1_rna.fastq.gz


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt12_2_ribo.fastq.gz


echo '### Processing "zt12_2_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt12_2_rna.fastq.gz


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt18_1_ribo.fastq.gz


echo '### Processing "zt18_1_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt18_1_rna.fastq.gz


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt18_2_ribo.fastq.gz


echo '### Processing "zt18_2_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt18_2_rna.fastq.gz


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt21_1_ribo.fastq.gz


echo '### Processing "zt21_1_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt21_1_rna.fastq.gz


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt21_2_ribo.fastq.gz


echo '### Processing "zt21_2_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt21_2_rna.fastq.gz


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt3_1_ribo.fastq.gz


echo '### Processing "zt3_1_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt3_1_rna.fastq.gz


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt3_2_ribo.fastq.gz


echo '### Processing "zt3_2_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt3_2_rna.fastq.gz


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt6_1_ribo.fastq.gz


echo '### Processing "zt6_1_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt6_1_rna.fastq.gz


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt6_2_ribo.fastq.gz


echo '### Processing "zt6_2_rna" ###'
fastqc -o report_fastqc/fastq_dra -t 16 data_source/fastq_dra/zt6_2_rna.fastq.gz


mkdir -p report_fastqc/fastq_qf


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt0_1_ribo.fastq.gz


echo '### Processing "zt0_1_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt0_1_ribo_fail.fastq.gz


echo '### Processing "zt0_1_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt0_1_rna.fastq.gz


echo '### Processing "zt0_1_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt0_1_rna_fail.fastq.gz


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt0_2_ribo.fastq.gz


echo '### Processing "zt0_2_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt0_2_ribo_fail.fastq.gz


echo '### Processing "zt0_2_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt0_2_rna.fastq.gz


echo '### Processing "zt0_2_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt0_2_rna_fail.fastq.gz


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt12_1_ribo.fastq.gz


echo '### Processing "zt12_1_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt12_1_ribo_fail.fastq.gz


echo '### Processing "zt12_1_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt12_1_rna.fastq.gz


echo '### Processing "zt12_1_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt12_1_rna_fail.fastq.gz


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt12_2_ribo.fastq.gz


echo '### Processing "zt12_2_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt12_2_ribo_fail.fastq.gz


echo '### Processing "zt12_2_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt12_2_rna.fastq.gz


echo '### Processing "zt12_2_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt12_2_rna_fail.fastq.gz


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt18_1_ribo.fastq.gz


echo '### Processing "zt18_1_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt18_1_ribo_fail.fastq.gz


echo '### Processing "zt18_1_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt18_1_rna.fastq.gz


echo '### Processing "zt18_1_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt18_1_rna_fail.fastq.gz


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt18_2_ribo.fastq.gz


echo '### Processing "zt18_2_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt18_2_ribo_fail.fastq.gz


echo '### Processing "zt18_2_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt18_2_rna.fastq.gz


echo '### Processing "zt18_2_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt18_2_rna_fail.fastq.gz


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt21_1_ribo.fastq.gz


echo '### Processing "zt21_1_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt21_1_ribo_fail.fastq.gz


echo '### Processing "zt21_1_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt21_1_rna.fastq.gz


echo '### Processing "zt21_1_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt21_1_rna_fail.fastq.gz


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt21_2_ribo.fastq.gz


echo '### Processing "zt21_2_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt21_2_ribo_fail.fastq.gz


echo '### Processing "zt21_2_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt21_2_rna.fastq.gz


echo '### Processing "zt21_2_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt21_2_rna_fail.fastq.gz


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt3_1_ribo.fastq.gz


echo '### Processing "zt3_1_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt3_1_ribo_fail.fastq.gz


echo '### Processing "zt3_1_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt3_1_rna.fastq.gz


echo '### Processing "zt3_1_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt3_1_rna_fail.fastq.gz


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt3_2_ribo.fastq.gz


echo '### Processing "zt3_2_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt3_2_ribo_fail.fastq.gz


echo '### Processing "zt3_2_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt3_2_rna.fastq.gz


echo '### Processing "zt3_2_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt3_2_rna_fail.fastq.gz


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt6_1_ribo.fastq.gz


echo '### Processing "zt6_1_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt6_1_ribo_fail.fastq.gz


echo '### Processing "zt6_1_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt6_1_rna.fastq.gz


echo '### Processing "zt6_1_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt6_1_rna_fail.fastq.gz


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt6_2_ribo.fastq.gz


echo '### Processing "zt6_2_ribo_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt6_2_ribo_fail.fastq.gz


echo '### Processing "zt6_2_rna" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt6_2_rna.fastq.gz


echo '### Processing "zt6_2_rna_fail" ###'
fastqc -o report_fastqc/fastq_qf -t 16 data_preproc/fastq_qf/zt6_2_rna_fail.fastq.gz


mkdir -p report_fastqc/fastq_cl_adaptor


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt0_1_ribo.fastq.gz


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt0_2_ribo.fastq.gz


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt12_1_ribo.fastq.gz


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt12_2_ribo.fastq.gz


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt18_1_ribo.fastq.gz


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt18_2_ribo.fastq.gz


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt21_1_ribo.fastq.gz


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt21_2_ribo.fastq.gz


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt3_1_ribo.fastq.gz


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt3_2_ribo.fastq.gz


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt6_1_ribo.fastq.gz


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_adaptor -t 16 data_preproc/fastq_cl_adaptor/zt6_2_ribo.fastq.gz


mkdir -p report_fastqc/fastq_cl_umi


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt0_1_ribo.fastq.gz


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt0_2_ribo.fastq.gz


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt12_1_ribo.fastq.gz


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt12_2_ribo.fastq.gz


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt18_1_ribo.fastq.gz


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt18_2_ribo.fastq.gz


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt21_1_ribo.fastq.gz


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt21_2_ribo.fastq.gz


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt3_1_ribo.fastq.gz


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt3_2_ribo.fastq.gz


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt6_1_ribo.fastq.gz


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/fastq_cl_umi -t 16 data_preproc/fastq_cl_umi/zt6_2_ribo.fastq.gz


mkdir -p report_fastqc/fastq_rm_ncRNA


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt0_1_ribo.fastq.gz


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt0_2_ribo.fastq.gz


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt12_1_ribo.fastq.gz


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt12_2_ribo.fastq.gz


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt18_1_ribo.fastq.gz


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt18_2_ribo.fastq.gz


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt21_1_ribo.fastq.gz


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt21_2_ribo.fastq.gz


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt3_1_ribo.fastq.gz


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt3_2_ribo.fastq.gz


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt6_1_ribo.fastq.gz


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_ncRNA -t 16 data_preproc/fastq_rm_ncRNA/zt6_2_ribo.fastq.gz


mkdir -p report_fastqc/fastq_rm_marker


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt0_1_ribo.fastq.gz


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt0_2_ribo.fastq.gz


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt12_1_ribo.fastq.gz


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt12_2_ribo.fastq.gz


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt18_1_ribo.fastq.gz


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt18_2_ribo.fastq.gz


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt21_1_ribo.fastq.gz


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt21_2_ribo.fastq.gz


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt3_1_ribo.fastq.gz


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt3_2_ribo.fastq.gz


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt6_1_ribo.fastq.gz


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/fastq_rm_marker -t 16 data_preproc/fastq_rm_marker/zt6_2_ribo.fastq.gz


mkdir -p report_fastqc/mapped_by_star


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt0_1_ribo.sort.bam


echo '### Processing "zt0_1_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt0_1_rna.sort.bam


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt0_2_ribo.sort.bam


echo '### Processing "zt0_2_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt0_2_rna.sort.bam


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt12_1_ribo.sort.bam


echo '### Processing "zt12_1_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt12_1_rna.sort.bam


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt12_2_ribo.sort.bam


echo '### Processing "zt12_2_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt12_2_rna.sort.bam


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt18_1_ribo.sort.bam


echo '### Processing "zt18_1_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt18_1_rna.sort.bam


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt18_2_ribo.sort.bam


echo '### Processing "zt18_2_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt18_2_rna.sort.bam


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt21_1_ribo.sort.bam


echo '### Processing "zt21_1_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt21_1_rna.sort.bam


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt21_2_ribo.sort.bam


echo '### Processing "zt21_2_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt21_2_rna.sort.bam


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt3_1_ribo.sort.bam


echo '### Processing "zt3_1_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt3_1_rna.sort.bam


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt3_2_ribo.sort.bam


echo '### Processing "zt3_2_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt3_2_rna.sort.bam


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt6_1_ribo.sort.bam


echo '### Processing "zt6_1_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt6_1_rna.sort.bam


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt6_2_ribo.sort.bam


echo '### Processing "zt6_2_rna" ###'
fastqc -o report_fastqc/mapped_by_star -t 16 data_preproc/mapped_by_star/zt6_2_rna.sort.bam


mkdir -p report_fastqc/umi_dedup


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt0_1_ribo.sort.bam


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt0_2_ribo.sort.bam


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt12_1_ribo.sort.bam


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt12_2_ribo.sort.bam


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt18_1_ribo.sort.bam


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt18_2_ribo.sort.bam


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt21_1_ribo.sort.bam


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt21_2_ribo.sort.bam


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt3_1_ribo.sort.bam


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt3_2_ribo.sort.bam


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt6_1_ribo.sort.bam


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/umi_dedup -t 16 data_preproc/umi_dedup/zt6_2_ribo.sort.bam


mkdir -p report_fastqc/bam_filtered_len


echo '### Processing "zt0_1_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt0_1_ribo.sort.bam


echo '### Processing "zt0_2_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt0_2_ribo.sort.bam


echo '### Processing "zt12_1_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt12_1_ribo.sort.bam


echo '### Processing "zt12_2_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt12_2_ribo.sort.bam


echo '### Processing "zt18_1_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt18_1_ribo.sort.bam


echo '### Processing "zt18_2_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt18_2_ribo.sort.bam


echo '### Processing "zt21_1_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt21_1_ribo.sort.bam


echo '### Processing "zt21_2_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt21_2_ribo.sort.bam


echo '### Processing "zt3_1_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt3_1_ribo.sort.bam


echo '### Processing "zt3_2_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt3_2_ribo.sort.bam


echo '### Processing "zt6_1_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt6_1_ribo.sort.bam


echo '### Processing "zt6_2_ribo" ###'
fastqc -o report_fastqc/bam_filtered_len -t 16 data_preproc/bam_filtered_len/zt6_2_ribo.sort.bam

