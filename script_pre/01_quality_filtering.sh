#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/fastq_qf

echo '### Processing "zt0_1_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt0_1_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt0_1_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt0_1_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt0_1_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt0_1_ribo_fastp.json


echo '### Processing "zt0_2_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt0_2_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt0_2_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt0_2_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt0_2_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt0_2_ribo_fastp.json


echo '### Processing "zt3_1_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt3_1_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt3_1_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt3_1_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt3_1_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt3_1_ribo_fastp.json


echo '### Processing "zt3_2_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt3_2_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt3_2_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt3_2_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt3_2_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt3_2_ribo_fastp.json


echo '### Processing "zt6_1_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt6_1_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt6_1_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt6_1_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt6_1_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt6_1_ribo_fastp.json


echo '### Processing "zt6_2_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt6_2_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt6_2_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt6_2_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt6_2_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt6_2_ribo_fastp.json


echo '### Processing "zt12_1_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt12_1_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt12_1_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt12_1_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt12_1_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt12_1_ribo_fastp.json


echo '### Processing "zt12_2_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt12_2_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt12_2_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt12_2_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt12_2_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt12_2_ribo_fastp.json


echo '### Processing "zt18_1_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt18_1_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt18_1_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt18_1_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt18_1_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt18_1_ribo_fastp.json


echo '### Processing "zt18_2_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt18_2_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt18_2_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt18_2_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt18_2_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt18_2_ribo_fastp.json


echo '### Processing "zt21_1_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt21_1_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt21_1_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt21_1_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt21_1_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt21_1_ribo_fastp.json


echo '### Processing "zt21_2_ribo" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt21_2_ribo.fastq.gz \
  --out1 data_preproc/fastq_qf/zt21_2_ribo.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt21_2_ribo_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt21_2_ribo_fastp.html
mv fastp.json data_preproc/fastq_qf/zt21_2_ribo_fastp.json


echo '### Processing "zt0_1_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt0_1_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt0_1_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt0_1_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt0_1_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt0_1_rna_fastp.json


echo '### Processing "zt0_2_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt0_2_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt0_2_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt0_2_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt0_2_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt0_2_rna_fastp.json


echo '### Processing "zt3_1_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt3_1_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt3_1_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt3_1_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt3_1_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt3_1_rna_fastp.json


echo '### Processing "zt3_2_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt3_2_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt3_2_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt3_2_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt3_2_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt3_2_rna_fastp.json


echo '### Processing "zt6_1_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt6_1_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt6_1_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt6_1_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt6_1_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt6_1_rna_fastp.json


echo '### Processing "zt6_2_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt6_2_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt6_2_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt6_2_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt6_2_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt6_2_rna_fastp.json


echo '### Processing "zt12_1_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt12_1_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt12_1_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt12_1_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt12_1_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt12_1_rna_fastp.json


echo '### Processing "zt12_2_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt12_2_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt12_2_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt12_2_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt12_2_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt12_2_rna_fastp.json


echo '### Processing "zt18_1_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt18_1_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt18_1_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt18_1_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt18_1_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt18_1_rna_fastp.json


echo '### Processing "zt18_2_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt18_2_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt18_2_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt18_2_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt18_2_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt18_2_rna_fastp.json


echo '### Processing "zt21_1_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt21_1_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt21_1_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt21_1_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt21_1_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt21_1_rna_fastp.json


echo '### Processing "zt21_2_rna" ###'
fastp \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 15 \
  --disable_length_filtering \
  --disable_adapter_trimming \
  --disable_trim_poly_g \
  --thread 8 \
  --in1 data_preproc/fastq_cat/zt21_2_rna.fastq.gz \
  --out1 data_preproc/fastq_qf/zt21_2_rna.fastq.gz \
  --failed_out data_preproc/fastq_qf/zt21_2_rna_fail.fastq.gz
mv fastp.html data_preproc/fastq_qf/zt21_2_rna_fastp.html
mv fastp.json data_preproc/fastq_qf/zt21_2_rna_fastp.json


