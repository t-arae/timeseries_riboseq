#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/bam_coverage/bam_psite
echo '### Processing "zt0_1_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt0_1_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt0_1_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt0_2_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt0_2_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt0_2_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt12_1_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt12_1_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt12_1_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt12_2_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt12_2_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt12_2_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt18_1_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt18_1_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt18_1_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt18_2_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt18_2_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt18_2_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt21_1_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt21_1_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt21_1_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt21_2_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt21_2_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt21_2_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt3_1_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt3_1_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt3_1_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt3_2_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt3_2_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt3_2_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt6_1_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt6_1_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt6_1_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

echo '### Processing "zt6_2_ribo" ###'
bamCoverage \
  -b data_preproc/bam_psite/zt6_2_ribo.sort.bam \
  -o data_preproc/bam_coverage/bam_psite/zt6_2_ribo_none.bw \
  --binSize 1 \
  --numberOfProcessors 8 \
  --normalizeUsing None \
  --verbose

