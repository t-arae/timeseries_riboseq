#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_source/fastq_dra

echo '### Processing "zt0_1_ribo" ###'
curl -C - -o zt0_1_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479933/DRR495762.fastq.bz2
bzip2 -dc zt0_1_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt0_1_ribo.fastq.gz
rm zt0_1_ribo.fastq.bz2


echo '### Processing "zt0_2_ribo" ###'
curl -C - -o zt0_2_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479934/DRR495763.fastq.bz2
bzip2 -dc zt0_2_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt0_2_ribo.fastq.gz
rm zt0_2_ribo.fastq.bz2


echo '### Processing "zt3_1_ribo" ###'
curl -C - -o zt3_1_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479935/DRR495764.fastq.bz2
bzip2 -dc zt3_1_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt3_1_ribo.fastq.gz
rm zt3_1_ribo.fastq.bz2


echo '### Processing "zt3_2_ribo" ###'
curl -C - -o zt3_2_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479936/DRR495765.fastq.bz2
bzip2 -dc zt3_2_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt3_2_ribo.fastq.gz
rm zt3_2_ribo.fastq.bz2


echo '### Processing "zt6_1_ribo" ###'
curl -C - -o zt6_1_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479937/DRR495766.fastq.bz2
bzip2 -dc zt6_1_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt6_1_ribo.fastq.gz
rm zt6_1_ribo.fastq.bz2


echo '### Processing "zt6_2_ribo" ###'
curl -C - -o zt6_2_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479938/DRR495767.fastq.bz2
bzip2 -dc zt6_2_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt6_2_ribo.fastq.gz
rm zt6_2_ribo.fastq.bz2


echo '### Processing "zt12_1_ribo" ###'
curl -C - -o zt12_1_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479939/DRR495768.fastq.bz2
bzip2 -dc zt12_1_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt12_1_ribo.fastq.gz
rm zt12_1_ribo.fastq.bz2


echo '### Processing "zt12_2_ribo" ###'
curl -C - -o zt12_2_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479940/DRR495769.fastq.bz2
bzip2 -dc zt12_2_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt12_2_ribo.fastq.gz
rm zt12_2_ribo.fastq.bz2


echo '### Processing "zt18_1_ribo" ###'
curl -C - -o zt18_1_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479941/DRR495770.fastq.bz2
bzip2 -dc zt18_1_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt18_1_ribo.fastq.gz
rm zt18_1_ribo.fastq.bz2


echo '### Processing "zt18_2_ribo" ###'
curl -C - -o zt18_2_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479942/DRR495771.fastq.bz2
bzip2 -dc zt18_2_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt18_2_ribo.fastq.gz
rm zt18_2_ribo.fastq.bz2


echo '### Processing "zt21_1_ribo" ###'
curl -C - -o zt21_1_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479943/DRR495772.fastq.bz2
bzip2 -dc zt21_1_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt21_1_ribo.fastq.gz
rm zt21_1_ribo.fastq.bz2


echo '### Processing "zt21_2_ribo" ###'
curl -C - -o zt21_2_ribo.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479944/DRR495773.fastq.bz2
bzip2 -dc zt21_2_ribo.fastq.bz2 | gzip > data_source/fastq_dra/zt21_2_ribo.fastq.gz
rm zt21_2_ribo.fastq.bz2


echo '### Processing "zt0_1_rna" ###'
curl -C - -o zt0_1_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479945/DRR495774.fastq.bz2
bzip2 -dc zt0_1_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt0_1_rna.fastq.gz
rm zt0_1_rna.fastq.bz2


echo '### Processing "zt0_2_rna" ###'
curl -C - -o zt0_2_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479946/DRR495775.fastq.bz2
bzip2 -dc zt0_2_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt0_2_rna.fastq.gz
rm zt0_2_rna.fastq.bz2


echo '### Processing "zt3_1_rna" ###'
curl -C - -o zt3_1_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479947/DRR495776.fastq.bz2
bzip2 -dc zt3_1_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt3_1_rna.fastq.gz
rm zt3_1_rna.fastq.bz2


echo '### Processing "zt3_2_rna" ###'
curl -C - -o zt3_2_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479948/DRR495777.fastq.bz2
bzip2 -dc zt3_2_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt3_2_rna.fastq.gz
rm zt3_2_rna.fastq.bz2


echo '### Processing "zt6_1_rna" ###'
curl -C - -o zt6_1_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479949/DRR495778.fastq.bz2
bzip2 -dc zt6_1_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt6_1_rna.fastq.gz
rm zt6_1_rna.fastq.bz2


echo '### Processing "zt6_2_rna" ###'
curl -C - -o zt6_2_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479950/DRR495779.fastq.bz2
bzip2 -dc zt6_2_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt6_2_rna.fastq.gz
rm zt6_2_rna.fastq.bz2


echo '### Processing "zt12_1_rna" ###'
curl -C - -o zt12_1_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479951/DRR495780.fastq.bz2
bzip2 -dc zt12_1_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt12_1_rna.fastq.gz
rm zt12_1_rna.fastq.bz2


echo '### Processing "zt12_2_rna" ###'
curl -C - -o zt12_2_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479952/DRR495781.fastq.bz2
bzip2 -dc zt12_2_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt12_2_rna.fastq.gz
rm zt12_2_rna.fastq.bz2


echo '### Processing "zt18_1_rna" ###'
curl -C - -o zt18_1_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479953/DRR495782.fastq.bz2
bzip2 -dc zt18_1_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt18_1_rna.fastq.gz
rm zt18_1_rna.fastq.bz2


echo '### Processing "zt18_2_rna" ###'
curl -C - -o zt18_2_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479954/DRR495783.fastq.bz2
bzip2 -dc zt18_2_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt18_2_rna.fastq.gz
rm zt18_2_rna.fastq.bz2


echo '### Processing "zt21_1_rna" ###'
curl -C - -o zt21_1_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479955/DRR495784.fastq.bz2
bzip2 -dc zt21_1_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt21_1_rna.fastq.gz
rm zt21_1_rna.fastq.bz2


echo '### Processing "zt21_2_rna" ###'
curl -C - -o zt21_2_rna.fastq.bz2 https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016840/DRX479956/DRR495785.fastq.bz2
bzip2 -dc zt21_2_rna.fastq.bz2 | gzip > data_source/fastq_dra/zt21_2_rna.fastq.gz
rm zt21_2_rna.fastq.bz2


