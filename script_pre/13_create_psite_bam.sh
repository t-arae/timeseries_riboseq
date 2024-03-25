#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/bam_psite

echo '### Processing "zt0_1_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt0_1_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt0_1_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt0_1_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt0_1_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt0_1_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt0_1_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt0_2_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt0_2_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt0_2_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt0_2_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt0_2_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt0_2_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt0_2_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt12_1_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt12_1_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt12_1_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt12_1_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt12_1_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt12_1_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt12_1_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt12_2_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt12_2_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt12_2_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt12_2_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt12_2_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt12_2_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt12_2_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt18_1_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt18_1_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt18_1_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt18_1_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt18_1_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt18_1_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt18_1_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt18_2_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt18_2_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt18_2_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt18_2_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt18_2_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt18_2_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt18_2_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt21_1_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt21_1_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt21_1_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt21_1_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt21_1_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt21_1_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt21_1_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt21_2_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt21_2_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt21_2_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt21_2_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt21_2_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt21_2_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt21_2_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt3_1_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt3_1_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt3_1_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt3_1_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt3_1_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt3_1_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt3_1_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt3_2_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt3_2_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt3_2_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt3_2_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt3_2_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt3_2_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt3_2_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt6_1_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt6_1_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt6_1_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt6_1_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt6_1_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt6_1_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt6_1_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

echo '### Processing "zt6_2_ribo" ###'
bedtools bamtobed -bed12 -i data_preproc/bam_filtered_len/zt6_2_ribo.sort.bam > temp.bed12

awk '($11==32){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -18 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt6_2_ribo_32.bed
awk '($11==33){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -13 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt6_2_ribo_33.bed
awk '($11==34){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -19 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt6_2_ribo_34.bed
awk '($11==35){print}' temp.bed12 | \
  bedtools bed12tobed6 -i stdin | \
  bedtools slop -s -l -14 -r -20 -g idx/idx_star/chrNameLength.txt -i - > data_preproc/bam_psite/zt6_2_ribo_35.bed
cat data_preproc/bam_psite/*.bed | \
  sort -k1,1 -k2,2n | \
  bedtools bedtobam -g idx/idx_star/chrNameLength.txt > data_preproc/bam_psite/zt6_2_ribo.sort.bam
rm data_preproc/bam_psite/*.bed
rm temp.bed12

samtools index data_preproc/bam_psite/zt0_1_ribo.sort.bam
samtools index data_preproc/bam_psite/zt0_2_ribo.sort.bam
samtools index data_preproc/bam_psite/zt12_1_ribo.sort.bam
samtools index data_preproc/bam_psite/zt12_2_ribo.sort.bam
samtools index data_preproc/bam_psite/zt18_1_ribo.sort.bam
samtools index data_preproc/bam_psite/zt18_2_ribo.sort.bam
samtools index data_preproc/bam_psite/zt21_1_ribo.sort.bam
samtools index data_preproc/bam_psite/zt21_2_ribo.sort.bam
samtools index data_preproc/bam_psite/zt3_1_ribo.sort.bam
samtools index data_preproc/bam_psite/zt3_2_ribo.sort.bam
samtools index data_preproc/bam_psite/zt6_1_ribo.sort.bam
samtools index data_preproc/bam_psite/zt6_2_ribo.sort.bam

