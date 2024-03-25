#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/bam_transcriptome

java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt0_1_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt0_1_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt0_2_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt0_2_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt3_1_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt3_1_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt3_2_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt3_2_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt6_1_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt6_1_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt6_2_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt6_2_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt12_1_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt12_1_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt12_2_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt12_2_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt18_1_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt18_1_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt18_2_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt18_2_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt21_1_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt21_1_ribo.tr.bam


java \
  -Xmx512M \
  -jar ubu-1.2-jar-with-dependencies.zip \
  sam-xlate \
  --single \
  --reverse \
  --bed data_modified/gff_gtf/transcript.bed12 \
  --in data_preproc/umi_dedup/zt21_2_ribo.sort.bam \
  --out data_preproc/bam_transcriptome/zt21_2_ribo.tr.bam



samtools sort -@ 8 data_preproc/bam_transcriptome/zt0_1_ribo.tr.bam > data_preproc/bam_transcriptome/zt0_1_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt0_1_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt0_1_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt0_2_ribo.tr.bam > data_preproc/bam_transcriptome/zt0_2_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt0_2_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt0_2_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt3_1_ribo.tr.bam > data_preproc/bam_transcriptome/zt3_1_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt3_1_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt3_1_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt3_2_ribo.tr.bam > data_preproc/bam_transcriptome/zt3_2_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt3_2_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt3_2_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt6_1_ribo.tr.bam > data_preproc/bam_transcriptome/zt6_1_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt6_1_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt6_1_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt6_2_ribo.tr.bam > data_preproc/bam_transcriptome/zt6_2_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt6_2_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt6_2_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt12_1_ribo.tr.bam > data_preproc/bam_transcriptome/zt12_1_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt12_1_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt12_1_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt12_2_ribo.tr.bam > data_preproc/bam_transcriptome/zt12_2_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt12_2_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt12_2_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt18_1_ribo.tr.bam > data_preproc/bam_transcriptome/zt18_1_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt18_1_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt18_1_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt18_2_ribo.tr.bam > data_preproc/bam_transcriptome/zt18_2_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt18_2_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt18_2_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt21_1_ribo.tr.bam > data_preproc/bam_transcriptome/zt21_1_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt21_1_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt21_1_ribo.tr.bam


samtools sort -@ 8 data_preproc/bam_transcriptome/zt21_2_ribo.tr.bam > data_preproc/bam_transcriptome/zt21_2_ribo.tr.sort.bam
samtools index data_preproc/bam_transcriptome/zt21_2_ribo.tr.sort.bam
rm data_preproc/bam_transcriptome/zt21_2_ribo.tr.bam


