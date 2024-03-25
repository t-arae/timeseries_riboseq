#!/bin/bash
set -eu
set -o pipefail

mkdir -p idx/idx_bt2_ncRNA/idx_bt2_ncRNA

bowtie2-build --threads 8 -f misc/fasta/Araport11_GFF3_genes_transposons.201606_ncRNA_sort.fasta idx/idx_bt2_ncRNA/idx_bt2_ncRNA

mkdir -p data_preproc/fastq_rm_ncRNA

echo '### Processing "zt0_1_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt0_1_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt0_1_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt0_1_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt0_1_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt0_1_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt0_1_ribo_idxstats.txt


echo '### Processing "zt0_2_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt0_2_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt0_2_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt0_2_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt0_2_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt0_2_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt0_2_ribo_idxstats.txt


echo '### Processing "zt3_1_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt3_1_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt3_1_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt3_1_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt3_1_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt3_1_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt3_1_ribo_idxstats.txt


echo '### Processing "zt3_2_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt3_2_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt3_2_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt3_2_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt3_2_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt3_2_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt3_2_ribo_idxstats.txt


echo '### Processing "zt6_1_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt6_1_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt6_1_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt6_1_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt6_1_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt6_1_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt6_1_ribo_idxstats.txt


echo '### Processing "zt6_2_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt6_2_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt6_2_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt6_2_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt6_2_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt6_2_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt6_2_ribo_idxstats.txt


echo '### Processing "zt12_1_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt12_1_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt12_1_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt12_1_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt12_1_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt12_1_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt12_1_ribo_idxstats.txt


echo '### Processing "zt12_2_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt12_2_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt12_2_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt12_2_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt12_2_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt12_2_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt12_2_ribo_idxstats.txt


echo '### Processing "zt18_1_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt18_1_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt18_1_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt18_1_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt18_1_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt18_1_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt18_1_ribo_idxstats.txt


echo '### Processing "zt18_2_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt18_2_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt18_2_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt18_2_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt18_2_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt18_2_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt18_2_ribo_idxstats.txt


echo '### Processing "zt21_1_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt21_1_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt21_1_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt21_1_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt21_1_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt21_1_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt21_1_ribo_idxstats.txt


echo '### Processing "zt21_2_ribo" ###'
gzip -dc data_preproc/fastq_cl_umi/zt21_2_ribo.fastq.gz | \
  bowtie2 -L 20 -p 8 -t --quiet -x idx/idx_bt2_ncRNA/idx_bt2_ncRNA -U - | \
  samtools view -uS > temp.bam
samtools view -@ 8 -u -f 0x4 temp.bam | \
  samtools fastq - | \
  gzip > data_preproc/fastq_rm_ncRNA/zt21_2_ribo.fastq.gz
samtools view -@ 8 -u -F 0x4 temp.bam | samtools sort > data_preproc/fastq_rm_ncRNA/zt21_2_ribo_al.sort.bam
rm temp.bam
samtools index -@ 8 data_preproc/fastq_rm_ncRNA/zt21_2_ribo_al.sort.bam
samtools idxstats -@ 8 data_preproc/fastq_rm_ncRNA/zt21_2_ribo_al.sort.bam > data_preproc/fastq_rm_ncRNA/zt21_2_ribo_idxstats.txt


