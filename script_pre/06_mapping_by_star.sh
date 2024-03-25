#!/bin/bash
set -eu
set -o pipefail

mkdir -p idx/idx_star

STAR \
  --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --genomeSAindexNbases 12 \
  --genomeFastaFiles misc/fasta/BSgenome.Athaliana.TAIR.TAIR9.fasta \
  --sjdbOverhang 41 \
  --sjdbGTFfile misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gff \
  --sjdbGTFtagExonParentTranscript Parent \
  --sjdbInsertSave All

mkdir -p idx/idx_star_rna

STAR \
  --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --genomeSAindexNbases 12 \
  --genomeFastaFiles misc/fasta/BSgenome.Athaliana.TAIR.TAIR9.fasta \
  --sjdbOverhang 49 \
  --sjdbGTFfile misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gff \
  --sjdbGTFtagExonParentTranscript Parent \
  --sjdbInsertSave All

mkdir -p data_preproc/mapped_by_star

echo '### Processing "zt0_1_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt0_1_ribo.fastq.gz) \
  --outFileNamePrefix zt0_1_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt0_1_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt0_1_ribo.sort.bam
  rm zt0_1_riboAligned.out.bam
  mv zt0_1_riboLog.final.out data_preproc/mapped_by_star/zt0_1_ribo.final.log
  mv zt0_1_riboLog.out data_preproc/mapped_by_star/zt0_1_ribo.log
  mv zt0_1_riboSJ.out.tab data_preproc/mapped_by_star/zt0_1_ribo_sj.tsv
  rm -rf zt0_1_ribo_STARgenome
  rm zt0_1_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt0_1_ribo.sort.bam
fi


echo '### Processing "zt0_2_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt0_2_ribo.fastq.gz) \
  --outFileNamePrefix zt0_2_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt0_2_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt0_2_ribo.sort.bam
  rm zt0_2_riboAligned.out.bam
  mv zt0_2_riboLog.final.out data_preproc/mapped_by_star/zt0_2_ribo.final.log
  mv zt0_2_riboLog.out data_preproc/mapped_by_star/zt0_2_ribo.log
  mv zt0_2_riboSJ.out.tab data_preproc/mapped_by_star/zt0_2_ribo_sj.tsv
  rm -rf zt0_2_ribo_STARgenome
  rm zt0_2_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt0_2_ribo.sort.bam
fi


echo '### Processing "zt3_1_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt3_1_ribo.fastq.gz) \
  --outFileNamePrefix zt3_1_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt3_1_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt3_1_ribo.sort.bam
  rm zt3_1_riboAligned.out.bam
  mv zt3_1_riboLog.final.out data_preproc/mapped_by_star/zt3_1_ribo.final.log
  mv zt3_1_riboLog.out data_preproc/mapped_by_star/zt3_1_ribo.log
  mv zt3_1_riboSJ.out.tab data_preproc/mapped_by_star/zt3_1_ribo_sj.tsv
  rm -rf zt3_1_ribo_STARgenome
  rm zt3_1_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt3_1_ribo.sort.bam
fi


echo '### Processing "zt3_2_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt3_2_ribo.fastq.gz) \
  --outFileNamePrefix zt3_2_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt3_2_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt3_2_ribo.sort.bam
  rm zt3_2_riboAligned.out.bam
  mv zt3_2_riboLog.final.out data_preproc/mapped_by_star/zt3_2_ribo.final.log
  mv zt3_2_riboLog.out data_preproc/mapped_by_star/zt3_2_ribo.log
  mv zt3_2_riboSJ.out.tab data_preproc/mapped_by_star/zt3_2_ribo_sj.tsv
  rm -rf zt3_2_ribo_STARgenome
  rm zt3_2_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt3_2_ribo.sort.bam
fi


echo '### Processing "zt6_1_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt6_1_ribo.fastq.gz) \
  --outFileNamePrefix zt6_1_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt6_1_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt6_1_ribo.sort.bam
  rm zt6_1_riboAligned.out.bam
  mv zt6_1_riboLog.final.out data_preproc/mapped_by_star/zt6_1_ribo.final.log
  mv zt6_1_riboLog.out data_preproc/mapped_by_star/zt6_1_ribo.log
  mv zt6_1_riboSJ.out.tab data_preproc/mapped_by_star/zt6_1_ribo_sj.tsv
  rm -rf zt6_1_ribo_STARgenome
  rm zt6_1_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt6_1_ribo.sort.bam
fi


echo '### Processing "zt6_2_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt6_2_ribo.fastq.gz) \
  --outFileNamePrefix zt6_2_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt6_2_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt6_2_ribo.sort.bam
  rm zt6_2_riboAligned.out.bam
  mv zt6_2_riboLog.final.out data_preproc/mapped_by_star/zt6_2_ribo.final.log
  mv zt6_2_riboLog.out data_preproc/mapped_by_star/zt6_2_ribo.log
  mv zt6_2_riboSJ.out.tab data_preproc/mapped_by_star/zt6_2_ribo_sj.tsv
  rm -rf zt6_2_ribo_STARgenome
  rm zt6_2_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt6_2_ribo.sort.bam
fi


echo '### Processing "zt12_1_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt12_1_ribo.fastq.gz) \
  --outFileNamePrefix zt12_1_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt12_1_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt12_1_ribo.sort.bam
  rm zt12_1_riboAligned.out.bam
  mv zt12_1_riboLog.final.out data_preproc/mapped_by_star/zt12_1_ribo.final.log
  mv zt12_1_riboLog.out data_preproc/mapped_by_star/zt12_1_ribo.log
  mv zt12_1_riboSJ.out.tab data_preproc/mapped_by_star/zt12_1_ribo_sj.tsv
  rm -rf zt12_1_ribo_STARgenome
  rm zt12_1_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt12_1_ribo.sort.bam
fi


echo '### Processing "zt12_2_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt12_2_ribo.fastq.gz) \
  --outFileNamePrefix zt12_2_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt12_2_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt12_2_ribo.sort.bam
  rm zt12_2_riboAligned.out.bam
  mv zt12_2_riboLog.final.out data_preproc/mapped_by_star/zt12_2_ribo.final.log
  mv zt12_2_riboLog.out data_preproc/mapped_by_star/zt12_2_ribo.log
  mv zt12_2_riboSJ.out.tab data_preproc/mapped_by_star/zt12_2_ribo_sj.tsv
  rm -rf zt12_2_ribo_STARgenome
  rm zt12_2_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt12_2_ribo.sort.bam
fi


echo '### Processing "zt18_1_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt18_1_ribo.fastq.gz) \
  --outFileNamePrefix zt18_1_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt18_1_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt18_1_ribo.sort.bam
  rm zt18_1_riboAligned.out.bam
  mv zt18_1_riboLog.final.out data_preproc/mapped_by_star/zt18_1_ribo.final.log
  mv zt18_1_riboLog.out data_preproc/mapped_by_star/zt18_1_ribo.log
  mv zt18_1_riboSJ.out.tab data_preproc/mapped_by_star/zt18_1_ribo_sj.tsv
  rm -rf zt18_1_ribo_STARgenome
  rm zt18_1_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt18_1_ribo.sort.bam
fi


echo '### Processing "zt18_2_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt18_2_ribo.fastq.gz) \
  --outFileNamePrefix zt18_2_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt18_2_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt18_2_ribo.sort.bam
  rm zt18_2_riboAligned.out.bam
  mv zt18_2_riboLog.final.out data_preproc/mapped_by_star/zt18_2_ribo.final.log
  mv zt18_2_riboLog.out data_preproc/mapped_by_star/zt18_2_ribo.log
  mv zt18_2_riboSJ.out.tab data_preproc/mapped_by_star/zt18_2_ribo_sj.tsv
  rm -rf zt18_2_ribo_STARgenome
  rm zt18_2_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt18_2_ribo.sort.bam
fi


echo '### Processing "zt21_1_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt21_1_ribo.fastq.gz) \
  --outFileNamePrefix zt21_1_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt21_1_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt21_1_ribo.sort.bam
  rm zt21_1_riboAligned.out.bam
  mv zt21_1_riboLog.final.out data_preproc/mapped_by_star/zt21_1_ribo.final.log
  mv zt21_1_riboLog.out data_preproc/mapped_by_star/zt21_1_ribo.log
  mv zt21_1_riboSJ.out.tab data_preproc/mapped_by_star/zt21_1_ribo_sj.tsv
  rm -rf zt21_1_ribo_STARgenome
  rm zt21_1_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt21_1_ribo.sort.bam
fi


echo '### Processing "zt21_2_ribo" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_rm_marker/zt21_2_ribo.fastq.gz) \
  --outFileNamePrefix zt21_2_ribo

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt21_2_riboAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt21_2_ribo.sort.bam
  rm zt21_2_riboAligned.out.bam
  mv zt21_2_riboLog.final.out data_preproc/mapped_by_star/zt21_2_ribo.final.log
  mv zt21_2_riboLog.out data_preproc/mapped_by_star/zt21_2_ribo.log
  mv zt21_2_riboSJ.out.tab data_preproc/mapped_by_star/zt21_2_ribo_sj.tsv
  rm -rf zt21_2_ribo_STARgenome
  rm zt21_2_riboLog.progress.out
  samtools index data_preproc/mapped_by_star/zt21_2_ribo.sort.bam
fi

mkdir -p data_preproc/mapped_by_star

echo '### Processing "zt0_1_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt0_1_rna.fastq.gz) \
  --outFileNamePrefix zt0_1_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt0_1_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt0_1_rna.sort.bam
  rm zt0_1_rnaAligned.out.bam
  mv zt0_1_rnaLog.final.out data_preproc/mapped_by_star/zt0_1_rna.final.log
  mv zt0_1_rnaLog.out data_preproc/mapped_by_star/zt0_1_rna.log
  mv zt0_1_rnaSJ.out.tab data_preproc/mapped_by_star/zt0_1_rna_sj.tsv
  rm -rf zt0_1_rna_STARgenome
  rm zt0_1_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt0_1_rna.sort.bam
fi


echo '### Processing "zt0_2_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt0_2_rna.fastq.gz) \
  --outFileNamePrefix zt0_2_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt0_2_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt0_2_rna.sort.bam
  rm zt0_2_rnaAligned.out.bam
  mv zt0_2_rnaLog.final.out data_preproc/mapped_by_star/zt0_2_rna.final.log
  mv zt0_2_rnaLog.out data_preproc/mapped_by_star/zt0_2_rna.log
  mv zt0_2_rnaSJ.out.tab data_preproc/mapped_by_star/zt0_2_rna_sj.tsv
  rm -rf zt0_2_rna_STARgenome
  rm zt0_2_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt0_2_rna.sort.bam
fi


echo '### Processing "zt3_1_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt3_1_rna.fastq.gz) \
  --outFileNamePrefix zt3_1_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt3_1_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt3_1_rna.sort.bam
  rm zt3_1_rnaAligned.out.bam
  mv zt3_1_rnaLog.final.out data_preproc/mapped_by_star/zt3_1_rna.final.log
  mv zt3_1_rnaLog.out data_preproc/mapped_by_star/zt3_1_rna.log
  mv zt3_1_rnaSJ.out.tab data_preproc/mapped_by_star/zt3_1_rna_sj.tsv
  rm -rf zt3_1_rna_STARgenome
  rm zt3_1_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt3_1_rna.sort.bam
fi


echo '### Processing "zt3_2_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt3_2_rna.fastq.gz) \
  --outFileNamePrefix zt3_2_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt3_2_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt3_2_rna.sort.bam
  rm zt3_2_rnaAligned.out.bam
  mv zt3_2_rnaLog.final.out data_preproc/mapped_by_star/zt3_2_rna.final.log
  mv zt3_2_rnaLog.out data_preproc/mapped_by_star/zt3_2_rna.log
  mv zt3_2_rnaSJ.out.tab data_preproc/mapped_by_star/zt3_2_rna_sj.tsv
  rm -rf zt3_2_rna_STARgenome
  rm zt3_2_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt3_2_rna.sort.bam
fi


echo '### Processing "zt6_1_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt6_1_rna.fastq.gz) \
  --outFileNamePrefix zt6_1_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt6_1_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt6_1_rna.sort.bam
  rm zt6_1_rnaAligned.out.bam
  mv zt6_1_rnaLog.final.out data_preproc/mapped_by_star/zt6_1_rna.final.log
  mv zt6_1_rnaLog.out data_preproc/mapped_by_star/zt6_1_rna.log
  mv zt6_1_rnaSJ.out.tab data_preproc/mapped_by_star/zt6_1_rna_sj.tsv
  rm -rf zt6_1_rna_STARgenome
  rm zt6_1_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt6_1_rna.sort.bam
fi


echo '### Processing "zt6_2_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt6_2_rna.fastq.gz) \
  --outFileNamePrefix zt6_2_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt6_2_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt6_2_rna.sort.bam
  rm zt6_2_rnaAligned.out.bam
  mv zt6_2_rnaLog.final.out data_preproc/mapped_by_star/zt6_2_rna.final.log
  mv zt6_2_rnaLog.out data_preproc/mapped_by_star/zt6_2_rna.log
  mv zt6_2_rnaSJ.out.tab data_preproc/mapped_by_star/zt6_2_rna_sj.tsv
  rm -rf zt6_2_rna_STARgenome
  rm zt6_2_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt6_2_rna.sort.bam
fi


echo '### Processing "zt12_1_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt12_1_rna.fastq.gz) \
  --outFileNamePrefix zt12_1_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt12_1_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt12_1_rna.sort.bam
  rm zt12_1_rnaAligned.out.bam
  mv zt12_1_rnaLog.final.out data_preproc/mapped_by_star/zt12_1_rna.final.log
  mv zt12_1_rnaLog.out data_preproc/mapped_by_star/zt12_1_rna.log
  mv zt12_1_rnaSJ.out.tab data_preproc/mapped_by_star/zt12_1_rna_sj.tsv
  rm -rf zt12_1_rna_STARgenome
  rm zt12_1_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt12_1_rna.sort.bam
fi


echo '### Processing "zt12_2_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt12_2_rna.fastq.gz) \
  --outFileNamePrefix zt12_2_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt12_2_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt12_2_rna.sort.bam
  rm zt12_2_rnaAligned.out.bam
  mv zt12_2_rnaLog.final.out data_preproc/mapped_by_star/zt12_2_rna.final.log
  mv zt12_2_rnaLog.out data_preproc/mapped_by_star/zt12_2_rna.log
  mv zt12_2_rnaSJ.out.tab data_preproc/mapped_by_star/zt12_2_rna_sj.tsv
  rm -rf zt12_2_rna_STARgenome
  rm zt12_2_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt12_2_rna.sort.bam
fi


echo '### Processing "zt18_1_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt18_1_rna.fastq.gz) \
  --outFileNamePrefix zt18_1_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt18_1_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt18_1_rna.sort.bam
  rm zt18_1_rnaAligned.out.bam
  mv zt18_1_rnaLog.final.out data_preproc/mapped_by_star/zt18_1_rna.final.log
  mv zt18_1_rnaLog.out data_preproc/mapped_by_star/zt18_1_rna.log
  mv zt18_1_rnaSJ.out.tab data_preproc/mapped_by_star/zt18_1_rna_sj.tsv
  rm -rf zt18_1_rna_STARgenome
  rm zt18_1_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt18_1_rna.sort.bam
fi


echo '### Processing "zt18_2_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt18_2_rna.fastq.gz) \
  --outFileNamePrefix zt18_2_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt18_2_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt18_2_rna.sort.bam
  rm zt18_2_rnaAligned.out.bam
  mv zt18_2_rnaLog.final.out data_preproc/mapped_by_star/zt18_2_rna.final.log
  mv zt18_2_rnaLog.out data_preproc/mapped_by_star/zt18_2_rna.log
  mv zt18_2_rnaSJ.out.tab data_preproc/mapped_by_star/zt18_2_rna_sj.tsv
  rm -rf zt18_2_rna_STARgenome
  rm zt18_2_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt18_2_rna.sort.bam
fi


echo '### Processing "zt21_1_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt21_1_rna.fastq.gz) \
  --outFileNamePrefix zt21_1_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt21_1_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt21_1_rna.sort.bam
  rm zt21_1_rnaAligned.out.bam
  mv zt21_1_rnaLog.final.out data_preproc/mapped_by_star/zt21_1_rna.final.log
  mv zt21_1_rnaLog.out data_preproc/mapped_by_star/zt21_1_rna.log
  mv zt21_1_rnaSJ.out.tab data_preproc/mapped_by_star/zt21_1_rna_sj.tsv
  rm -rf zt21_1_rna_STARgenome
  rm zt21_1_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt21_1_rna.sort.bam
fi


echo '### Processing "zt21_2_rna" ###'
STAR \
  --runThreadN 8 \
  --genomeDir idx/idx_star_rna \
  --alignIntronMin 20 \
  --alignIntronMax 3000 \
  --outFilterScoreMinOverLread 0.1 \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 10 \
  --outSAMtype BAM Unsorted \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --readFilesIn <(gzip -dc data_preproc/fastq_qf/zt21_2_rna.fastq.gz) \
  --outFileNamePrefix zt21_2_rna

if [ $? = 0 ]; then
  samtools view -uh -F 256 zt21_2_rnaAligned.out.bam | \
    samtools sort -@ 8 > data_preproc/mapped_by_star/zt21_2_rna.sort.bam
  rm zt21_2_rnaAligned.out.bam
  mv zt21_2_rnaLog.final.out data_preproc/mapped_by_star/zt21_2_rna.final.log
  mv zt21_2_rnaLog.out data_preproc/mapped_by_star/zt21_2_rna.log
  mv zt21_2_rnaSJ.out.tab data_preproc/mapped_by_star/zt21_2_rna_sj.tsv
  rm -rf zt21_2_rna_STARgenome
  rm zt21_2_rnaLog.progress.out
  samtools index data_preproc/mapped_by_star/zt21_2_rna.sort.bam
fi


