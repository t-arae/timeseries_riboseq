#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/readcount/count_rna_exon

echo '### Processing "zt0_1_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt0_1_rna_gene_counts.txt \
  data_preproc/bam_psite/zt0_1_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt0_1_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt0_2_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt0_2_rna_gene_counts.txt \
  data_preproc/bam_psite/zt0_2_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt0_2_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt3_1_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt3_1_rna_gene_counts.txt \
  data_preproc/bam_psite/zt3_1_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt3_1_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt3_2_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt3_2_rna_gene_counts.txt \
  data_preproc/bam_psite/zt3_2_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt3_2_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt6_1_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt6_1_rna_gene_counts.txt \
  data_preproc/bam_psite/zt6_1_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt6_1_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt6_2_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt6_2_rna_gene_counts.txt \
  data_preproc/bam_psite/zt6_2_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt6_2_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt12_1_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt12_1_rna_gene_counts.txt \
  data_preproc/bam_psite/zt12_1_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt12_1_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt12_2_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt12_2_rna_gene_counts.txt \
  data_preproc/bam_psite/zt12_2_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt12_2_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt18_1_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt18_1_rna_gene_counts.txt \
  data_preproc/bam_psite/zt18_1_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt18_1_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt18_2_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt18_2_rna_gene_counts.txt \
  data_preproc/bam_psite/zt18_2_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt18_2_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt21_1_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt21_1_rna_gene_counts.txt \
  data_preproc/bam_psite/zt21_1_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt21_1_rna_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt21_2_rna" ###'
featureCounts \
  -O \
  -M \
  -t exon \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  -o data_preproc/readcount/count_rna_exon/zt21_2_rna_gene_counts.txt \
  data_preproc/bam_psite/zt21_2_rna.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_rna_exon/zt21_2_rna_gene_counts.log \
  3>&2 2>&1 1>&3



