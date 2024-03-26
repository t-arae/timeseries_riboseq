#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/readcount/count_ribo_central_cds_psite

echo '### Processing "zt0_1_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt0_1_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt0_1_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt0_1_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt0_2_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt0_2_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt0_2_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt0_2_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt3_1_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt3_1_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt3_1_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt3_1_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt3_2_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt3_2_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt3_2_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt3_2_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt6_1_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt6_1_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt6_1_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt6_1_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt6_2_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt6_2_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt6_2_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt6_2_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt12_1_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt12_1_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt12_1_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt12_1_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt12_2_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt12_2_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt12_2_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt12_2_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt18_1_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt18_1_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt18_1_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt18_1_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt18_2_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt18_2_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt18_2_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt18_2_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt21_1_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt21_1_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt21_1_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt21_1_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



echo '### Processing "zt21_2_ribo" ###'
featureCounts \
  -O \
  -M \
  -t CDS \
  -g gene_id \
  -T 4 \
  -a misc/gff_gtf/Araport11_central_cds.gtf \
  -o data_preproc/readcount/count_ribo_central_cds_psite/zt21_2_ribo_gene_counts.txt \
  data_preproc/bam_psite/zt21_2_ribo.sort.bam \
  3>&2 2>&1 1>&3 | tee data_preproc/readcount/count_ribo_central_cds_psite/zt21_2_ribo_gene_counts.log \
  3>&2 2>&1 1>&3



