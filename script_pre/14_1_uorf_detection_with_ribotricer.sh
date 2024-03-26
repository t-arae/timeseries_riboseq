#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_preproc/ribotricer_out

ribotricer prepare-orfs \
  --gtf misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf \
  --fasta misc/fasta/BSgenome.Athaliana.TAIR.TAIR9.fasta \
  --prefix data_preproc/ribotricer_out/araport11 \
  --min_orf_length 3 \
  --start_codons ATG \
  --stop_codons TAG,TAA,TGA \
  --longest

