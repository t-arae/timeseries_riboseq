#!/bin/bash
set -eu
set -o pipefail

mkdir -p data_modified/gff_gtf

gtfToGenePred misc/gff_gtf/Araport11_GFF3_genes_transposons.201606_mod.gtf data_modified/gff_gtf/transcript.GenePred
genePredToBed -tab data_modified/gff_gtf/transcript.GenePred temp.txt
tr -d ';' < temp.txt > data_modified/gff_gtf/transcript.bed12

