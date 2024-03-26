#!/bin/bash
set -eu
set -o pipefail

mkdir -p analysis/phase_analysis_ejtk/ribo_uorf_psite

eJTK-CalcP.py \
  --filename analysis/phase_analysis_ejtk/ribo_uorf_psite/data.txt \
  --waveform analysis/phase_analysis_ejtk/waveform_cosine.txt \
  --period analysis/phase_analysis_ejtk/period24.txt \
  --phase analysis/phase_analysis_ejtk/phases_00-21_by3.txt \
  --asymmetry analysis/phase_analysis_ejtk/asymmetries_03-21_by3.txt \
  --prefix ribo_uorf_psite

