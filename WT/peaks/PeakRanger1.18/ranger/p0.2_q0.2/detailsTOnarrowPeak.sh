#!/bin/bash

for i in 1 3 8
do
  tail -n +24 WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_details > WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead
  grep -v 'chloroplast\|mitochondria' WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead > WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead_tmp1
  awk 'BEGIN {OFS="\t"}; {$1 = "chr"$1; print}' WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead_tmp1 > WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead_tmp2
  awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $8, $8, $8, $9, $6, $7, $5}' WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead_tmp2 > WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead_tmp3
  awk 'BEGIN {OFS="\t"}; {$4 = "."; $5 = "."; $6 = "."; $10 = $10-$2; print}' WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead_tmp3 > WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2.narrowPeak
  rm WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_noHead*
done
