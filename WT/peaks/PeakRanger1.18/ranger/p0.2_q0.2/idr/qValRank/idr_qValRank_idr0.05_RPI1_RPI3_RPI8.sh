# Perform IDR analysis to select a reproducible set of peaks across biological replicates

peakDir=/projects/ajt200/BAM_masters/SPO11-oligo/WT/SPO11_peaks/PeakRanger1.18/ranger/p0.2_q0.2

for i in 3 8
do
  idr --samples $peakDir/WT_SPO11oligos_RPI1_peaks_peakranger_ranger_p0.2_q0.2_TreadsNormCreads.narrowPeak \
                $peakDir/WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2_TreadsNormCreads.narrowPeak \
      --input-file-type narrowPeak --output-file-type narrowPeak \
      --output-file WT_SPO11oligos_RPI1_RPI${i}_idrValues_qValRank_idr0.05 \
      --rank q.value --idr-threshold 0.05 --plot \
      --log-output-file WT_SPO11oligos_RPI1_RPI${i}_idrValues_qValRank_idr0.05.log
  awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' WT_SPO11oligos_RPI1_RPI${i}_idrValues_qValRank_idr0.05 > WT_SPO11oligos_RPI1_RPI${i}_idrValues_qValRank_idr0.05.narrowPeak
done

idr --samples WT_SPO11oligos_RPI1_RPI3_idrValues_qValRank_idr0.05.narrowPeak \
              WT_SPO11oligos_RPI1_RPI8_idrValues_qValRank_idr0.05.narrowPeak \
    --input-file-type narrowPeak --output-file-type narrowPeak \
    --output-file WT_SPO11oligos_RPI1_RPI3_RPI8_idrValues_qValRank_idr0.05 \
    --rank q.value --idr-threshold 0.05 --plot \
    --log-output-file WT_SPO11oligos_RPI1_RPI3_RPI8_idrValues_qValRank_idr0.05.log
awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' WT_SPO11oligos_RPI1_RPI3_RPI8_idrValues_qValRank_idr0.05 > WT_SPO11oligos_RPI1_RPI3_RPI8_idrValues_qValRank_idr0.05.narrowPeak

