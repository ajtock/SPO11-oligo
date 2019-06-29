#!/bin/bash

csmit -m 30G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_MYC_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1" & sleep 20;
csmit -m 30G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep2" & sleep 20;
csmit -m 30G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep1" & sleep 20;
csmit -m 30G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab.bed MNase" & sleep 20;
csmit -m 30G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1" & sleep 20;
csmit -m 30G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI3" & sleep 20;
csmit -m 30G -c 1 "./peak_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI8" & sleep 20;

wait

