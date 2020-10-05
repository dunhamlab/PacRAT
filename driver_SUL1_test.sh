#!/bin/bash
#$ -S /bin/bash
#$ -l mfree=16G -pe serial 4 -l h_rt=0:24:0:0
#$ -cwd
#$ -N msa_pacbio

# for Clara! CHANGEEEE
# -o /net/dunham/vol2/Clara/projects/Subassembly_PB_IndelCorrection/test
# -e /net/dunham/vol2/Clara/projects/Subassembly_PB_IndelCorrection/test

# module load python/3.5.2
# module load muscle/3.8.31
# module load EMBOSS/6.4.0

# for Clara :) 
# python msa_pacbio.py -d /net/dunham/vol2/Clara/projects/Subassembly_PB_IndelCorrection/test \
#     -o SUL1_test100_seq_barcodes_aligned.txt \
#     --highQual SUL1_test100_combined_minQ0_assignment.tsv \
#     --inputSeqs SUL1_test100_seq_barcodes_filtered.txt \
#     -c 1 -t 0.7

# for Cindy :)
#python msa_pacbio.py -d /net/dunham/vol2/Cindy/pacbio_git/200518_subassembly_pairwise/Subassembly_PB_IndelCorrection/test -o SUL1_test100_seq_barcodes_aligned.txt --highQual SUL1_test100_combined_minQ0_assignment.tsv --inputSeqs SUL1_test100_seq_barcodes_filtered.txt -c 2 -t 0.6
#/net/gs/vol3/software/modules-sw/EMBOSS/6.6.0/Linux/CentOS7/x86_64/EMBOSS-EMBOSS-6.6.0/emboss

python msa_pacbio.py -d ./test \
	-o SUL1_test100_seq_barcodes_aligned.txt \
	--highQual SUL1_test100_combined_minQ0_assignment.tsv \
	--inputSeqs SUL1_test100_seq_barcodes_filtered.txt -c 2 -t 0.6 \
	-m /net/gs/vol3/software/modules-sw/muscle/3.8.31/Linux/RHEL6/x86_64/bin/muscle \
	-n /net/gs/vol3/software/modules-sw/EMBOSS/6.6.0/Linux/CentOS7/x86_64/EMBOSS-EMBOSS-6.6.0/emboss/needle

# Cindy Desktop
# python msa_pacbio.py -d ./test \
# 	-o SUL1_test100_seq_barcodes_aligned.txt \
# 	--highQual SUL1_test100_combined_minQ0_assignment.tsv \
# 	--inputSeqs SUL1_test100_seq_barcodes_filtered.txt -c 2 -t 0.6 -v \
# 	-m ../../muscle/muscle -n ../../emboss/EMBOSS-6.6.0/emboss/needle
