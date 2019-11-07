#!/bin/bash
#$ -S /bin/bash
#$ -l mfree=16G -l h_rt=0:24:0:0
#$ -cwd
#$ -N dbl_msa_pacbio
#$ -o /net/dunham/vol2/Cindy/pacbio_git/Subassembly_PB_IndelCorrection/test
#$ -e /net/dunham/vol2/Cindy/pacbio_git/Subassembly_PB_IndelCorrection/test

# for Clara! CHANGEEEE
# -o /net/dunham/vol2/Clara/projects/Subassembly_PB_IndelCorrection/test
# -e /net/dunham/vol2/Clara/projects/Subassembly_PB_IndelCorrection/test

# module load python/2.7.3 #3.5.2
# module load numpy/1.7.0
# module load biopython/1.63
# module load muscle/3.8.31

module load python/3.5.2
module load muscle/3.8.31

# for Clara :) 
#python msa_pacbio.py -d /net/dunham/vol2/Clara/projects/Subassembly_PB_IndelCorrection/test -o SUL1_test100_seq_barcodes_aligned.txt --highQual SUL1_test100_combined_minQ0_assignment.tsv --inputSeqs SUL1_test100_seq_barcodes_filtered.txt

# for Cindy :)
python msa_pacbio.py -d /net/dunham/vol2/Cindy/pacbio_git/Subassembly_PB_IndelCorrection/test -o SUL1_test100_seq_barcodes_aligned.txt --highQual SUL1_test100_combined_minQ0_assignment.tsv --inputSeqs SUL1_test100_seq_barcodes_filtered.txt -v
