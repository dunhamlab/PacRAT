#!/bin/bash
#$ -S /bin/bash
#$ -l mfree=4G -pe serial 16 -l h_rt=0:24:0:0
#$ -cwd
#$ -N msa_pacbio

# example error/output file locations: for Clara! CHANGEEEE
# -o /net/dunham/vol2/Clara/projects/Subassembly_PB_IndelCorrection/test
# -e /net/dunham/vol2/Clara/projects/Subassembly_PB_IndelCorrection/test


SOURCECONDA=$(conda info --base)
if [ -z "$SOURCECONDA" ]
then
	echo "Conda not available. Please install Anaconda or Miniconda https://docs.conda.io/en/latest/miniconda.html"
	exit 1
fi

ENVCHECK=$(conda env list | grep "msa_ccs")
source $SOURCECONDA/etc/profile.d/conda.sh
if [ -z "$ENVCHECK" ]
then
	echo "msa_ccs environment not installed! installing now..."
	conda env create --file msaccs_env.yml
	echo "msa_ccs environment installed"
else
	echo "msa_ccs environment is installed!"
fi

conda activate msa_ccs
echo "Environment activated"

# for running in cluster
python msa_pacbio.py -d ./test \
	-o SUL1_test100_subassembly.txt \
	--highQual SUL1_test100_combined_minQ0_assignment.tsv \
	--inputSeqs SUL1_test100_seq_barcodes_filtered.txt -c 2 -t 0.6 \
	-m /net/gs/vol3/software/modules-sw/muscle/3.8.31/Linux/CentOS7/x86_64/bin/muscle \
	-n /net/gs/vol3/software/modules-sw/EMBOSS/6.6.0/Linux/CentOS7/x86_64/bin/needle


# example for running on desktop
# python msa_pacbio.py -d ./test \
# 	-o SUL1_test100_seq_barcodes_aligned.txt \
# 	--highQual SUL1_test100_combined_minQ0_assignment.tsv \
# 	--inputSeqs SUL1_test100_seq_barcodes_filtered.txt -c 2 -t 0.6 -v \
# 	-m ../../muscle/muscle -n ../../emboss/EMBOSS-6.6.0/emboss/needle
