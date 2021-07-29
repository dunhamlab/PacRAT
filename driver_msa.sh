#!/bin/bash
#$ -S /bin/bash
#$ -l mfree=4G -pe serial 16 -l h_rt=0:48:0:0
#$ -cwd
#$ -N PacRAT

# For CentOS7 on the cluster
module load python/3.7.7
module load numpy/1.19.2
module load biopython/1.77
module load joblib/0.15.1

# for running in GS cluster
python msa_pacbio.py -d ./output -o SUL1_test100_barcode_variant_map.txt \
	--highQual input/SUL1_test100_combined_minQ0_assignment.tsv \
	--inputSeqs input/SUL1_test100_seq_barcodes_filtered.txt \
	-c 2 -t 0.6 -v \
	-m /net/gs/vol3/software/modules-sw/muscle/3.8.31/Linux/CentOS7/x86_64/bin/muscle \
	-n /net/gs/vol3/software/modules-sw/EMBOSS/6.6.0/Linux/CentOS7/x86_64/bin/needle

# Running locally with conda environment
# SOURCECONDA=$(conda info --base)
# if [ -z "$SOURCECONDA" ]
# then
#	echo "Conda not available. Please install Anaconda or Miniconda https://docs.conda.io/en/latest/miniconda.html"
#	exit 1
# fi

# ENVCHECK=$(conda env list | grep "msa_ccs")
# source $SOURCECONDA/etc/profile.d/conda.sh
# if [ -z "$ENVCHECK" ]
# then
# 	echo "msa_ccs environment not installed! installing now..."
#	conda env create --file msaccs_env.yml
#	echo "msa_ccs environment installed"
# else
#	echo "msa_ccs environment is installed!"
# fi

# conda activate msa_ccs
# echo "Environment activated"

# example for running locally
# python msa_pacbio.py -d ./output -o SUL1_test100_barcode_variant_map.txt \
# 	--highQual input/SUL1_test100_combined_minQ0_assignment.tsv \
# 	--inputSeqs input/SUL1_test100_seq_barcodes_filtered.txt \
#	-c 2 -t 0.6 -v \
# 	-m ../../muscle/muscle -n ../../emboss/EMBOSS-6.6.0/emboss/needle
