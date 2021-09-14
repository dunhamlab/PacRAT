#!/bin/bash
#$ -S /bin/bash
#$ -l mfree=4G -pe serial 16 -l h_rt=0:48:0:0
#$ -cwd
#$ -N PacRAT

## There are two sections in this script; please be sure to comment/uncomment as appropriate for your situation






## ******* Section 1: For CentOS7 on the Genome Sciences cluster ******* #
## Be sure the comment this section out if you are running locally or on a different cluster environment
module load python/3.7.7
module load numpy/1.19.2
module load biopython/1.77
module load joblib/0.15.1

python msa_pacbio.py -d ./output -o H2B_barcode_variant_map_msa.txt \
	--highQual ../input/H2B_highQual_seqs_1000_barcodes.tsv \
	--inputSeqs ../input/H2B_reads_1000_barcodes.txt \
	-c 1 -t 0.6 -s \
	-m /net/gs/vol3/software/modules-sw/muscle/3.8.31/Linux/CentOS7/x86_64/bin/muscle \
	-n /net/gs/vol3/software/modules-sw/EMBOSS/6.6.0/Linux/CentOS7/x86_64/bin/needle






## ******* Section 2: Running locally with conda environment ******* #
## This is a sample for running PacRAT locally. We have here a skeleton of necessary steps you need to take in order to run PacRAT. Customize as needed.
## Please look at README file for directions on installing MUSCLE and NEEDLE. You must specify the location of the software.
## Uncomment this section if you're choosing to run locally

## Validate that Anaconda/Miniconda are installed
# SOURCECONDA=$(conda info --base)
# if [ -z "$SOURCECONDA" ]
# then
#	echo "Conda not available. Please install Anaconda or Miniconda https://docs.conda.io/en/latest/miniconda.html"
#	exit 1
# fi

## Install and activate conda environment
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

## example for running locally
## change directory location for the options -m and -n
## files specified in --highQual and --inputSeqs are for the test case; be sure to change it for your specific library.
## unzip H2B_seq_barcodes.txt.gz
# python msa_pacbio.py -d ./output -o H2B_barcode_variant_map_msa.txt \
# 	--highQual ../input/H2B_combined_minQ0_assignment.tsv \
# 	--inputSeqs ../input/H2B_seq_barcodes.txt \
#	-c 1 -t 0.6 -s \
# 	-m ../../muscle/muscle -n ../../emboss/EMBOSS-6.6.0/emboss/needle
