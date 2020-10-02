# Subassembly_PB_IndelCorrection

Correcting indels to PacBio reads for subassembly

#### Running instructions
qsub ./driver_SUL1_test.sh

#### Running locally

Scripts are now edited to run locally. To run on the cluster, download a previous version. 
The msa_pacbio.py script utilizes the multiple sequence aligner, MUSCLE (https://www.drive5.com/muscle/downloads.htm) and Needle through EMBOSS (ftp://emboss.open-bio.org/pub/EMBOSS/). 
To install these, extract the zipped MUSCLE file `tar -zxvf muscle_filename_here.tar.gz`. The MUSCLE software is ready to run as soon as it is unzipped.
For the EMBOSS file, unzip it `tar -zxvf emboss_filename_here.tar.gz`. Through terminal, go to the unzipped EMBOSS directory and type `./configure`. When that is finished, type `make` (this may take 5-10 minutes). The software needed through EMBOSS is called needle, which is located in `EMBOSS-versionX/emboss/needle`.

In your driver script, be sure to specify the location of each software.

