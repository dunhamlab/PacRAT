# Subassembly_PB_IndelCorrection

Correcting indels to PacBio reads for subassembly

#### Running instructions
qsub ./driver_SUL1_test.sh

#### Running locally

Scripts are now edited to run locally. To run on the cluster, download an older version (before 10/2/2020).  
To set up your environment, install Anaconda and ensure it is working on your computer. Load the environment `conda env create --file msaccs_env.yml`  
The msa_pacbio.py script utilizes the multiple sequence aligner, MUSCLE (https://www.drive5.com/muscle/downloads.htm) and Needle through EMBOSS (ftp://emboss.open-bio.org/pub/EMBOSS/).  
To install MUSCLE, extract the zipped MUSCLE file `tar -zxvf muscle_filename_here.tar.gz`. The MUSCLE software is ready to run as soon as it is unzipped.  
For the EMBOSS file, unzip it `tar -zxvf emboss_filename_here.tar.gz`. Through terminal, go to the unzipped EMBOSS directory and type `./configure`. When that is finished, type `make` (this may take 5-10 minutes). The software needed through EMBOSS is called needle, which is located in `EMBOSS-versionX/emboss/needle`.  
  
In your driver script, be sure to specify the location of each software.
  
| Option | Description |
|--------|-------------|
| **-d**,**--directory**	| Specify working directory where intermediate and output files will be located |
| **-o**,**--out** |	Specify the name of the final output file  (default = Seq_barcodes_aligned.txt) |
| **--highQual** | File of barcode-variant association, where the variant is the highest quality read  |
| **--inputSeqs** | Raw barcode, variant, and quality of sequencies |
| **-c**,**--cutoff** |	Minimum number of CCS reads needed in order to retain reads associated with specific barcodes (default = 2) |
| **-t**,**--threshold** |	Minimum frequency threshold for calling consensus reads (default = 0.7) |
| **-v**,**--verbose** |	Debug output |
| **-m**,**--muscle** | Location of compiled/extracted MUSCLE program |
| **-n**,**--needle** | Location of compiled Needle program |
