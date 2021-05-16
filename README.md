# Subassembly_PB_IndelCorrection

Correcting indels to PacBio reads for subassembly

#### Running instructions on GS Cluster
qsub ./driver_SUL1_test.sh

#### Running locally

Scripts are now edited to run locally. To run on the cluster, download an older version (before 10/2/2020).  
To set up your environment, install [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) and ensure it is working on your computer. Load the environment `conda env create --file msaccs_env.yml`  


(Updated 5/15/2021)
The msa_pacbio.py script utilizes the multiple sequence aligner, [MUSCLE](https://www.drive5.com/muscle/downloads.htm) and Needle through 
[EMBOSS](http://emboss.sourceforge.net/docs/adminguide/node7.html).  
To install MUSCLE, extract the zipped MUSCLE file `tar -zxvf muscle_filename_here.tar.gz`. The MUSCLE software is ready to run as soon as it is unzipped. For the EMBOSS file, unzip it `tar -zxvf emboss_filename_here.tar.gz`. Through terminal, go to the unzipped EMBOSS directory and type `./configure`. When that is finished, type `make` (this may take 5-10 minutes). The software needed through EMBOSS is called needle, which is located in `EMBOSS-versionX/emboss/needle`.  


In your driver script, be sure to specify the location of each software.
  
| Option | Description |
|--------|-------------|
| **-d**,**--directory**	| Specify working directory where intermediate and output files will be located |
| **-o**,**--out** |	Specify the name of the final output file  (default = Seq_barcodes_aligned.txt) |
| **--highQual** | File of barcode-variant association, where the variant is the highest quality read  |
| **--inputSeqs** | Raw barcode, variant, and quality of sequencies |
| **-c**,**--cutoff** |	Minimum number of CCS reads needed in order to retain reads associated with specific barcodes (default = 2) |
| **-t**,**--threshold** |	Minimum frequency threshold for calling consensus reads (default = 0.7) |
| **-s**,**--stats**  | Option to generate alignment stats. Currently outputs below_threshold_Ncount.txt file which returns the barcodes with ambiguous sites and how many ambiguous sites were present total. |
| **-v**,**--verbose** |	Debug output |
| **-m**,**--muscle** | Location of compiled/extracted MUSCLE program |
| **-n**,**--needle** | Location of compiled Needle program |
| **--cont** | If program is disrupted or aborts, this feature will allow user to continue with unprocessed reads. Previously processed reads will not be reprocessed |

**Input files:**

This script only requires two input files. The input files for both `--highQual` and `--inputSeqs` file should be tab-delimited file, where the first column is the barcode and the second column is the associated read.

**Output files:**
  * The output file, dictated by the file name you specify for `--out`, is a tab-delimited file, where the first column is the barcode and the second column is the aligned associated read. 
  * The script will also produce a file called `progress_file.txt`. This file is generated when a barcode is finished being processed. If your script is killed in the middle of a run, you can include the `--cont` option the next time you run the driver script, and it will not reprocess barcodes that were already processed in the last run.
  * The `barcodes_below_cutoff.txt` file returns barcodes that were not processed because it did not meet the read cutoff specified by `-c`. You can double check that the script is working properly by adding together the number of barcodes in your output file (from `--out`) and the number of barcodes in the `barcodes_below_cutoff.txt` file.
  * The `below_threshold_Ncount.txt` file returns the number of barcodes that contained ambiguous sites if the `--stats` option is included. Sites are determined to be ambiguous if nucleotides in the same position do not pass the majority threshold specified by the `--threshold` parameter. If the `--stats` option is not included, the `below_threshold_Ncount.txt` file will be deleted.
  * `run_times.tsv` returns the start time, end time, and number of cores for a particular job. This is appended to the file.
 

