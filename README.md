# PacBio Read Alignment Tool (PacRAT)

![pacRAT_final_small](https://user-images.githubusercontent.com/50029049/129951417-ef55a258-e6ba-40c8-8e97-eb1aca1e8d01.png)

Improving PacBio barcode-variant mapping (subassembly) using multiple sequence alignment.

### Implementation

To use, git clone or download and unzip the PacRAT code. 

`git clone https://github.com/dunhamlab/PacRAT.git`

#### Running on UW GS Cluster

From the PacRAT directory, run: 

`qsub ./driver_pacrat.sh`

to run PacRAT with the example data provided in this repository. To run with different data, update the `driver_pacrat.sh` script with the location of the input and output files, described below. 

#### Running locally

Scripts can run locally, although we recommend using a cluster/job submission system to optimize memory usage. 
To set up your environment, install [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) and ensure it is working on your computer. Once Anaconda is installed, you should be able to run the `driver_pacrat.sh` script. If you run the Python script on its own (without using the `driver_pacrat.sh` script, be sure to load the environment `conda env create --file pacrat_env.yml` and `conda activate pacrat_env` beforehand to ensure all packages are installed, compatible, and activated.

In your driver script (you can use `driver_pacrat.sh` to run your program; be sure to comment out the appropriate sections), specify the location of the input files (`--highQual` and `--inputSeqs`), as well as your working directory (`-d`) and output file (`-o`).

**Parameter Descriptions**
  
| Option | Description |
|--------|-------------|
| **-d**,**--directory**	| Specify working directory where intermediate and output files will be located |
| **-o**,**--out** |	Specify the name of the final output file  (default = H2B_barcode_variant_map_msa.txt) |
| **--highQual** | File of barcode-variant association, where the variant is the highest quality read  |
| **--inputSeqs** | Raw barcode, variant, and quality of sequences |
| **-c**,**--cutoff** |	Minimum number of CCS reads needed in order to retain reads associated with specific barcodes (default = 1) |
| **-t**,**--threshold** |	Minimum frequency threshold for calling consensus reads (default = 0.6) |
| **-s**,**--stats**  | Option to generate alignment stats. Currently outputs below_threshold_Ncount.txt, barcodes_below_cutoff.txt files, and ccs_count_per_barcode.txt; see below for details |
| **-v**,**--verbose** |	Print verbose debug output |
| **--cont** | If program is disrupted or aborts, enabling this feature will allow user to continue with unprocessed reads. Previously processed reads will not be reprocessed |
| **-r**,**--rmint** | Removes intermediate alignment files |

**Input files:**

This script requires two input files. The input files for both `--highQual` and `--inputSeqs` file should be tab-delimited file, where the first column is the barcode and the second column is the associated read. You can generate both these files following the pipeline described in the [AssemblyByPacBio repository](https://github.com/shendurelab/AssemblyByPacBio). The `--inputSeqs` file will be generated from the `extractBarcodeInsertPairs_moreQC.py` script, and the `--highQual` file will be generated from the `unifyAssignment.py` script.

**Output files:**
  * The output file, dictated by the file name you specify for `--out`, is a tab-delimited file, where the first column is the barcode and the second column is the aligned associated read. 
  * The script will also produce a file called `progress_file.txt`. This file is generated when a barcode is finished being processed. If your script is killed in the middle of a run, you can include the `--cont` option the next time you run the driver script, and it will not reprocess barcodes that were already processed in the last run.
  * The `barcodes_below_cutoff.txt` file lists barcodes that were not processed because they did not meet the read cutoff specified by `--cutoff`, if the `--stats` option is included. You can double check that the script is working properly by adding together the number of barcodes in your output file (from `--out`) and the number of barcodes in the `barcodes_below_cutoff.txt` file. If the `--stats` option is not included, the `barcodes_below_cutoff.txt` file will not be generated.
  * The `below_threshold_Ncount.txt` file returns the number of barcodes that contained ambiguous sites if the `--stats` option is included. Sites are determined to be ambiguous if nucleotides in the same position do not pass the majority threshold specified by the `--threshold` parameter. If the `--stats` option is not included, the `below_threshold_Ncount.txt` file will not be generated.
  * The `ccs_count_per_barcode.txt` file returns the number of CCS reads associated with each barcode. File is created if the `--stats` option is included. If the `--stats` option is not included, the `ccs_count_per_barcode.txt` file will not be generated.
 
 **Example input/output files:**
 
 The example input files provided are located in the `input` folder and hold the `--highQual` and `--inputSeqs` files from a subset (1000 barcodes) of the simulated H2B library. These example files run in ~ 1 min on the GS cluster. 
