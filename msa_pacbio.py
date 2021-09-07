import os
import os.path
import numpy as np
from Bio import SeqIO
import glob
import sys
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from optparse import OptionParser
from datetime import datetime
from joblib import Parallel, delayed
import multiprocessing
import getpass

# Print timing info
startTime = str(datetime.now())
print("Starting time: "+startTime)

# Option Parser: required arguments are two input files (raw barcode-sequence file, and high quality sequence file)
parser = OptionParser()
parser.add_option("-d","--directory", dest="workdir", help="Working directory",default=os.getcwd(),type="string")
parser.add_option("-o","--out", dest="out", help="Output file",default="Seq_barcodes_aligned.txt",type="string")
parser.add_option("--highQual", dest="highQualFile", help="File of barcode-seq association, seq from highest quality read",type="string")
parser.add_option("--inputSeqs", dest="inputSeqsFile", help="Raw barcode, sequence, quality input sequences",type="string")
# Additional options
parser.add_option("-c","--cutoff", dest="cutoff", help="Minimum number of ccs reads for analysis",default=1,type="int")
parser.add_option("-t","--threshold", dest="thresh", help="Minimum threshold to determine consensus sequence",default=0.6,type="float")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_option("-m","--muscle", dest="muscle", help="Compiled MUSCLE program",default="./muscle",type="string")
parser.add_option("-n","--needle", dest="needle", help="Compiled NEEDLE program",default="./needle",type="string")
parser.add_option("--cont", dest="cont", help="Continue working after disrupted run",default=False,action="store_true")
parser.add_option("-s", "--stats", dest="stats", help="Get stats for barcodes that need realignment",default=False,action="store_true")
parser.add_option("-r","--rmint", dest="rm_intermediates",help="Removes intermediate files when finished",default=False,action="store_true")

(options, args) = parser.parse_args()
muscle_exe = options.muscle
needle_exe = options.needle

os.chdir(options.workdir) # change working directory to output folder

# **************** Read input files ******************************************************* #
print("Program started, reading input files")
# read original assignments into dictionary
hq_dict = {}
assignments = open(options.highQualFile, "r") # min_Q0_assignment.tsv
for line in assignments:
	line = line.upper()
	paired_bcread = line.strip().split()
	hq_dict[paired_bcread[0]] = paired_bcread[1]
assignments.close()
totalBarcodes = len(hq_dict.keys())
if options.verbose: print(str(totalBarcodes) + " barcodes found in hq file")

# Read all ccs reads into dict: BC: [seq1, seq2...]. Seq quality pairs stored as tuples
read_dict = {}
reads = open(options.inputSeqsFile, "r") # seq_barcodes.txt
for line in reads: 
	line = line.upper()
	paired_bcread = line.strip().split()
	if paired_bcread[0] in read_dict:
		read_dict[paired_bcread[0]].append(paired_bcread[1])
	else:
		read_dict[paired_bcread[0]] = [paired_bcread[1]]
reads.close()
totalBarcodes2 = len(read_dict.keys())
if options.verbose: print(str(totalBarcodes2) + " barcodes found in other file") 
# ***************************************************************************************** #

# create intermediates directories in output folder
os.system("mkdir -p intermediates & mkdir -p intermediates/fasta & mkdir -p intermediates/alignments & mkdir -p intermediates/fasta_2 & mkdir -p intermediates/realignments") 

# ************* if needed to pick up where script left off (option --cont) **************** #
# deletes progress file if --cont is not included; removes barcodes from hq_dict if --cont is included
progress_file_name = "progress_file.txt" # create progress file that writes barcode each time loop_bcs function is run
if options.cont:
	outputfile = open(options.out,"a")
	if os.path.exists(progress_file_name):
		prog_file = open(progress_file_name,"r+")
		remove_keys = [] # list of keys to remove from dict because it's already been processed
		for line in prog_file:
			line = line.strip()
			if line != "":
				remove_keys.append(line)
		# remove barcodes in hq_dict if it's already been processed through loop_bcs function
		for bc in remove_keys:
			if bc in hq_dict:
				hq_dict.pop(bc, None)
		print("Continue option enabled, " + str(len(remove_keys)) + " barcodes already processed, continuing.")
	else: # If no progress file exists, print error and exit
		sys.exit("Continue option specified but no progress file found, exiting")

else: # if --continue is not included, create/rewrite the progress file from the last run.
	outputfile = open(options.out, "w")
	prog_file = open(progress_file_name, "w")
	
# ***************************************************************************************** #
	
# ************** Stats option: open/append relevant files (option -s, --stats) ************ #
cutoff_filename = "barcodes_below_cutoff.txt" # barcodes that miss the cutoff
threshold_filename = "below_threshold_Ncount.txt" # barcodes that have an ambigious site in sequence
if options.stats:
	if options.verbose: print("Creating stats files")
	if options.cont:
		cutoff_file = open(cutoff_filename, "a")
		threshold_file = open(threshold_filename, "a")
	else:
		cutoff_file = open(cutoff_filename, "w")
		threshold_file = open(threshold_filename, "w")
# ***************************************************************************************** #

# **************** Main alignment/Consensus function ************************************** #
# Align all seuqences with the same barcode. If consensus found, use that sequence. 
# If not, realign to high quality sequence
def loop_bcs(key):
	bc_entry = read_dict[key] #list of sequences
	#create fasta file for each barcode: 
	fasta_file_name = os.path.join("intermediates/fasta/" + key + ".fasta") 
	fasta_file = open(fasta_file_name, "w")
	i = 0       
	for item in bc_entry: #add each read for a particular barcode in fasta format
		fasta_file.write(">" + key + "_" + str(i) + "\n")
		fasta_file.write(item+"\n")
		i += 1
	fasta_file.close()
	if options.verbose: print("Made fasta file for BC: " + str(key))

	# only align if there are at least CUTOFF ccs reads
	if len(bc_entry) >= options.cutoff:
		if len(bc_entry) == 1: #special case, don't need to align here
			outputfile.write(key+"\t"+bc_entry[0]+"\n")
			
		else: 
			#align files together - first alignment
			aln_file_name = "intermediates/alignments/" + key + ".aln"
			#muscle system call here, write to output file
			muscle_cline = MuscleCommandline(muscle_exe, input=fasta_file_name, out=aln_file_name)
			stdout, stderr = muscle_cline(fasta_file_name)
			if options.verbose: print("Passed cutoff, made first alignment " + str(key))

			#get consensus: 
			consensus = ""
			if os.path.exists(aln_file_name): #probably should have an else here that throws an error, but it would be better to check that the muscle stderr is empty
				alignment = AlignIO.read(aln_file_name, 'fasta')
				summary_align = AlignInfo.SummaryInfo(alignment)
				consensus = summary_align.gap_consensus(threshold=options.thresh,  ambiguous='N')
				consensus = str(consensus)
				consensus = consensus.replace("-","") 
			else: sys.exit("Could not find alignment file for BC: ", str(key))

			#if N's: realign (pairwise aligner w/in python) to highest qual, and find consensus from that
			if 'N' in consensus:
				# List barcodes that have ambigious sites, if --stats is included
				if options.stats:
					Ncount = consensus.count("N")
					threshold_file.write(key + "\t" + str(Ncount) + "\n")
					threshold_file.flush()
					if options.verbose: print(key + " barcode has " + str(Ncount) + " ambiguous sites.")
				
				#write 1st consensus and HQ read to new file
				int_file_name_2 = os.path.join("intermediates/fasta_2/" + key + ".fasta")
				fasta_2 = open(int_file_name_2,"w")
				fasta_2.write(">"+key+"\n"+consensus+"\n")
				fasta_2.close()
				fasta_hq_name = os.path.join("intermediates/fasta_2/" + key + "_hq.fasta")
				fasta_hq = open(fasta_hq_name,"w")
				fasta_hq.write(">"+key+"_hq\n"+hq_dict[key]+"\n")
				fasta_hq.close()
				aln_file_name_2 = "intermediates/realignments/"+key+".aln"
				cmd = needle_exe + " " + int_file_name_2 + " " + fasta_hq_name + " -auto -gapopen 10 -gapextend 0.5 -outfile " + aln_file_name_2 + " -aformat fasta"
				os.system(cmd)

				#consensus of new alignment file (mostly same from previous script)
				alignment_2 = list(SeqIO.parse(aln_file_name_2,"fasta"))
				consensus_seq = str(alignment_2[0].seq)
				hq_seq = str(alignment_2[1].seq)
				finalSeq = ""
				for i in range(len(consensus_seq)):
					if consensus_seq[i] == "N":
						finalSeq = finalSeq+hq_seq[i]
					else:
						finalSeq = finalSeq + consensus_seq[i]
				consensus = finalSeq
				consensus = consensus.replace("-","")		
				outputfile.write(key+"\t"+consensus+"\n")
				if options.verbose: print("Realigned and got new consensus for BC: " + str(key))
	
			#if no Ns: write consensus to output file
			else:
				outputfile.write(key+"\t"+consensus+"\n")
		
	else: # generates a file of barcodes that did not meet the minimum number of reads (under option -c)
		if options.stats:
			cutoff_file.write(key+"\n")
			cutoff_file.flush()

	outputfile.flush()
	# delete fasta files
	if os.path.exists("intermediates/fasta/" + key + ".fasta"):
		os.remove("intermediates/fasta/" + key + ".fasta")
	if os.path.exists("intermediates/fasta_2/" + key + ".fasta"):
		os.remove("intermediates/fasta_2/" + key + ".fasta")
	if os.path.exists("intermediates/fasta_2/" + key + "_hq.fasta"):
		os.remove("intermediates/fasta_2/" + key + "_hq.fasta")
	# if --rmint, remove intermediate alignment files
	if options.rm_intermediates: 
		if os.path.exists("intermediates/alignments/" + key + ".aln"):
			os.remove("intermediates/alignments/" + key + ".aln")
		if os.path.exists("intermediates/realignments/"+key+".aln"):
			os.remove("intermediates/realignments/"+key+".aln")
	prog_file.write(key+"\n")
	if options.verbose: print("Wrote " + key + " to progress file")
	prog_file.flush()
# ***************************************************************************************** #

# **************** Parallelization and print ********************************************** #
num_cores = multiprocessing.cpu_count()
print("Starting alignments with",str(num_cores), "cores")
results = Parallel(n_jobs=(num_cores),prefer="threads")(delayed(loop_bcs)(key) for key in hq_dict)
print("All alignments completed, generating output files")

# removes intermediates folder if --rmint is flagged
if options.rm_intermediates:
	os.system("rm -r intermediates")
else:
	os.system("rm -r intermediates/fasta*")
	
# close output files
outputfile.close()
prog_file.close()
if options.stats:
	cutoff_file.close()
	threshold_file.close()

endTime = str(datetime.now())
print("Ending time: "+endTime)
# ***************************************************************************************** #
