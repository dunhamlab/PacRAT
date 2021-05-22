import os
import os.path
import numpy as np
from Bio import SeqIO
import glob
import sys
from Bio.Align.Applications import MuscleCommandline
#from StringIO import StringIO # for Python 2
from io import StringIO # for Python 3
from Bio import AlignIO
from Bio.Align import AlignInfo
from optparse import OptionParser
from datetime import datetime
from joblib import Parallel, delayed
import multiprocessing
import getpass


startTime = str(datetime.now())
print("Starting time: "+startTime)

# Option Parser: required arguments are two input files (raw barcode-sequence file, and high quality sequence file)
parser = OptionParser()
parser.add_option("-d","--directory", dest="workdir", help="Working directory",default=os.getcwd(),type="string")
parser.add_option("-o","--out", dest="out", help="Output file",default="Seq_barcodes_aligned.txt",type="string")
parser.add_option("--highQual", dest="highQualFile", help="File of barcode-seq association, seq from highest quality read",type="string")
parser.add_option("--inputSeqs", dest="inputSeqsFile", help="Raw barcode, sequence, quality input sequences",type="string")
# Additional options
parser.add_option("-c","--cutoff", dest="cutoff", help="Minimum number of ccs reads for analysis",default=2,type="int")
parser.add_option("-t","--threshold", dest="thresh", help="Minimum threshold to determine consensus sequence",default=0.7,type="float")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_option("-m","--muscle", dest="muscle", help="Compiled MUSCLE program",default="./muscle",type="string")
parser.add_option("-n","--needle", dest="needle", help="Compiled NEEDLE program",default="./needle",type="string")
parser.add_option("--cont", dest="cont", help="Continue working after disrupted run",default=False,action="store_true")
parser.add_option("-s", "--stats", dest="stats", help="Get stats for barcodes that need realignment",default=False,action="store_true")
parser.add_option("-r","--rmint", dest="rm_intermediates",help="Removes intermediate files when finished",default=False,action="store_true")

(options, args) = parser.parse_args()
muscle_exe = options.muscle
needle_exe = options.needle

# **************** Read input files ******************************************************* #
if options.verbose: print("Reading barcodes + reads file...")
# read original assignments into dictionary
hq_dict = {}
assignments = open(options.highQualFile, "r") # min_Q0_assignment.tsv
for line in assignments:
	line = line.upper()
	paired_bcread = line.strip().split()
	hq_dict[paired_bcread[0]] = paired_bcread[1]
assignments.close()
if options.verbose: print("Done reading HQ barcodes.")

# Read all ccs reads into dict: BC: [seq1, seq2...]. Seq quality pairs stored as tuples
if options.verbose: print("Reading all PB reads...")
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
totalBarcodes = len(hq_dict.keys())
if options.verbose: print(str(totalBarcodes) + " barcodes found in hq file")
totalBarcodes2 = len(read_dict.keys())
if options.verbose: print(str(totalBarcodes2) + " barcodes found in other file") 
# ***************************************************************************************** #

os.chdir(options.workdir) # change working directory to output folder

#create intermediates directories in output folder
os.system("mkdir -p intermediates & mkdir -p intermediates/fasta & mkdir -p intermediates/alignments & mkdir -p intermediates/fasta_2 & mkdir -p intermediates/realignments") 

if options.cont:
	outputfile = open(options.out,"a")
else:
	outputfile = open(options.out, "w+")

# **************** if needed to pick up where script broke (option --cont) **************** #
# only run if --cont is included
# deletes progress file if --cont is not included; removes barcodes from hq_dict if --cont is included
# appends to barcode cutoff file if --cont is included; opens & wipes barcode cutoff file if not
progress_file_name = "progress_file.txt" # create progress file that writes barcode each time loop_bcs function is run
under_cutoff_bcs_name = "barcodes_below_cutoff.txt" # barcodes that miss the cutoff
if options.cont:
	if os.path.exists(progress_file_name):
		prog_file = open(progress_file_name,"r")
		remove_keys = [] # list of keys to remove from dict because it's already been processed
		for line in prog_file:
			line = line.strip()
			if line != "":
				remove_keys.append(line)
	# remove barcodes in hq_dict if it's already been processed through loop_bcs function
		for bc in remove_keys:
			if bc in hq_dict:
				hq_dict.pop(bc, None)
		totalBarcodes3 = len(hq_dict.keys())
		new_totalBarcodes = totalBarcodes - totalBarcodes3
		if options.verbose: print(str(new_totalBarcodes) + " barcodes removed from hq file -- already processed.")
	if os.path.exists(under_cutoff_bcs_name):
		cutoff_bcs_file = open(under_cutoff_bcs_name, "a")
else: # if --continue is not included, delete the progress file from the last run.
	if os.path.exists(progress_file_name):
		os.remove(progress_file_name)
		if options.verbose: print("Deleted old progress file")
	cutoff_bcs_file = open(under_cutoff_bcs_name,"w+")
	if options.verbose: print("Created files to track barcodes that don't meet threshold")

# creates or appends to progress file
if not os.path.exists(progress_file_name):
	progress_file = open(progress_file_name,"w+")
else:
	progress_file = open(progress_file_name, "a")
# ***************************************************************************************** #
	
# **************** for analyzing barcodes that don't meet threshold (option -s) *********** #
threshold_analyzed_bcs = []
if not options.cont:
	if os.path.exists("below_threshold_Ncount.txt"):
		os.remove("below_threshold_Ncount.txt")
else:
	if os.path.exists("below_threshold_Ncount.txt"):
		threshold_file_out = open("below_threshold_Ncount.txt","r")
		for line in threshold_file_out:
			line = line.strip().split("\t")
			threshold_analyzed_bcs.append(line[0])
			

def threshold_analysis(_barcode,_consensus):
	threshold_file_out_name = "below_threshold_Ncount.txt"
	if not os.path.exists(threshold_file_out_name):
		threshold_file_out = open(threshold_file_out_name,"w+")
	else:
		threshold_file_out = open(threshold_file_out_name, "a")
	count_Ns = _consensus.count("N")
	threshold_file_out.write(_barcode + "\t" + str(count_Ns) + "\n")
	threshold_file_out.close()
	return(count_Ns)
# ***************************************************************************************** #

# **************** Main alignment/Consensus function ************************************** #
# Align all seuqences with the same barcode. If consensus found, use that sequence. 
# If not, realign to high quality sequence
def loop_bcs(key):
	bc_entry = read_dict[key] #list of sequences
	#create fasta file for each barcode: 
	int_file_name = os.path.join("intermediates/fasta/" + key + ".fasta") 
	if not os.path.isfile(int_file_name): 
		intermediate_file = open(int_file_name, "w+")
		i = 0       
		for item in bc_entry: #add each read for a particular barcode in fasta format
			intermediate_file.write(">" + key + "_" + str(i) + "\n")
			intermediate_file.write(item+"\n")
			i = i+1
		intermediate_file.close()
	if options.verbose: print("made fasta file " + str(key))

	# only align if there are at least CUTOFF ccs reads
	if len(bc_entry) >= options.cutoff:
		if len(bc_entry) == 1: #special case, don't need to align here
			outputfile.write(key+"\t"+bc_entry[0]+"\n")
			#consensus_dict[key] = bc_entry[0]
			
		else: 
			#align files together - first alignment
			aln_file_name = "intermediates/alignments/" + key + ".aln"
			#muscle system call here, write to output file
			muscle_cline = MuscleCommandline(muscle_exe, input=int_file_name, out=aln_file_name)
			stdout, stderr = muscle_cline(int_file_name)
			if options.verbose: print("passed cutoff, made first alignment " + str(key))

			#get consensus: 
			consensus = ""
			if os.path.exists(aln_file_name): #probably should have an else here that throws an error, but it would be better to check that the muscle stderr is empty
				alignment = AlignIO.read(aln_file_name, 'fasta')
				summary_align = AlignInfo.SummaryInfo(alignment)
				consensus = summary_align.gap_consensus(threshold=options.thresh,  ambiguous='N')
				consensus = str(consensus)
				consensus = consensus.replace("-","") 
				if options.verbose: print("got consensus 1 " + str(key))
		
			#if N's: realign (pairwise aligner w/in python) to highest qual, and find consensus from that
			if 'N' in consensus:
				# analyzes barcodes that don't meet threshold, if -s is included
				if options.stats:
					if key not in threshold_analyzed_bcs:
						countN = threshold_analysis(key,consensus)
						if options.verbose: print(key + " barcode has " + str(countN) + " ambiguous sites.")
				#write 1st consensus and HQ read to new file
				int_file_name_2 = os.path.join("intermediates/fasta_2/" + key + ".fasta")
				fasta_2 = open(int_file_name_2,"w+")
				fasta_2.write(">"+key+"\n"+consensus)
				fasta_2.close()
				fasta_hq_name = os.path.join("intermediates/fasta_2/" + key + "_hq.fasta")
				fasta_hq = open(fasta_hq_name,"w+")
				fasta_hq.write(">"+key+"_hq\n"+hq_dict[key])
				fasta_hq.close()
				aln_file_name_2 = "intermediates/realignments/"+key+".aln"
				cmd = needle_exe + " " + int_file_name_2 + " " + fasta_hq_name + " -auto -gapopen 10 -gapextend 0.5 -outfile " + aln_file_name_2 + " -aformat fasta"
				os.system(cmd)

				#consensus of new alignment file (mostly same from previous script)
				alignment_2 = list(SeqIO.parse(aln_file_name_2,"fasta"))
				consensus_seq = str(alignment_2[0].seq)
				hq_seq = str(alignment_2[1].seq)
				finalSeq = ""
				lengthOfAlignment = len(consensus_seq)
				for i in range(lengthOfAlignment):
					if consensus_seq[i] == "N":
						finalSeq = finalSeq+hq_seq[i]
					else:
						finalSeq = finalSeq + consensus_seq[i]
				consensus = finalSeq
				consensus = consensus.replace("-","")		
				outputfile.write(key+"\t"+consensus+"\n")
				#consensus_dict[key] = consensus
				if options.verbose: print("realigned and got new consensus " + str(key))
	
			#if no Ns: write consensus to output file
			else:
				outputfile.write(key+"\t"+consensus+"\n")
				#consensus_dict[key] = consensus
		
			if options.verbose: print("got consensus " + str(key))
		outputfile.flush()
	else: # generates a file of barcodes that did not meet the minimum number of reads (under option -c)
		cutoff_bcs_file.write(key+"\n")
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
	progress_file.write(key+"\n")
	if options.verbose: print("Wrote " + key + " to progress file")
	progress_file.flush()
# ***************************************************************************************** #


# **************** Parallelization and print ********************************************** #
num_cores = multiprocessing.cpu_count()
print("Number of cores: " + str(num_cores))
results = Parallel(n_jobs=(num_cores),prefer="threads")(delayed(loop_bcs)(key) for key in hq_dict)

# removes intermediates folder if --rmint is flagged
if options.rm_intermediates:
	os.system("rm -r intermediates")
else:
	os.system("rm -r intermediates/fasta*")
	
# print all consensus sequences to outfile
#totalCons = len(consensus_dict.keys())
#print("Sequences found for " + str(totalCons) + " barcodes")
#for bc, value in consensus_dict.items():
#	outputfile.write(bc+"\t"+value+"\n")
# close output file  
outputfile.close()
progress_file.close()
cutoff_bcs_file.close()

endTime = str(datetime.now())
print("Ending time: "+endTime)

#remove later?
recordTime = open("run_times.tsv","a+")
recordTime.write(str(getpass.getuser())+"\t"+startTime+"\t"+endTime+"\t"+str(num_cores)+"\n")
recordTime.close()
# ***************************************************************************************** #
