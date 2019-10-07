import os
import os.path
import numpy as np
from Bio import SeqIO
import glob
import sys
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
from Bio.Align import AlignInfo

#todo in shell script (write driver script) load muscle/latest
muscle_exe = "/net/gs/vol3/software/modules-repo/RHEL6/muscle"

outputfile = open("test_output.txt", "w")

# TODO parse options using OptionParser 

print("Reading barcodes + reads file...")
# read original assignments into dictionary
bc_dict = {}
assignments = open(sys.argv[1], "r") # min_Q0_assignment.tsv
for line in assignments:
	paired_bcread = line.strip().split()
	bc_dict[paired_bcread[0]] = paired_bcread[1]
assignments.close()

# read all ccs reads into dict: BC: [(seq,qual), (seq,qual)...]
# seq quality pairs stored as tuples
read_dict = {}
reads = open(sys.argv[2], "r") # seq_barcodes.txt
for line in reads: 
	paired_bcread = line.strip().split()
	if paired_bcread[0] in read_dict:
		read_dict[paired_bcread[0]].append((paired_bcread[1],paired_bcread[2]))
	else:
		read_dict[paired_bcread[0]] = [(paired_bcread[1],paired_bcread[2])]
reads.close()
totalBarcodes = len(bc_dict.keys()
print(str(totalBarcodes)) + " barcodes found"

# loop through all barcodes
for key in bc_dict:
	int_file_name = os.path.join("intermediates/fasta/" + key + ".fasta") 
	if not os.path.isfile(int_file_name):	
		intermediate_file = open(int_file_name, "w+")
		i = 0
		for item in read_dict[key]: #add each read for a particular barcode in fasta format
			intermediate_file.write(">" + key + "_" + str(i) + "\n")
			intermediate_file.write(item[0]+"\n")
			i = i+1
		intermediate_file.close()
	# check if at least CUTOFF number of ccs reads here (i >= CUTOFF) default 3?
	
	aln_file_name = "intermediates/alignments/" + key + ".fasta" #.afa suffix?
	#muscle system call here, write to output file
	#TODO: options of which aligner to use
	#shell("clustalo -i {output.intfiles} -o {output.alnfiles}")
	shell("muscle -in "+int_file_name+" -out "+aln_file_name) #string cat before?
	
	#get consensus: 
	alignment = AlignIO.read(aln_file_name, 'fasta')
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = summary_align.dumb_consensus(threshold=0.5,  ambiguous='N') #threshold: default 0.7
	#consensusOutput = consensus.replace("-","") #not sure if there will be gaps in this one
	
	consensusCount = 0
	#if N's: realign (pairwise aligner w/in python) to highest qual
	if 'N' in consensus: 
		#TODO stuff
	
	#if no Ns: write consensus to output file
	else: 
		outputfile.write("key\t"+consensus+"\n")
		consensusCount += 1
	
#print stats on how many had consensus, etc
print(str(consensusCount)+"of "+ str(totalBarcodes)+" barcodes had a consensus sequence")

#close output file
