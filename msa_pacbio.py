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
from optparse import OptionParser
from datetime import datetime


muscle_exe = "/net/gs/vol3/software/modules-sw/muscle/3.8.31/Linux/RHEL6/x86_64/bin/muscle"

# TO RUN SUBPROCESS:
# import subprocess
# cmd = "muscle -in /net/gs/vol1/home/ccyeh/cindy_dunham/SUL1_alleleLibraryCompetition/results/19-08-29_pacbio_run2_np3_muscle_CDS_only/intermediates/consensus_1/TTTTTAAATA.fasta -out test.fasta"
# subprocess.call(cmd,shell=True)
# temp file resources: https://docs.python.org/3/library/tempfile.html#tempfile.mkdtemp

# Option Parser
parser = OptionParser()
parser.add_option("-d","--directory", dest="workdir", help="Working directory",default=os.getcwd(),type="string")
parser.add_option("-o","--out", dest="out", help="Output file",default="Seq_barcodes_aligned.txt",type="string")
parser.add_option("--highQual", dest="highQualFile", help="File of barcode-seq association, seq from highest quality read",type="string")
parser.add_option("--inputSeqs", dest="inputSeqsFile", help="Raw barcode, sequence, quality input sequences",type="string")

#TODO implement these
#parser.add_option("-c","--cutoff", dest="cutoff", help="Minimum number of ccs reads for analysis",default=2,type="int")
#parser.add_option("-t","--threshold", dest="thresh", help="Minimum threshold to determine consensus sequence",default=0.7,type="float")
#parser.add_option("-a","--aligner", dest="aligner", help="Muscle or clustal",default="muscle",type="string") #todo how can i make it from a list?
#parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
(options, args) = parser.parse_args()

os.chdir(options.workdir)
#os.chdir("/net/dunham/vol2/Cindy/pacbio_git/test_pacbio") #this is temporary, just a test

#create intermediates directories
os.system("mkdir -p intermediates & mkdir -p intermediates/fasta & mkdir -p intermediates/alignments & mkdir -p intermediates/fasta_2 & mkdir -p intermediates/realignments") 

outputfile = open(options.out, "w+")

print("Reading barcodes + reads file...")
#print("Time: "+str(datetime.now())
# read original assignments into dictionary
hq_dict = {}
assignments = open(options.highQualFile, "r") # min_Q0_assignment.tsv
for line in assignments:
	paired_bcread = line.strip().split()
	hq_dict[paired_bcread[0]] = paired_bcread[1]
assignments.close()
print("Done reading HQ barcodes.")

# read all ccs reads into dict: BC: [(seq,qual), (seq,qual)...]
# seq quality pairs stored as tuples
print("Reading all PB reads...")
read_dict = {}
reads = open(options.inputSeqsFile, "r") # seq_barcodes.txt
for line in reads: 
	paired_bcread = line.strip().split()
	if paired_bcread[0] in read_dict:
		read_dict[paired_bcread[0]].append((paired_bcread[1],paired_bcread[2]))
	else:
		read_dict[paired_bcread[0]] = [(paired_bcread[1],paired_bcread[2])]
reads.close()
totalBarcodes = len(hq_dict.keys())
print(str(totalBarcodes) + " barcodes found")

# GIANT FOR LOOP HERE
# loop through all barcodes
consensus_dict = {}
consensusCount = 0

for key in hq_dict:
	#create fasta file for each barcode
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
	
	#align files together - first alignment
	aln_file_name = "intermediates/alignments/" + key + ".aln" # aligned file suffix - i think we can make this whatever - it'll be in fasta format
	#muscle system call here, write to output file
	if len(read_dict[key]) > 1:
		muscle_cline = MuscleCommandline(muscle_exe, input=int_file_name, out=aln_file_name)
		stdout, stderr = muscle_cline(int_file_name) #not sure that we need this line at all?
	
	#get consensus: 
	consensus = ""
	if os.path.exists(aln_file_name):
		alignment = AlignIO.read(aln_file_name, 'fasta')
		summary_align = AlignInfo.SummaryInfo(alignment)
		consensus = summary_align.dumb_consensus(threshold=0.5,  ambiguous='N') #threshold: default 0.7
		consensus = str(consensus)
		consensus = consensus.replace("-","") 
	
		#consensusCount = 0 # moving this outside the loop, i think...??
		
	#if N's: realign (pairwise aligner w/in python) to highest qual, and find consensus from that
	if 'N' in consensus:
		#write 1st consensus and HQ read to new file
		int_file_name_2 = os.path.join("intermediates/fasta_2/" + key + ".fasta")
		fasta_2 = open(int_file_name_2,"w+")
		fasta_2.write(">"+key+"\n"+consensus+"\n"+">"+key+"_hq\n"+hq_dict[key]+"\n")
		fasta_2.close()
		aln_file_name_2 = "intermediates/realignments/"+key+".aln"
		muscle_cline_2 = MuscleCommandline(muscle_exe, input=int_file_name_2, out=aln_file_name_2)
		stdout, stderr = muscle_cline_2(int_file_name_2)

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
		consensus_dict[key] = consensus
		
		outputfile.write(key+"\t"+consensus+"\n")
		consensusCount += 1
		
	#if no Ns: write consensus to output file
	else:
		if len(consensus) > 0: # this is only if there were reads to align. maybe this is something that can be optional
			consensus_dict[key] = consensus
			outputfile.write(key+"\t"+consensus+"\n")
			consensusCount += 1

# not sure why this line was necessary?
#	if len(consensus) > 0:
#		print(consensus)
		
#print stats on how many had consensus, etc
print(str(consensusCount)+" of "+ str(totalBarcodes)+" barcodes had a consensus sequence")

#close output file   AGCAGCTGCTGGCTAAGCTAGC
outputfile.close()
#print("Time: "+str(datetime.now())
