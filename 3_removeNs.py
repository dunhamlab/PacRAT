# this file removes any residual ambiguity. which hopefully shouldn't be present at this point
import os
import os.path
import numpy as np
from Bio import SeqIO
import glob
import sys

alignedFile = open("intermediates/alignments/"+sys.argv[1]+".aln", "r")
lengthFile = alignedFile.read()
line1 = alignedFile.readline()
alignedFile.seek(0)
if len(lengthFile) > 0 and ">" in line1:
	reads = list(SeqIO.parse(alignedFile, "fasta"))
	alignedFile.close()
	finalSeq = ""
	lengthOfAlignment = len(str(reads[0].seq))
	consensusSeq = str(reads[0].seq)
	highQualSeq = str(reads[1].seq)
	for i in range(lengthOfAlignment):
		if consensusSeq[i] == "N":
			finalSeq = finalSeq + highQualSeq[i]
		else:
			finalSeq = finalSeq+consensusSeq[i]
	finalSeq = finalSeq.replace("-","")
	consensusFile=open("intermediates/final_consensus/"+sys.argv[1]+".fasta", "w+")
	name = str(reads[0].id)
	consensusFile.write(">"+name+"\n"+finalSeq)
	consensusFile.close()
else:
	alignedFile.close()
	if os.path.exists("intermediates/final_consensus/"+sys.argv[1]+".fasta"):
		os.utime("intermediates/final_consensus/"+sys.argv[1]+".fasta", None)
	else:
		open("intermediates/final_consensus/"+sys.argv[1]+".fasta","a+").close()
	#subassembledFile = open(output.final_subassembly, "a+")
