# this script retrieves the consensus of the aligned barcode sequences, and then puts it into a fasta file with the highest quality read
import os
import os.path
import numpy as np
from Bio import SeqIO
import glob
import sys

alnFile = open("intermediates/alignments/"+sys.argv[1]+".aln", "r")
lengthFile = alnFile.read()
alnFile.seek(0)
if len(lengthFile) > 0:
	indiv_reads = list(SeqIO.parse(alnFile, "fasta")) #parse file in seq and ID
	seqLength = len(indiv_reads[0].seq) #what is the sequence length?
	consensusDict = {}
	for k in range(seqLength-1): #each key goes from 0 to seqLength-1. each key corresponds to an empty list
		consensusDict[k] = []
	for reads in indiv_reads: # for each individual read of the same barcode
		seq = reads.seq # get seq
		name = reads.id # and the name
		for i in range(len(str(seq))-1): #for each letter position in the read
			consensusDict[i].append(seq[i]) #add letter to dictionary by position in string
	consensusSeqList = []
	nucleotides = ["A","G", "C", "T", "-"]
	for n in range(len(consensusDict)): #for each item in dictionary (all letters at a position), pick the most frequent one
		countDict = {"A":0, "G":0, "C":0, "T":0, "-":0}
		letterAtPosition = ""
		for p in consensusDict[n]:
			countDict[p] = countDict[p] + 1
		maxVal = 0
		whichNtMax = []
		for nt in nucleotides:
			if countDict[nt] > maxVal:
				whichNtMax = []
				maxVal = countDict[nt]
				whichNtMax.append(nt)
			elif countDict[nt] == maxVal:
				whichNtMax.append(nt)
		if len(whichNtMax) == 1:
			letterAtPosition = whichNtMax[0]
		elif len(whichNtMax) > 1: # if there's ambiguity, make that nt at the position "N"
			letterAtPosition = "N"
		consensusSeqList.append(letterAtPosition)
	outputSeq = ''.join(consensusSeqList) # final sequence
	outputSeq = outputSeq.replace("-","") # remove gaps
	alnFile.close()
	subFile = open("intermediates/consensus_1/"+sys.argv[1]+".fasta", "w+")
	name = name.split("_")
	name = name[0]
	subFile.write(">"+name+"\n"+outputSeq+"\n") #write sequence to fasta file
	if "N" in outputSeq: #if "N" is in the file, add the highest qual read for that barcode to same fasta file
		highQualRead = open(sys.argv[2],"r")
		for line in highQualRead:
			line = line.strip().split()
			if line[0] == name:
				bc = line[0]
				highQualRead = line[1]
				subFile.write(">"+name+"_highqual\n"+highQualRead+"\n")
		subFile.close()
	else: # if it has no "Ns," add the same fasta file to the final_consensus folder. 
		fileConsensusName = os.path.join("intermediates/final_consensus/"+name+".fasta")
		fileConsensus = open(fileConsensusName, "w+")
		fileConsensus.write(">"+name+"\n"+outputSeq+"\n")
		fileConsensus.close()
		#not sure if this is needed
# 		fileConsensusName = os.path.join("intermediates/consensus_2/"+name+".fasta")
# 		fileConsensus = open(fileConsensusName, "w+")
# 		fileConsensus.write(">"+name+"\n"+outputSeq+"\n")
# 		fileConsensus.close()
else:
	open("intermediates/consensus_1/"+sys.argv[1]+".fasta", "a+").close()
alnFile.close()
