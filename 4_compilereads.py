import os
import os.path
import numpy as np
from Bio import SeqIO
import glob
import sys

# this script creates a file of the final barcode:variant map.

barcode_file = open(sys.argv[1], "r")

subassembled_file_name = "intermediates/bc_final_map.fasta"
subassembled_file_tsv_name = "intermediates/bc_final_map.tsv"
if not os.path.isfile(subassembled_file_name) and not os.path.isfile(subassembled_file_tsv_name):
	subassembled_file = open(subassembled_file_name,"w+")
	subassembled_file_tsv = open(subassembled_file_tsv_name,"w+")
	unique_barcodes = []
	for line in barcode_file:
		line = line.split()
		unique_barcodes.append(line[0])
	barcode_file.close()
	unique_barcodes=list(set(unique_barcodes))
	print("Done  getting barcodes. " + str(len(unique_barcodes)) + " unique barcodes found.")
	
	i=0
	for barcode in unique_barcodes:
		bc_filename = "intermediates/final_consensus/" + barcode + ".fasta"
		bc_file=open(bc_filename,"r")
		for line in bc_file:
			line = line.strip()
			if line[0] == ">":
				line = line.split("_")
				subassembled_file.write(line[0]+"\n")
				subassembled_file_tsv.write(line[0][1:]+"\t")
			else:
				subassembled_file.write(line+"\n")
				subassembled_file_tsv.write(line+"\n")
		bc_file.close()
		i += 1
	print("Number of barcodes: " + str(i))
	
	subassembled_file.close()
	subassembled_file_tsv.close()
	
#********************* reverse complement ********************* #
subassembled_file = open(subassembled_file_name, "r")
read_seqs = list(SeqIO.parse(subassembled_file, "fasta"))
subassembled_file_rc = open("intermediates/bc_final_map_revcomp.fasta", "w+")
bc_map_final = open("intermediates/bc_final_map_revcomp.tsv","w+")
# aa sequence from reverse complemented sequence
subassembled_file_aa = open("./intermediates/bc_map_aa.fasta", "w+")
unique_allele = []
for sequence in read_seqs:
	unique_allele.append(str(sequence.seq))
	raw_seq = sequence.seq
	raw_seq = raw_seq.reverse_complement()
# 	if len(raw_seq) % 3 == 1:
# 		raw_seq = raw_seq + "NN"
# 	elif len(raw_seq) % 3 == 2:
# 		raw_seq = raw_seq + "N"
	if not raw_seq.startswith("A"):
		raw_seq = "A" + raw_seq
	subassembled_file_rc.write(">"+str(sequence.id)+"\n"+str(raw_seq)+"\n")
	bc_map_final.write(str(sequence.id)+"\t"+str(raw_seq)+"\n")
	aa_seq = raw_seq.translate()
	subassembled_file_aa.write(">"+str(sequence.id)+"\n"+str(aa_seq)+"\n")



