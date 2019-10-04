import os
import os.path
import numpy as np
from Bio import SeqIO
import glob
import sys

print("Reading barcodes + reads file...")
bc_dict = {}
reads = open(sys.argv[1], "r")
for line in reads: #make dict with barcode as key and list of reads as the value
	paired_bcread = line.split()
	if paired_bcread[0] in bc_dict:
		bc_dict[paired_bcread[0]].append(paired_bcread[1])
	else:
		bc_dict[paired_bcread[0]] = []
		bc_dict[paired_bcread[0]].append(paired_bcread[1])
reads.close()

for key in bc_dict:
	int_file_name = os.path.join("intermediates/fasta/" + key + ".fasta") #make a new file in /intermediates/fasta called {barcode}.fasta
	if not os.path.isfile(int_file_name):	
		intermediate_file = open(int_file_name, "w+")
		i = 0
		for item in bc_dict[key]: #add each read for a particular barcode in fasta format
			intermediate_file.write(">" + key + "_" + str(i) + "\n")
			intermediate_file.write(item+"\n")
			i = i+1
		intermediate_file.close()
print("Done reading file and making intermediate fastas.")
print(str(len(bc_dict.keys())) + " barcodes found")

# BARCODES=list(set(bc_dict.keys()))
# print("got barcodes")