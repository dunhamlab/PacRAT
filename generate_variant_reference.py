import sys
import random

# Input sequence of gene use to generate variant reference (no headers)
inputSeq=open(sys.argv[1], "r").readline().strip()
# Output file name of tab-delimited barcode/variant file, each line contains barcode and associated variant sequence (no plasmid flanking sequences included)
output_var_file = open(sys.argv[2], "w")
# Output file name of fasta-formatted variant file, each fasta entry contains full variant-barcode sequence and flanking plsamid sequence
output_fasta_file = open(sys.argv[3], "w")
# Desired barcode coverage (5 is five barcodes per variant)
covg = int(sys.argv[4])
# Output file name of WT reference, includes flanking sequences, but does not include barcode. Reference used for alignment
wt_ref_file = open(sys.argv[5], "w")

print("Gene lenth: ", len(inputSeq))

# generate n length of reference sequences
def generate_random(len):
   base=["A","T","C","G"]
   random_seq=""
   for i in range(0,len):
      random_base=random.choice(base)
      random_seq+=random_base
   return random_seq

# list of possible codons using NNK oligos
NNK_table={
        'AAG': 'K', 'ATG': 'M', 'ACG': 'T', 'AGG': 'R',
        'TAG': '_', 'TTG': 'L', 'TCG': 'S', 'TGG': 'W',
        'CAG': 'Q', 'CTG': 'L', 'CCG': 'P', 'CGG': 'R',
        'GAG': 'E', 'GTG': 'V', 'GCG': 'A', 'GGG': 'G',
        'AAT': 'N', 'ATT': 'I', 'ACT': 'T', 'AGT': 'S',
        'TAT': 'Y', 'TTT': 'F', 'TCT': 'S', 'TGT': 'C',
        'CAT': 'H', 'CTT': 'L', 'CCT': 'P', 'CGT': 'R',
        'GAT': 'D', 'GTT': 'V', 'GCT': 'A', 'GGT': 'G'}

# generate variants using NNK table
variant_list=[]
tempSeq=""
for i in range(0,len(inputSeq),3):
   codon_left=inputSeq[0:i]
   codon_right=inputSeq[i+3:]
   for j in NNK_table:
      tempSeq=codon_left+j+codon_right
      variant_list.append(tempSeq)

print("Number of variants: ", len(variant_list))

# Use the same sequences for all plasmid flanking sequences
plasmid_seq=['TGGAGTGTTCGGACCAGGAGGAGATTCGTCGAGACTGGCAAGCTGGACGACCTAATCAAGTACGACACCACTACGAGAAAGACACCACGAGCCATTAGG','GCCTATTCGGAGAAATCTATACGCGAAGAGGTCTGTACTACGGATATCAAACCGCTGTTAACAGTTACTATTGATACAGGGACGATTAGTCTGGATGCG','TTCCATTCGTTCACCAAGAACTTCAACCCAAAGTGTACCATGATATACGTAAGTTTGTCCATGTACGGAGGCAAGCTCTCCACTGGGAGTCCTTCTACT']

wildtype_seq=plasmid_seq[0]+inputSeq+plasmid_seq[1]+plasmid_seq[2]
print("Wildtype length (with flanking regions, but without barcode): ",len(wildtype_seq))
wt_ref_file.write(">" + "wildtype" + "\n" + wildtype_seq)
wt_ref_file.close()

# generate barcodes for each variant
barcodeVariantDict = {}
for i in range(0,covg): # generate coverage number of barcodes per variant
   count = 0
   for var in variant_list:
      barcode=generate_random(18)
      while barcode in barcodeVariantDict: #if barcode not unique, generate another
         barcode=generate_random(18)
         print("replicate barcode found")
      ref=plasmid_seq[0]+var+plasmid_seq[1]+barcode+plasmid_seq[2]
      barcodeVariantDict[barcode] = ref
      output_var_file.write(barcode+"\t"+var+"\n")
      output_fasta_file.write(">"+"Variant"+str(count)+"\n")
      output_fasta_file.write(ref+"\n")
      count += 1

print("Number of unique barcodes: ",len(barcodeVariantDict))
output_var_file.close()
output_fasta_file.close()
