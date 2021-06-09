import sys
import gzip

# Read in two bc-variant tsv files to compare, can be gzipped
# Read each file into a dictionary
file1 = sys.argv[1]
file1_dict = {}
if file1[-3:] == ".gz":
    fh1 = gzip.open(sys.argv[1],"r")
    for line in fh1:
        info = line.decode().strip().split()
        barcode = info[0]
        if barcode in file1_dict:
            file1_dict[barcode]+=info[1]
            print("duplicate barcode found: ", barcode)
        else :
            file1_dict[barcode]=info[1]
else:
    fh1 = open(sys.argv[1], "r")
    for line in fh1:
        info = line.strip().split()
        barcode = info[0]
        if barcode in file1_dict:
            file1_dict[barcode]+=info[1]
            print("duplicate barcode found: ", barcode)
        else :
            file1_dict[barcode]=info[1]

file2 = sys.argv[2]
file2_dict = {}
if file2[-3:] == ".gz":
    fh2 = gzip.open(sys.argv[2], "r")
    for line in fh2:
        info = line.decode().strip().split()
        barcode = info[0]
        if barcode in file2_dict:
            file2_dict[barcode]+=info[1]
            print("duplicate barcode found: ", barcode)
        else :
            file2_dict[barcode]=info[1]
else:
    fh2 = open(sys.argv[2], "r")
    for line in fh2:
        info = line.strip().split()
        barcode = info[0]
        if barcode in file2_dict:
            file2_dict[barcode]+=info[1]
            print("duplicate barcode found: ", barcode)
        else :
            file2_dict[barcode]=info[1]

# Print stats
print("Number of barcodes in file 1: ", str(len(file1_dict)))
print("Number of barcodes in file 2: ", str(len(file2_dict)))

# Function to compare two bc-variant dictionaries and find the number of indels
def compare_maps(dict1, dict2):
    match = 0
    diff = 0
    indel = 0
    for barcode in dict1:
        if barcode in dict2: ## barcodes are in both files
            barcode1 = dict1[barcode]
            barcode2 = dict2[barcode]
            if barcode1 == barcode2: # compare
                match +=1
            else:
                diff += 1
                if len(barcode1) != len(barcode2): indel += 1

    return match, diff, indel

matched,different,indels= compare_maps(file1_dict, file2_dict)

# Print stats
print("Number of barcodes with matching sequences: ", str(matched))
print("Number of barcodes with different sequences: "+ str(different)+" ("+str(indels)+" barcodes with indels)")
