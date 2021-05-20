import gzip, sys

input_file = sys.argv[1] # must be a gzipped fastq
output_file = sys.argv[2] # writes output in fastq.gz format
output_file_handle = gzip.open(output_file, "wb") # write out the output file in gzipped format
passes_list=[] # list of floats holding pass numbers
with gzip.open(input_file, 'r') as fh:
    try:
        fh.read(1)
    except OSError:
        print('input_file is not a valid gzip file by OSError')
    else:
        i=0 #total number of unfiltered reads
        record = []
        for line in fh:
            record.append(line.decode().strip())
            if len(record) == 4:
                i+=1
                tags = record[0].split(";")
                passes=float(tags[6].split("=")[1]) # specific to simlord output formatting
                if passes>=3.0: #only write to output fastq if read has at least 3 passes
                    passes_list.append(passes)
                    readName=tags[0].split("/")[0]
                    variantName=tags[3].split("=")[1]
                    read_info=readName + "/" + variantName + "\n"
                    output_file_handle.write(read_info.encode())
                    output_file_handle.write(record[1].encode())
                    output_file_handle.write("\n".encode())
                    output_file_handle.write(record[2].encode())
                    output_file_handle.write("\n".encode())
                    output_file_handle.write(record[3].encode())
                    output_file_handle.write("\n".encode())
                record = []

        # print stats
        filtered_reads = len(passes_list)
        print(str(filtered_reads)," out of ",str(i), " reads have at least 3 passes")
        average_passes = sum(passes_list) / float(filtered_reads)
        print("Average pass number for filtered reads is: ",str(round(average_passes,2)))

fh.close()
output_file_handle.close()


