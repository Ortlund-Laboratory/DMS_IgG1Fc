# But our barcodes are at the R2 end.
# We therefore need to convert R2 sequences to R1 via reverse transcription
# Alasdair Keith script (if any problems, email adkeith@emory.edu)
#################################################################

import gzip
import csv

# Open barcode_runs.csv to identify which fastq.gz files we need to convert
# These are found in last column.
with open('./data/barcode_runs.csv','r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader)
    for row in reader:
        last_col = row[-1]

# Open fastq.gz files, maintaining compression to avoid memory problems
        with gzip.open(last_col,'rb+') as file:
# counter is used to determine what to do with each of the four lines
# of a fastq file.
            counter=3
            for line in file:
                counter+=1
                if (counter % 4 == 0):
# Need to relabel sequences from 2 to 1, since we're now looking at R1, not R2
                    write1 = line.replace(b'2:N:0', b'1:N:0')
                    write1 = write1.strip()
                elif (counter % 4 == 1):
# Reverse transcription of barcodes occurs in this branch
                    reversedstring=reversed(line)
                    rs=bytes(reversedstring)
                    write2 = rs.replace(b'A',b'X').replace(b'T',b'A').replace(b'X',b'T').replace(b'C',b'Y').replace(b'G',b'C').replace(b'Y',b'G')
                    write3 = b'\n'
                elif (counter % 4 == 2):
                    write4 = line
                    write4 = write4.strip()
                elif (counter % 4 == 3):
# Need to reverse accuracies of reads as well
                    reversedacc=reversed(line)
                    ra=bytes(reversedacc)
                    write5 = ra
                    write6 = b'\n'
                    last_col_file = last_col.split('/')
# Write to new file, prefixed with rt_
                    with gzip.open('rt_'+last_col_file[-1],'ab+') as outfile:
                        outfile.write(write1)
                        outfile.write(write2)
                        outfile.write(write3)
                        outfile.write(write4)
                        outfile.write(write5)
                        outfile.write(write6)
                    outfile.close()
                else:
                    break

csvfile.close()

# Now we need to generate a new barcode_runs file
# prefixed with rt_, which includes all the fastq.gz
# files from before, only now they are all prefixed with rt_ also
# ie this is the list of our new reverse transcribed sequences.

with open('./data/barcode_runs.csv','r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rowcounter=0
    for row in reader:
        if (rowcounter==0):
            with open('./data/rt_barcode_runs.csv','w') as outcsvfile:
                outcsvfile.write(','.join(row))
                outcsvfile.write('\n')
            outcsvfile.close()
        else:
            last_col = row[-1]
            last_col_file = last_col.split('/')
            newfilefolder = last_col_file[0:-1]
            newfilename = 'rt_'+last_col_file[-1]
            newline = row[0:7]
            with open('./data/rt_barcode_runs.csv','a') as outcsvfile:
                outcsvfile.write(','.join(newline))
                outcsvfile.write(',')
                outcsvfile.write('/'.join(newfilefolder))
                outcsvfile.write('/')
                outcsvfile.write(newfilename)
                outcsvfile.write('\n')
        rowcounter+=1
