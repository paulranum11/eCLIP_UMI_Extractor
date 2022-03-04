#!/usr/bin/env python

#####
#Description reframe the barcode extraction script to do the following.
#1. Improve speed by storing barcodes by position instead of string searching for them.
#2. Support sequence error handling by filtering and correcting barcodes on the fly using hamming distance.  
#   Filtering this way will improve speed by reducing reads carried forward to downstream steps.
#3. Store matepair information.
#4. Reduce the number of output files created. Depricate the single .fastq file per cell output format.
#####

import sys
import os
import argparse
import itertools
import json

##########
# Example use
##########
# python eCLIP_ExtractUMI.py -f path.to.my.fastq.F -r path.to.my.fastq.R -o ./ -b 10000 -v true

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--inputFastqF', required=False, help='Input a forward .fastq format file')
parser.add_argument('-r', '--inputFastqR', required=False, help='Input a reverse .fastq format file')
parser.add_argument('-o', '--outputDir', required=False, help='Name of the output directory used to store the output.fastq')
parser.add_argument('-b', '--bin', required=False, help='Number of reads to process before saving to disc. Binning helps accomodate large input files')
parser.add_argument('-v', '--verbose', required=False, action='store_true', help='Provide -v flag to turn on verbose progress reporting for each bin', default=False)
parser.add_argument('-l', '--linesInInputFastq', required=False, help='provide the length in number of lines of the input fastq. This can be obtained using wc -l myfastq.fastq')
parser.add_argument('-s', '--split_num', required=False, help='provide the number of segments you would like to split the fastq into for binning purposes')
args = parser.parse_args()

# perform fastq losless binning calculation
split_num = int(args.split_num)
fastq_file = str(args.inputFastqF)
lengthFastq = int(args.linesInInputFastq)
    
# Get input fastq file dimensions
print("Getting input fastqF file dimensions")
length_fastq = int(lengthFastq)
print("Fastq length is " + str(length_fastq))
split_size = length_fastq / int(split_num)
print("The split size before tuning is " + str(split_size))
for i in range(100):
    if (split_size % 4 != 0):
        print("Trying split size " + str(split_size))
        split_size = int(split_size) + 1
    if (split_size % 4 != 0):
        print("WARNING!!! Unable to split input fastq without read loss, please try again with a different bin size.")
    print("The split size after tuning is " + str(split_size))

######
# Build a class to store information from each read.
######
print("Creating classes to store FastqReads")
class FastQReadF():
    def __init__(self, name, read, quality, lineNumber, umi):
        self.name = name
        self.read = read
        self.quality = quality
        self.lineNumber = lineNumber
        self.umi = umi

    def display_read(self):
        print("name = " + str(self.name) + "\n" \
            "read = " + str(self.read) + "\n" \
            "quality = " + str(self.quality) + "\n" \
            "lineNumber = " + str(self.lineNumber) + "\n", end='')

    def return_fastq(self):
        return(str(self.name.rstrip()) + "\n" \
            + str(self.read) + "\n" \
            + "+" + "\n" \
            + str(self.quality) + "\n")

    def return_umi(self):
        return(str(self.umi)

class FastQReadR():
    def __init__(self, name, read, quality, lineNumber):
        self.name = name
        self.read = read
        self.quality = quality
        self.lineNumber = lineNumber

    def display_read(self):
        print("name = " + str(self.name) + "\n" \
            "read = " + str(self.read) + "\n" \
            "quality = " + str(self.quality) + "\n" \
            "lineNumber = " + str(self.lineNumber) + "\n", end='')

    def return_fastq(self):
        return(str(self.name.rstrip()) + "\n" \
            + str(self.read) + "\n" \
            + "+" + "\n" \
            + str(self.quality) + "\n")


######
# Create bin parameters
######
# Create a value used to bin the dataset
binIterator = split_size
print("binIterator is set to " + str(binIterator))
# Define workingBin
counter = 0
readCounterFinal = 0
workingBin = counter + binIterator


#######
# Set up a method to simultaneously iterate through forward and reverse fastqs in bins
#######
print("Beginning iteration through forward and reverse fastqs")

# Get the number of lines in the input file
linesInInputFastq = args.linesInInputFastq
print("The linesInInputFastq value is set to " + str(linesInInputFastq))

######
# Step1: Iterate through input fastqs in bins.
######
print("Beginning step 1")
bin_counter = 0

for i in range(0,int(linesInInputFastq),int(binIterator)):
    print("i is currently " + str(i))
    
    #define empty dictionaries (also clears previously filled dictionaries from earlier iteration)
    readsF = {}
    readsR = {}

    
    with open(args.inputFastqF, "r") as fastqF:
        start_ct = 0
        completeReadCounter = 0
        line_ct1 = 0
        read_counter = 0
        starterF = 0
        #print("Processing Forward Reads")
        for line in itertools.islice(fastqF, i, int(i + binIterator)):
            #print("Printing line count " + str(line_ct1))
            if (starterF == 0 and line.startswith('@') == False):
                continue
            if (starterF == 0 and line.startswith('@') == True):
                starterF += 1
                lineName=str(line[0:].rstrip())
                lineNameSplit = lineName.split()
                readIDF = lineNameSplit[0]
                completeReadCounter += 1
                line_ct1 += 1
                continue
            if (line_ct1 % 4 == 0 and line.startswith('@')):
                lineName=str(line[0:].rstrip())
                lineNameSplit = lineName.split()
                readIDF = lineNameSplit[0]
                #print("returning IDF " + readIDF)
                completeReadCounter += 1
            if (line_ct1 % 4 == 1):
                lineReadF=str(line[0:].rstrip())
                completeReadCounter += 2
            if (line_ct1 % 4 == 3):
                lineQuality=str(line[0:].rstrip())
                completeReadCounter += 3
            if (completeReadCounter == 6):
                #print("Read completed appending to dict")
                processedRead = FastQReadF(name = str("@" + lineReadF[0:8] + ":" + lineNameSplit[0] + " " + lineNameSplit[1]), \
                    read = str(lineReadF[8:len(lineReadF)]), \
                    quality = str(lineQuality[8:len(lineQuality)]), \
                    lineNumber = line_ct1, \
                    umi = lineReadF[0:8])
                readsF[readIDF]=processedRead
                del(readIDF)
                completeReadCounter = 0
                read_counter += 1
            if (starterF == 1):
                line_ct1 += 1         



    with open(args.inputFastqR, "r") as fastqR:
        start_ct = 0
        completeReadCounter = 0
        read_counter = 0
        line_ct1 = 0
        starterR = 0
        #print("Processing Reverse Reads")
        for line in itertools.islice(fastqR, i, int(i + binIterator)):
            if (starterR == 0 and line.startswith('@') == False):
                continue
            if (starterR == 0 and line.startswith('@') == True):
                starterR += 1
                lineName=str(line[0:].rstrip())
                lineNameSplit = lineName.split()
                readIDR = lineNameSplit[0]
                completeReadCounter += 1
                line_ct1 += 1
                continue
            if (line_ct1 % 4 == 0 and line.startswith('@')):
                lineName=str(line[0:].rstrip())
                lineNameSplit = lineName.split()
                readIDR = lineNameSplit[0]
                #print("returning IDR " + readIDR)
                completeReadCounter += 1
            if (line_ct1 % 4 == 1):
                lineReadR=str(line[0:].rstrip())
                completeReadCounter += 2
            if (line_ct1 % 4 == 3):
                lineQuality=str(line[0:].rstrip())
                completeReadCounter += 3
            if (completeReadCounter == 6):
                #print("Read completed appending to dict")
                readF_umi = str(readsF[readIDR].return_umi())
                processedRead = FastQReadR(name = str("@" + readF_umi + ":" + lineNameSplit[0] + " " + lineNameSplit[1]), \
                    read = str(lineReadR), \
                    quality = str(lineQuality), \
                    lineNumber = line_ct1)
                readsR[readIDR]=processedRead
                del(readIDR)
                completeReadCounter = 0
                read_counter += 1
            if (starterR == 1):
                line_ct1 += 1



    ######
    # Step2: Write readF and readR to a .fastq file in outputDir directory.
    ######
    #        if not os.path.exists(outputDir):
    #            os.makedirs(outputDir)
    # Write the stored reads to disc by appending to the file "MergedCells_1.fastq"
    print("writing UMI extracted reads to disc")
    input_filenameF = str(args.inputFastqF)
    input_filenameR = str(args.inputFastqR)
    input_filenameF_2 = input_filenameF.split(".")
    input_filenameR_2 = input_filenameR.split(".")
    outfileF = open(str(input_filenameF_2[0] + 'umi_extracted'), 'a')
    outfileR = open(str(input_filenameR_2[0] + 'umi_extracted'), 'a')
    for key in set(readsF.keys()):
        #print("printing the key " + key)
        #readsF[key].display_read()
        #readsR[key].display_read()
        outfileF.write(str(readsF[key].return_fastq()))
        outfileR.write(str(readsR[key].return_fastq()))
        #print(str(readsF[key].return_fastq()))  
    outfileF.close() 
    outfileR.close() 

    bin_counter += 1
    sys.stdout.flush()
