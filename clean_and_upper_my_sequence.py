#!/usr/bin/env python3

import os
import re
import argparse
import gzip

import glob
from pathlib import Path

from Bio import SeqIO # to write sequences in fasta format
from Bio.Blast.Applications import *
from Bio.Blast import NCBIXML
from Bio.Seq import * # funcions like reverse_complement, translate, etc
from Bio.SeqRecord import *  ####====> to get sequence in fasta format

# *--------------------------------------- Parsing arguments --------------------------------------
# creating parser. 
parser = argparse.ArgumentParser()

# Adding arguments
## input sequence will be stored as a list by using nargs='+'
parser.add_argument(
    '-i', '--input',
    nargs='+', # 1 or more
    type=str,
    help='Fasta files with sequences that have gaps (-) and are in lower case.'
    )

# output directory
parser.add_argument(
    '-o', '--output_directory',
    type=str,
    help="Name of the directory where all the output files will be stored. Don't add the last / to the path"
    )

# Parsing arguments
args = parser.parse_args()

# *** ______ Storing arguments in variables ______ ***
input_files_list = args.input
output_dir = args.output_directory

# *---------------------------------------* Main script *---------------------------------------*

counter = 0

for input_file in input_files_list:
    
    counter += 1
    print("\n")
    print(f"File number {counter}: {Path(input_file).name}")
    print("\n")
    
    # This list has to be resetted to empty before starting a new file
    seq_record_list = []
    
    with open(file = input_file, mode = "r") as fh:
        # Read multifasta file with parse function
        records = SeqIO.parse(handle = fh, format = "fasta")
        
        for record in records:
            my_dna = record.seq
            my_dna_corrected = my_dna.ungap("-").upper()
            
            dna_length = len(my_dna_corrected)
            
            seq_record = SeqRecord(seq = my_dna_corrected, 
                                id = record.id, 
                                description = f"Length: {str(dna_length)}")
            
            seq_record_list.append(seq_record)
            
        # After finishing iterating trought all the records in the current file
    # The file is already closed.
    
    # Saving the list of sequences to a fasta file
    output_file = f"{output_dir}/{Path(input_file).stem}_cleaned-upper.fasta"
    
    SeqIO.write(sequences = seq_record_list, 
                handle = output_file, 
                format = "fasta")
    
    # To the next input file in the list ==>


print(f"The number of files processed was: {counter}")
print(f"Your files can be found in: {output_dir}")
print("Your files have been processed and your sequences are in uppercase and clean now!")
print("\n")
