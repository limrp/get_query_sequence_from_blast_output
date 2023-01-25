#!/usr/bin/env python3

# get_mlst_from_blastn_results.py

# *--------------------------------------------------------------------------------------------------
# | PROGRAM NAME: get_mlst.py
# | DATE: 14/01/2023 
# | VERSION: 2
# | CREATED BY: Lila Maciel Rodriguez Perez
# | PROJECT FILE: ../my_scripts/get_query_sequence_hit_from_blast_output
# *--------------------------------------------------------------------------------------------------
# | PURPOSE: Extract the query sequence that was the best hit with the MLST sequence of the database
# *--------------------------------------------------------------------------------------------------
# | USAGE:  extract_from_genbank.py -ft 'CDS' -name 'pbp2' -i data/*.gbk -o outputs/pbp2.fasta 
# |         extract_from_genbank.py -ft 'protein' -name 'pbp2' -i data/*.gbk -o outputs/pbp2.fasta 
# |         python get_mlst.py -db mlst.fasta -db_name mlst -i genomes/*.fasta -od res/
# |         python3 ../get_mlst.py -db ../efp_gene_26695_pubmlst.fa -db_name efp_26695 -i ../*.fasta
# *-------------------------------------------------------------------------------------------------
# | WARNING: 
# *-------------------------------------------------------------------------------------------------

# *-------------------------------------  Libraries ------------------------------------------------*
import os
import re
import argparse
import pprint
import gzip
from pprint import pprint

import glob
from pathlib import Path

from Bio import SeqIO # to write sequences in fasta format
from Bio.Blast.Applications import *
from Bio.Blast import NCBIXML
from Bio.Seq import * # funcions like reverse_complement, translate, etc
from Bio.SeqRecord import *  ####====> to get sequence in fasta format

# *------------------------------* Paths to programs in the user's pc *----------------------------* 
my_makeblastdb = "/home/koala/bin/makeblastdb"
my_blastn = "/home/koala/anaconda3/bin/blastn"

# *--------------------------------------* Parsing arguments *-------------------------------------*
# @---- Creating Parser. Adding description
parser = argparse.ArgumentParser()

# @---- Adding arguments
parser.add_argument('-db' , '--db_file', 
                    type = str, 
                    help = 'Fasta file with sequences to build the blast db.'
                    )

parser.add_argument('-db_name', '--database_name',
                    type = str, 
                    help = 'Name or alias of the blast db.'
                    )

parser.add_argument('-i', '--input_query_sequences', 
                    nargs = '+', 
                    type = str, 
                    help = 'Fasta files with sequences to blast against the blast db.'
                    )
## nargs = '+' tells me there will be more than 1 input sequence 
## and all the strings / input sequences will be stored in a list 

parser.add_argument('-o', '--output_filename', 
                    type=str, 
                    help='Name of your final fasta file'
                    )

# @---- Parsing Arguments
args = parser.parse_args()

# *** ______ Storing Parsed Arguments (args) in variables ______ *** #
subject_sequences = args.db_file
my_db_name = args.database_name
genomes_files_list = args.input_query_sequences
out_filename = args.output_filename
blastn_outdir = f"blast/results/blastn/{my_db_name}"


print("\n")
print(f"** Subject sequences to build db: {subject_sequences}")
print(f"** My db name: {my_db_name}")
print(f"** Genome list: {genomes_files_list}")
print(f"** The number of genomes to analyze is {len(genomes_files_list)}")
print(f"** Your final output fasta would be called: {out_filename}")
print("\n")

# *--------------------------------------* Defining functions *------------------------------------*

def make_my_blast_db(path_makeblastdb, sequences_for_db, db_name):
    """Creates a custom Blast Databse

    Args:
        path_makeblastdb (string): Path to program makeblastdb
        sequences_for_db (fasta file): Subject sequences to build the databse
        db_name (string): User defined name of the custom database

    Returns:
        tuple of 2 elements: 
        * the makeblastdb command that was build with the user's input.
        * the makeblastdb executable to actually run the command and generate
            the custom database. 
            Return a tuple of 2 elements: (standard output, standard error).
    """
    
    blast_db_command = NcbimakeblastdbCommandline(cmd = path_makeblastdb,
                                                dbtype = "nucl", 
                                                input_file = sequences_for_db, 
                                                out = f"blast/blast_db/{db_name}")
    # to print the command: blast_db_command
    # to execute the command: blast_db_command()
    # blast_db_command() gives as output another tuple: (std_out, std_err)
    return (blast_db_command, blast_db_command())

def check_and_create_blastDB(db_file_list, path_makeblastdb, sequences_for_db, db_name):
    """_summary_

    Args:
        db_file_list (_type_): _description_
        path_makeblastdb (_type_): _description_
        sequences_for_db (_type_): _description_
        db_name (_type_): _description_
    """
    
    c = 0
    # Check if all 3 db files already exist:
    while c < len(blast_db_files):
        
        file_path = blast_db_files[c]
        
        if Path(file_path).is_file():
            print(f"{file_path} exists!")
            print("\n")
            # if the file exists, keep counting
            c += 1
        else: # If even 1 of the db files doesn't exist:
            print(f"{file_path} doesn't exists! :(")
            print("\n")
            # Create db
            print("Creating Blast DB ...")
            blast_db = make_my_blast_db(
                path_makeblastdb = my_makeblastdb,
                sequences_for_db = subject_sequences,
                db_name = my_db_name
                                )
            # Print the command:
            print("Make blast DB command:")
            print(blast_db[0])
            
            # Execute the command
            # blast_db[1]
            # returns a tuple of (standard output, standard error)
            std_out = (blast_db[1])[0]
            std_err = (blast_db[1])[1]
            print(std_out)
            
            # Now to check if all 3 files as been created
            # reset c to check if now all the files exists
            c = 0
            print("Let's check again!")
            
def get_directory_name(file_path):
    
    try:
        dir_name = Path(file_path).parent.name
        dir_name = (dir_name.split("x"))[0]
        # if this gives an error (if there is no x in dir_name)
        # jump to the except block
    except:
        dir_name = Path(file_path).parent.name
        
    return dir_name

def blastn(path_blastn, output_name, query_file, db_name):
    """BLASTn function to run a query file against a custom blast DB.

    Args:
        path_blastn (string): Path to blastn executable in the users' PC
        output_name (string): Name of the final output
        query_file (string): Name of the query file to run blastn against the custom DB
        db_name (string): Name or alias of the Custom DB created by the user.

    Returns:
        tuple: A tuple with the command ran and the standard output as first and second elements.
    """
    # only for 1 query sequence
    blastn_cl = NcbiblastnCommandline(cmd = path_blastn, 
                                    out = output_name, 
                                    outfmt = 5, 
                                    query = query_file, 
                                    db = f"blast/blast_db/{db_name}", # blast_db already specified here!
                                    #evalue = 0.001
                                    evalue = 0.001
                                    )
    return (blastn_cl, blastn_cl())

def get_genome_name(result_file, db_name):
    """Get genome name from a particularly named file

    Args:
        result_file (string): This files must have the db name and the genome file separated by an _
        db_name (string): Name of the database

    Returns:
        string: The name of the genome that was blasted against the database.
    """
    result_file_name = Path(result_file).name
    
    regex_pattern = my_db_name + r"_(.*)\.xml"
    target_string = str(result_file_name)
    search_result = re.search(regex_pattern, target_string)
    
    return search_result[1]

def get_sequence_from_blast_result_file(result_file, genome, seq_directory, db_name):
    """
    Function to extract the query sequence from the best hit (match) to a sequence
    in the blast database in a blast record.

    Args:
        result_file (string): The blast result file in XML format
        genome (string): The genome name from which the MLST comes from.
        db_name (string): Blast DB named after the subject sequences, in this case, the MLSTs.

    Returns:
        List: Returns a list with a SeqRecord object that contains the MLST found in the query sequence given.
    """
    
    # Reset this value before starting a new results file
    largest_query_length = 0
    
    # List to store the SeqRecord objects (sequences) of a single result file
    seq_record_list = []
    
    with open(result_file, "r") as blast_result_file_handle:
        # Returns a generator object in which we can iterate throught the records
        blast_records = NCBIXML.parse(blast_result_file_handle)
        
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    
                    # HSP information necessary at this point:
                    query_seq_raw = hsp.query
                    query_seq_clean = query_seq_raw.replace("-", "")
                    current_query_length = len(query_seq_clean)
                    
                    if largest_query_length < current_query_length:
                        largest_query_length = current_query_length
                        
                        # Query information necessary at this point:
                        strand = hsp.strand
                        location = f"{hsp.query_start}-{hsp.query_end}"
                        query_len = len(query_seq_clean) # = largest_query_length
                        contig_name = blast_record.query
                        
                        try:
                            # if the contig name follows the pattern
                            # this will work
                            regex_pattern = r"([A-Z0-9]*_[A-Z0-9]*).*"
                            
                            search_result = re.search(
                                pattern = regex_pattern,
                                string = contig_name
                                )
                            
                            contig_name = search_result[1]
                            
                        except:
                            # If the contig name doesn't follow the pattern we will get an error.
                            # After the error, we will jump to this except block
                            # If the contig name doesn follow the pattern, it will keep its old value
                            contig_name = blast_record.query
                        
                        if strand == ('Plus', 'Plus'):
                            # MLST is in Plus sequence
                            # Store sequence as it is
                            # Generate SeqRecord object:
                            
                            seq_record = SeqRecord(
                                seq = Seq(query_seq_clean), 
                                id = f"{genome} | L: {seq_directory} | MLST: {my_db_name} | Contig: {contig_name} |",
                                description = f"Location: {location} | Strand: {strand} | Length: {query_len}"
                                )
                            
                        else: # any other mix of strands
                            # MLST is in Minus (complementary) sequence
                            # Reverse complement the sequence
                            # Before store it in a SeqRecord object:
                            query_seq_rc = reverse_complement(query_seq_clean)
                            #print(query_seq_rc)
                            # Generate SeqRecord object:
                            seq_record = SeqRecord(
                                seq = Seq(query_seq_rc),
                                id = f"{genome} | L: {seq_directory} | MLST: {my_db_name} | Contig: {contig_name} |",
                                description = f"Location: {location} | Strand: {strand} | Length: {query_len}"
                                )
        
        # After getting out of this loop
        # We have finished with this result file
        # the sequence stored in seq_record would be the largest in this result file
        # then, we should append the largest query in this result file
        seq_record_list.append(seq_record)
        
    return seq_record_list

# *-----------------------------------------* Main script *----------------------------------------*

print("\n")
print(
"""
*-------------------------------------------------------------------------------------------------
| PART 0: Creating directories
*-------------------------------------------------------------------------------------------------
""")
print("\n")

# List of directories
dir_list = ["blast/blast_db", 
            "blast/results/blastn", 
            f"blast/results/blastn/{my_db_name}", # blastn_outdir
            "blast/results/sequences"
            ]

# Test if a dir exist first
for dir_path in dir_list:
    if Path(dir_path).exists():
        print(f">> {dir_path} already exist!")
    else:
        # doesn't exist. Create it:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        print(f"@@ {dir_path} was created.")
        


print("\n")
print(
"""
*-------------------------------------------------------------------------------------------------
| PART 1: Create Blast Custom Database
*-------------------------------------------------------------------------------------------------
""")
print("\n")

blast_db_files = [f"blast/blast_db/{my_db_name}.nhr",
                  f"blast/blast_db/{my_db_name}.nin",
                  f"blast/blast_db/{my_db_name}.nsq"]


check_and_create_blastDB(db_file_list = blast_db_files, 
                         path_makeblastdb = my_makeblastdb, 
                         sequences_for_db = subject_sequences, 
                         db_name = my_db_name
                         )

print("\n")
print(
"""
*-------------------------------------------------------------------------------------------------
| PART 2: Blastn of Query sequence vs Blast Custom Database
*-------------------------------------------------------------------------------------------------
""")
print("\n")

counter = 0
# lineage_dict[genome] = directory name = lineage
lineage_dict = {}
blast_results_list = []

for genome in genomes_files_list:
    
    counter += 1
    print(f"File number {counter}: {Path(genome).name} vs {my_db_name}")
    
    genome_stem = Path(genome).stem
    
    # Some genomes starts with GCF and have very long names.
    # I want to save only the part that follows GCF as genome_stem.
    # will genome_stem change or will remain the same?
    # That depends on the file, if it follows the pattern (try) or not (except).
    try:
        print(">> Genomes that start with GCF")
        #genome_core
        #regex_pattern = r"(.*\.[0-9])_.*"
        regex_pattern = r"(GCF_[0-9]*\.[0-9]*)_.*"
        search_result = re.search(pattern = regex_pattern, string = genome_stem)
        print(f">> search_result[1]: {search_result[1]}")
        genome_stem = search_result[1]
        print(f">> Genome Stem: {genome_stem}")
        # if the genome_stem doesn't follow the pattern, it will trown an error
        # we have to catch the error
        # when an error ocurrs, jump to except block
    except:
        # Genomes that don't start with GCF
        print(">> There was an error in the try block, so we jumped to this block")
        #genome_stem = Path(genome).stem # >>>> up <<
        print(f">> genome stem: {genome_stem}")
        
    # Updating lineage dictionary, with genome_stem as keys and lineages (directory name) as values
    lineage_dict[genome_stem] = get_directory_name(file_path = genome)
    
    # Variables for blastn run
    tmp_file_path = f"/tmp/{genome_stem}.fasta"
    blastn_output = f"{blastn_outdir}/{my_db_name}_{genome_stem}.xml"
    # blastn_outdir = f"blast/results/blastn/{my_db_name}"
    
    # Adding result file to list of results:
    blast_results_list.append(blastn_output)
    
    if Path(tmp_file_path).is_file():
        # if tmp_file_path already exist
        # i.e. if the compressed file already has been decompressed in a previous run
        # run blastn using the decompressed file in the tmp directory
        
        print(f"file: {tmp_file_path}")
        print("The file already has been decompressed and it is in /tmp/ directory")
        
        run_blastn = blastn(
            path_blastn = my_blastn,
            output_name = blastn_output,
            query_file = tmp_file_path,
            db_name = my_db_name
            )
        
        # Print the command
        print("** The command line ran was:")
        print(run_blastn[0])
        # Execute / run the command
        run_blastn[1]
        
    # if the file is compressed / has not been decompressed yet
    # To deal with compressed files
    elif Path(genome).suffix == ".gz":
        
        print(f"@@ {genome} is compressed.")
        
        # tmp_file_path = f"/tmp/{genome_stem}.fasta" # genome_tmpfile_path
        
        # Openning 2 files: 1 to read and 1 to write:
        with gzip.open(filename = genome, mode = "rt") as rfh, open(tmp_file_path, "w") as wfh_temp:
            # openning gzipped file with rfh as file handle. 
            # To read its contents into a variable.
            # openning (creating) file genome_outfile_path with temp_wfh as file handle. 
            # To write into it.
            
            print(f"... Unzipping file ...")
            print("... Reading file content into a variable ...")
            upzipped_content = rfh.read()
            wfh_temp.write(upzipped_content)
            print(f"... Writing content to a uncompressed file in {tmp_file_path}.")
            print("... Remember to delete all the fasta files created in /tmp/ !!!!")
        
        
        run_blastn = blastn(
            path_blastn = my_blastn,
            output_name = blastn_output,
            query_file = tmp_file_path,
            db_name = my_db_name
            )
        
        # Print the command
        print("** The command line ran was:")
        print(run_blastn[0])
        # Execute / run the command
        run_blastn[1]
        
        # Adding result file to list of results:
        # blast_results_list.append(f"{blastn_outdir}/{my_db_name}_{genome_stem}.xml")
        
    else:
        # if file is not compressed:
        print(f"The file {genome} is not compressed.")
        
        run_blastn = blastn(
            path_blastn = my_blastn,
            output_name = blastn_output,
            query_file = genome,
            db_name = my_db_name
            )
        
        # Print the command
        print("** The command line ran was:")
        print(run_blastn[0])
        # Execute / run the command
        run_blastn[1]
        
        # Adding result file to list of results:
        # blast_results_list.append(f"{blastn_outdir}/{my_db_name}_{genome_stem}.xml")
        
    print("\n")

# blast_results_list => The output from Part 2 is the input for Part 3
print("Result List:")
pprint(blast_results_list)
print("\n")
print("Lineage dictionary:")
pprint(lineage_dict)

print("\n")
print(f">> Number of genome files analized: {counter}")

print("\n")
print(
"""
*-------------------------------------------------------------------------------------------------
| PART 3: Parse Blast results files and extract Query sequence Best Hit with DB
*-------------------------------------------------------------------------------------------------
""")
print("\n")

#blast_results_files = Path(blastn_outdir).glob("*.xml")

b = 0
sequences_list = []

# Based on query length, not on HSP alignment length
for blast_result_file in blast_results_list:
    
    # To count blastn results xml files
    b += 1
    #print(blast_result_file)
    
    # Get only the name of the result file
    result_file_name = Path(blast_result_file).name
    print(f"*** BLASTn result file number {b}: {result_file_name}")
    #print("\n")
    
    # Get the genome name from xml blastn result file
    genome_name = get_genome_name(
        result_file = blast_result_file,
        db_name = my_db_name
        )
    
    print(f"genome name: {genome_name}")
    print(f"Lineage: {lineage_dict[genome_name]}")
    
    # Parse Blastn results
    seq_record_list = get_sequence_from_blast_result_file(
        result_file = blast_result_file, 
        genome = genome_name, 
        seq_directory = lineage_dict[genome_name],
        db_name = my_db_name
        )
    
    print(f"Number of seq records in this list: {len(seq_record_list)}")
    
    # Was the sequence found or not in this result file?
    if len(seq_record_list) < 1:
        print(f"@ {my_db_name} was NOT FOUND in {genome_name}")
    else:
        print(f"@ {my_db_name} was found in {genome_name}")
        print(f"@ Number of desired sequences found in {genome_name}: {len(seq_record_list)}")
        sequences_list.extend(seq_record_list)
        print("\n")
        
    # to the next blast result file (xml) in the directory ==>

print("\n")
print(f">> The number of BLASTn results files (.xml) analized was {b}.")
print(f">> Number of {my_db_name} sequences found in all {b} files: {len(sequences_list)}")
print("\n")

# Writing sequences to a fasta file

final_sequence_file = f"./blast/results/sequences/{out_filename}"

SeqIO.write(sequences = sequences_list, 
            handle = final_sequence_file, 
            format = "fasta"
            )

print("If you want to keep only the first word in your sequences' ID, run:")
print("$ sed -re 's/^(>.*) \| L:.*/\1/' your_sequence.fa > your_sequence_shortIDs.fa")
print("If you want to sort your sequence in natural order (N), run:")
print("$ seqkit sort -N -i -2 my_seqs.fasta -o my_seqs.sorted.fasta")
print("\n")
print("Have fun with your sequences!!")
print("\n")