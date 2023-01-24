#!/usr/bin/env bash

# *------------------------------------------------------------------------------------------
# | PROGRAM NAME: concatenate_msa_alingments.sh
# | DATE: 16/02/22 
# | CREATED BY: Lila Maciel Rodriguez Perez 
# | PROJECT FILE: ~/biotools/my_scripts/phyloscripts 
# *------------------------------------------------------------------------------------------
# | PURPOSE: 
# | 1. Do MSA of each gene or protein 
# | 2. Concatenate the individual Multiple Sequence Alignments.
# | 3. Trimm individual MSAs if necessary (optional) 
# | (when they have different lengths)
# *------------------------------------------------------------------------------------------
# | USAGE: 
# *------------------------------------------------------------------------------------------
# | WARNING: 
# *------------------------------------------------------------------------------------------
echo
echo "Started at `date +%T`"
echo

echo """
# *----------------------------------------------------------------------------
# | PART 1: Looking for programs and variables, creating directories ...
# *----------------------------------------------------------------------------
"""
echo

# Programs
catfasta2phyml="/home/koala/biotools/catfasta2phyml/catfasta2phyml-v.1.2.0/catfasta2phyml.pl"
perl "$catfasta2phyml" -V

# Directories
mkdir -p msa/{edited_multifasta,individual_msa,supermatrix}

#out_dir1="msa/edited_multifasta"
out_dir2="msa/individual_msa"
out_dir3="msa/supermatrix"

#suffix="$1"

# 1. Edit sequence names inside each multifasta file
## Only keep gene or protein name


# 2. Multiple Sequence Alignment of individual genes or proteins

echo """
# *----------------------------------------------------------------------------
# | PART 2: Multiple Sequence Alignment of individual genes or proteins
# *----------------------------------------------------------------------------
"""
echo

#for file in *.{fasta,fa};
for file in *.fasta;
do
	#mafft -auto "$out_dir1"/"$file" > "$out_dir2"/"$file".msa.fasta
	#base=$(basename --suffix=".fasta" "$file")
	base=$(echo "$file" | sed -re 's/^(.*)_complete_.*/\1/')
	echo ">>>> Processing file of MLST: $base"
	echo "Doing MSA of $file"
	echo
	mafft --auto "$file" > "$out_dir2"/"$base".msa.fasta
done

# 3. Concatenate all the individual aligments

echo """
# *----------------------------------------------------------------------------
# | PART 3: Concatenate all the individual aligments
# *----------------------------------------------------------------------------
"""
echo

out_supermatrix="mlst.concatenated.msa.fasta"
out_partitions="mlst.partitions.txt"

echo "$out_supermatrix" 
echo "$out_partitions"
echo

# Maybe another way of running this program:
perl "$catfasta2phyml" --fasta \
	"$out_dir2"/*.msa.fasta \
	1> "$out_dir3"/"$out_supermatrix" 2> "$out_dir3"/"$out_partitions"

echo
#ls msa/*/ | egrep "*supermatrix*"
#ls msa/*/ | egrep "*partitions*"
   
echo

echo "Finished at `date +%T`"
echo
echo "Have fun with your super MSA!!!"
echo

