#!/bin/bash

#### POTYVIRUS ASSEMBLER ####

#This script takes paired-end illumina RNA-seq reads as input (in fastq format), uses BLAST to find the reads that align with known
#reference potyviral sequences, and then uses seqtk to substract the potyviral reads from the original fastq files. Using the substracted
#reads then it uses rnaSPAdes to assemble them. Finaly, the resulted viral contigs are annotated and are provided as output in a new file

### USAGE: bash potyvirus_assembler.sh <right reads in fastq format> <left reads in fastq format>
#   Example 1: bash potyvirus_assembler.sh Data/P_edulis_1.fastq Data/P_edulis_2.fastq
#   Example 2: bash potyvirus_assembler.sh reads1.fq reads2.fq

### INPUTS
# 2 fastq files, right and left reads, provided in arguments 1 and 2 of the command line respectively (see usage)

### MAIN OUTPUTS
# -'annotated_potyvirus_contigs.fasta' in the 'Results' folder with the the annotated viral contigs with the identity of the reference
# sequence to which they align, the contig and reference length, the alignment length and the percentage of identity\n"
# -'summary.txt' in the 'Results' folder containing a summarized table with the information for each viral contig

### OTHER OUTPUTS
# rnaSPAdes folder with all its respective files
# fastq files (right and left) in the 'Data' folder with only the viral reads
# 'Blast_results' folder with all the results from the BLASTS in the script
# 'DataBase' folder containing the fasta file with the reference potyviral sequences and the BLAST database files

### REQUIREMENTS
# seqtk (https://github.com/lh3/seqtk)
# local BLAST (https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/)
# rnaSPAdes (http://cab.spbu.ru/software/spades/)
# bash (Unix Shell) (Part of MacOS and Linux operating systems)

# NOTE: The script was written using a MacOS Mojave version 10.14.3 operating system. It should work without issues in any
# other MacOS operating systems. It should also work in Linux computers, although some minor modifications for commands such as 
# 'curl', 'sed', 'grep', 'awk' might be needed. 

#########################################################  SCRIPT #########################################################

#Usage function that prints the correct usage of the script  
display_usage() { 
	echo -e "USAGE: bash potyvirus_assembler.sh <right reads in fastq format> <left reads in fastq format>\nExample: bash potyvirus_assembler.sh reads1.fq reads2.fq\n" 
	} 

# If less than two arguments are given then display usage 
	if [  $# -le 1 ] 
	then 
		display_usage
		exit 1
	fi 

# If more than two arguments are given then display usage 
	if [  $# -ge 3 ] 
	then 
		display_usage
		exit 1
	fi 

#Indicates that file1 and file2 are variables that store the file name and path given in the command 
#line arguments 1 and 2 respectively
file1=$1
file2=$2

#Creates folders needed, in case they are not already present
if [ ! -d DataBase ]; then 
    mkdir DataBase
fi

if [ ! -d Temp  ]; then 
    mkdir Temp
fi

if [ ! -d Blast_results  ]; then 
    mkdir Blast_results
fi

if [ ! -d Data ]; then 
    mkdir Data
fi

if [ ! -d Results ]; then 
    mkdir Results
fi

#Prints a statement in the terminal indicating the files given by the user for the script
echo "Forward reads file given: $1"
echo "Reverse reads file given: $2"

#Step 1 donwloads from a Google Drive folder a file with the reference genome sequences of Potyvirus available in NCBI
echo "STEP 1: DOWNLOADING POTYVIRUS REFERENCE GENOMES DATABASE"
cd DataBase
curl -L -o Potyvirus.fasta "https://drive.google.com/uc?export=download&id=1oSLRLKN2un_KfpW2801k3VfVu9AEgvMo"
echo "SUCCESSFULLY DOWNLOADED DATABASE!"

#Step 2 creates with the sequences downloaded a BLAST database 
echo "STEP 2: CREATING A NUCLEOTIDE BLAST DATABASE"
makeblastdb -in Potyvirus.fasta -input_type fasta -dbtype nucl -out Potyvirus_DB > /dev/null 2>&1
echo "SUCCESSFULLY MADE DATABASE!"

cd ..

#Step 3 uses seqtk to convert one of the fastq files into fasta to be able to use it in BLAST
echo "STEP 3: ANALYZING READS AND SUBSTRACTING ALL THE POTYVIRAL READS (MIGHT TAKE A WHILE)"
seqtk seq -a "$1" > Temp/data_1.fasta

#Then, performs a blastn using the potyvirus database and the reads given by the user with an strict evalue to minimize 
#the presence of non-viral reads. The output is just the query (read) name in a list that will be used in the next command
blastn -query Temp/data_1.fasta -db DataBase/Potyvirus_db -outfmt "6 qseqid" -evalue 1e-15 -max_target_seqs 1 -out Blast_results/potyvirus_reads.txt > /dev/null 2>&1

#Now it uses seqtk to substract from the original .fastq files only the viral reads given by the list outputed by BLAST in the previous
#command
seqtk subseq "$1" Blast_results/potyvirus_reads.txt > Data/Potyvirus_1.fastq
seqtk subseq "$2" Blast_results/potyvirus_reads.txt > Data/Potyvirus_2.fastq

#Counts the total reads and the viral reads and prints a statement in the terminal indicating the number of each, and the percentage
#of viral reads in the sample
viral_reads=$(wc -l < Blast_results/potyvirus_reads.txt)
total_reads=$(wc -l < "$1")
total_reads_integer=$(echo $(( $total_reads/4)))
percentage_reads=$(echo $(($viral_reads*100)))
percentage_reads_2=$(echo $(($percentage_reads/$total_reads_integer)))

echo "The number of total reads is = $total_reads_integer"
echo "The number of viral reads is = $viral_reads ($percentage_reads_2 % of total reads)"

#Step 4 uses the viral-reads-only files obtained with seqtk as input for rnaSPAdes which will assemble them into contigs
echo "STEP 4: RUNNING rnaSPAdes TO ASSEMBLE VIRAL READS (MIGHT TAKE A WHILE)"
rnaspades.py -o rnaSPAdes_potyvirus -1 Data/Potyvirus_1.fastq -2 Data/Potyvirus_2.fastq > /dev/null 2>&1

echo "SUCCESSFULLY ASSEMBLED READS!"

#Takes the hard_filtered_transcripts from the rnaSPAdes output and BLASTs them against the potyviral database to select only the viral 
#contigs (as sometimes plant reads can still be present in low percentages and result in a contig assembly). This BLAST outputs in a 
#table the contig name, reference name, contig length, reference length, length of alignment and percentage of identity
echo "STEP 5: VERIFYING IF THE CONTIGS ASSEMBLED AND THEIR IDENTITY"
blastn -query rnaSPAdes_potyvirus/hard_filtered_transcripts.fasta -max_target_seqs 1 -db DataBase/Potyvirus_DB -outfmt "6 qseqid stitle qlen slen length pident" -evalue 1e-5 -out Blast_results/Blastn_after_assembly.txt > /dev/null 2>&1

#This BLAST outputs a list with only the names of the viral contigs to use in the next command
blastn -query rnaSPAdes_potyvirus/hard_filtered_transcripts.fasta -max_target_seqs 1 -db DataBase/Potyvirus_DB -outfmt "6 qseqid" -evalue 1e-5 -out Blast_results/Blastn_after_assembly_contig_names.txt > /dev/null 2>&1

#Using the list with the viral contigs from the previous command, seqtk substracts from the 'hard_filtered_transcripts.fasta', only 
#the viral contigs
seqtk subseq rnaSPAdes_potyvirus/hard_filtered_transcripts.fasta Blast_results/Blastn_after_assembly_contig_names.txt > Temp/viral_contigs.fasta

#Replaces spaces and commas in the BLAST results table to be able to manipulate it more easily
sed -e "s/,//g" -e "s/ /_/g" Blast_results/Blastn_after_assembly.txt > Blast_results/Blastn_after_assembly_formated.txt 

#Takes the BLAST table and outputs a summary table indicating for each viral contig the contig name, reference name, contig length, reference length, length of alignment and percentage of identity in a more "readable" manner
awk '{print $1" REFERENCE_SEQUENCE_WITH_BEST_HIT="$2",CONTIG_LENGTH="$3",REFERENCE_LENGTH="$4",ALIGNMENT_LENGTH="$5",PERCENT_IDENTITY="$6}' Blast_results/Blastn_after_assembly_formated.txt > Results/summary.txt

#Takes the summary table and the viral_contigs file and outputs an annotated fasta file where each sequence title is replaced by the
#information with the reference viral sequence to which the contig aligns to, the ontig length, reference length, length of alignment
#and percentage of identity along with the sequence of each contig. 
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{sub(/>/,">"a[$1]",ORIGINAL_CONTIG_NAME=",$1)}1' Results/summary.txt Temp/viral_contigs.fasta > Results/annotated_potyvirus_contigs.fasta

#Counts the number of total contigs produced by rnaSPAdes and the number of those contigs that are indeed viral contigs
number_contigs=$(grep -c ">" rnaSPAdes_potyvirus/hard_filtered_transcripts.fasta)
number_viral_contigs=$(grep -c ">" Results/annotated_potyvirus_contigs.fasta)

#Prints an statement indicating that the script is done, the number of total contigs and the number of viral contigs
printf "############################################################\nCOMPLETED! :)\nRESULTS ---->\nThe total number of contigs assembled by rnaSPAdes is = $number_contigs\nThe number of potyviral contigs is = $number_viral_contigs\n"

#Prints statements indicating the output files location and content of each
printf "############################################################\nThe file 'annotated_potyvirus_contigs.fasta' is located in the 'Results' folder. This file contains the annotated viral contigs with the identity of the reference sequence to which they align, the contig and reference length, the alignment length and the percentage of identity\n"

printf "############################################################\nThe file 'summary.txt' is located in the 'Results' folder. This file contains a summarized table with the information for each viral contig\n############################################################\n"

#Removes the Temporal folder
rm -r Temp
