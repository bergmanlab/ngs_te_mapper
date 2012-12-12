#!/usr/bin/bash 

#set up base directory for project 
projectdir="$PWD"
projectdir=$projectdir/example
#set up directories for input and output files
mkdir $projectdir/samples
mkdir $projectdir/samples/fastq
mkdir $projectdir/samples/fasta
mkdir $projectdir/reference
mkdir $projectdir/reference/te
mkdir $projectdir/reference/genome
mkdir $projectdir/analysis/
mkdir $projectdir/analysis/psl_te
mkdir $projectdir/analysis/fasta_aligned_te
mkdir $projectdir/analysis/psl_genome
mkdir $projectdir/analysis/bed_tsd
mkdir $projectdir/analysis/metadata
mkdir $projectdir/analysis/r_data_files
mkdir $projectdir/logo

#fetch and install reference genome
wget ftp://ftp.flybase.net/releases/current/dmel_r5.48/fasta/dmel-all-chromosome-r5.48.fasta.gz
gunzip dmel-all-chromosome-r5.48.fasta.gz
mv dmel-all-chromosome-r5.48.fasta $projectdir/reference/genome

#fetch and install reference TE set
wget http://www.fruitfly.org/data/p_disrupt/datasets/ASHBURNER/VERSION9.4.1.zip
unzip VERSION9.4.1.zip 
mv VERSION9.4.1/D_mel_transposon_sequence_set.fasta.v9.4.1 $projectdir/reference/te
rm -rf VERSION9.4.1*

#copy (test) input files into ngs input directory
cp example/sample.fasta $projectdir/samples/fasta/sample1.fasta
cp example/sample.fasta $projectdir/samples/fasta/sample2.fasta

#run ngs_te_mapper on all files ngs input directory
for input in $projectdir/samples/fasta/*
do 
sample=`basename $input`
 R --no-save < sourceCode/ngs_te_mapper.R $sample $projectdir 1 20
done


R --no-save < sourceCode/ngs_te_logo.R $projectdir 25
