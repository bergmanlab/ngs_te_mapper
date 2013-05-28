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
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz
tar -zxvf chromFa.tar.gz
rm chromFa.tar.gz
cat chrYHet.fa chrM.fa chr2L.fa chrX.fa chr3L.fa chr4.fa chr2R.fa chr3R.fa chrUextra.fa chr2RHet.fa chr2LHet.fa chr3LHet.fa chr3RHet.fa chrU.fa chrXHet.fa > dm3.fasta
mv dm3.fasta $projectdir/reference/genome
rm chrYHet.fa chrM.fa chr2L.fa chrX.fa chr3L.fa chr4.fa chr2R.fa chr3R.fa chrUextra.fa chr2RHet.fa chr2LHet.fa chr3LHet.fa chr3RHet.fa chrU.fa chrXHet.fa

#fetch and install reference TE set
wget http://www.fruitfly.org/data/p_disrupt/datasets/ASHBURNER/VERSION9.4.1.zip
unzip VERSION9.4.1.zip 
mv VERSION9.4.1/D_mel_transposon_sequence_set.fasta.v9.4.1 $projectdir/reference/te
rm -rf VERSION9.4.1*

#copy (test) input files into ngs input directory
cp example/sample1.fasta $projectdir/samples/fasta/sample1.fasta
cp example/sample2.fasta $projectdir/samples/fasta/sample2.fasta

#run ngs_te_mapper on all files ngs input directory
for input in $projectdir/samples/fasta/*
do 
sample=`basename $input`
 R --no-save < sourceCode/ngs_te_mapper.R $sample $projectdir 1 20
done

#run ngs_te_mapper on different files has if it was only one sample (for paired end)
#the names of the files have to be separated by ";"
R --no-save < sourceCode/ngs_te_mapper.R "sample1.fasta;sample2.fasta" $projectdir 1 20

R --no-save < sourceCode/ngs_te_logo.R $projectdir 25
