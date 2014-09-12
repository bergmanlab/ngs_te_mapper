#!/usr/bin/bash 
#set up base directory for project 
projectdir="$PWD"
projectdir=$projectdir
#set up directories for input and output files

#fetch and install reference genome
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz
tar -zxvf chromFa.tar.gz
rm chromFa.tar.gz
cat chrYHet.fa chrM.fa chr2L.fa chrX.fa chr3L.fa chr4.fa chr2R.fa chr3R.fa chrUextra.fa chr2RHet.fa chr2LHet.fa chr3LHet.fa chr3RHet.fa chrU.fa chrXHet.fa > dm3.fasta
mv dm3.fasta $projectdir/reference/genome
rm chrYHet.fa chrM.fa chr2L.fa chrX.fa chr3L.fa chr4.fa chr2R.fa chr3R.fa chrUextra.fa chr2RHet.fa chr2LHet.fa chr3LHet.fa chr3RHet.fa chrU.fa chrXHet.fa

#fetch TE set
wget https://raw.githubusercontent.com/cbergman/transposons/master/misc/D_mel_transposon_sequence_set.fa
mv D_mel_transposon_sequence_set.fa $projectdir/reference/te

#just until the new fastq files are not here 
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA018/SRA018188/SRX021061/SRR834530_1.fastq.bz2
bzip2 -d SRR834530_1.fastq.bz2
mv SRR834530_1.fastq $projectdir/example/sample1_1.fastq

wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA018/SRA018188/SRX021061/SRR834530_2.fastq.bz2
bzip2 -d SRR834530_2.fastq.bz2
mv SRR834530_2.fastq $projectdir/example/sample1_2.fastq

#run ngs_te_mapper on all files ngs input directory
#for input in $projectdir/samples/fasta/*
#do 
#R --no-save < sourceCode/ngs_te_mapper.R "--args sample='$input' genome='$projectdir/reference/genome/dm3.fasta' teFile='$projectdir/reference/te/D_mel_transposon_sequence_set.fa' output='$projectdir/analysis'"
#done

#run ngs_te_mapper on different files has if it was only one sample (for paired end)
#the names of the files have to be separated by ";"

#R --no-save < sourceCode/ngs_te_mapper.R --args "sample1_1.fasta;sample1_2.fasta" $projectdir 1 20 20

R --no-save < sourceCode/ngs_te_mapper.R "--args sample='$projectdir/example/sample1_1.fastq;$projectdir/example/sample1_2.fastq' genome='$projectdir/dm3.fasta' teFile='$projectdir/D_mel_transposon_sequence_set.fa' repeated=1 tolerance=20 tsd=20 output='$projectdir/analysis'"

R --no-save < sourceCode/ngs_te_logo.R "--args genome='$projectdir/dm3.fasta' output='$projectdir/analysis/logo' inputFolder='$projectdir/analysis/metadata' outputFile='$projectdir/analysis/allSamples.bed' window=25"

