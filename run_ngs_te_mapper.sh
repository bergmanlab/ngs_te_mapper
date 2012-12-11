!/usr/bin/bash 

projectdir="$PWD"

#mkdir $projectdir/samples
#mkdir $projectdir/samples/fastq
#mkdir $projectdir/samples/fasta
#mkdir $projectdir/reference
#mkdir $projectdir/reference/te
#mkdir $projectdir/reference/genome
#mkdir $projectdir/analysis/
#mkdir $projectdir/analysis/psl_te
#mkdir $projectdir/analysis/fasta_aligned_te
#mkdir $projectdir/analysis/psl_genome
#mkdir $projectdir/analysis/bed_tsd
#mkdir $projectdir/analysis/metadata
#mkdir $projectdir/analysis/r_data_files

#wget ftp://ftp.flybase.net/releases/current/dmel_r5.48/fasta/dmel-all-chromosome-r5.48.fasta.gz
#gunzip dmel-all-chromosome-r5.48.fasta.gz
#mv dmel-all-chromosome-r5.48.fasta $projectdir/reference/genome

#wget http://www.fruitfly.org/data/p_disrupt/datasets/ASHBURNER/VERSION9.4.1.zip
#unzip VERSION9.4.1.zip 
#mv VERSION9.4.1/D_mel_transposon_sequence_set.fasta.v9.4.1 $projectdir/reference/te
#rm -rf VERSION9.4.1*

cp sample.fasta $projectdir/samples/fasta/sample1.fasta
cp sample.fasta $projectdir/samples/fasta/sample2.fasta

for input in $projectdir/samples/fasta/*
do 
sample=`basename $input`
 R --no-save < ngs_te_mapper.R $sample $projectdir 1 20
done
