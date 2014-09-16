#!/usr/bin/bash 

# set up base directory for project 
projectdir="$PWD"
projectdir=$projectdir
mkdir $projectdir/example
# fetch D. melanogaster reference genome
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz
tar -zxvf chromFa.tar.gz
rm chromFa.tar.gz
cat chrYHet.fa chrM.fa chr2L.fa chrX.fa chr3L.fa chr4.fa chr2R.fa chr3R.fa chrUextra.fa chr2RHet.fa chr2LHet.fa chr3LHet.fa chr3RHet.fa chrU.fa chrXHet.fa > dm3.fasta

rm chrYHet.fa chrM.fa chr2L.fa chrX.fa chr3L.fa chr4.fa chr2R.fa chr3R.fa chrUextra.fa chr2RHet.fa chr2LHet.fa chr3LHet.fa chr3RHet.fa chrU.fa chrXHet.fa

# fetch D. melanogaster TE set
wget --no-check-certificate https://raw.githubusercontent.com/cbergman/transposons/master/misc/D_mel_transposon_sequence_set.fa

# fetch D. melanogaster WGS reads
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR834/SRR834530/SRR834530_1.fastq.gz
gunzip SRR834530_1.fastq.gz
mv SRR834530_1.fastq $projectdir/example/sample1_1.fastq

wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR834/SRR834530/SRR834530_2.fastq.gz
gunzip SRR834530_2.fastq.gz
mv SRR834530_2.fastq $projectdir/example/sample1_2.fastq

# run ngs_te_mapper and ngs_te_logo on multiple fastq files from a single sample.
# in this case, both files from a paired end run are being mapped independently, 
# with junction reads being pooled before TSD and logo identification

chmod u+x sourceCode/ngs_te_mapper.R 
sourceCode/ngs_te_mapper.R sample='$projectdir/example/sample1_1.fastq;$projectdir/example/sample1_2.fastq' genome=$projectdir/dm3.fasta teFile=$projectdir/D_mel_transposon_sequence_set.fa repeated=1 tolerance=20 tsd=20 output=$projectdir/analysis

chmod u+x sourceCode/ngs_te_logo.R 
sourceCode/ngs_te_logo.R  genome=$projectdir/dm3.fasta output=$projectdir/analysis/logo inputFolder=$projectdir/analysis/metadata outputFile=$projectdir/analysis/allSamples.bed window=25
