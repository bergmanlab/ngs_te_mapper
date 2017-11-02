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

Rscript --vanilla $projectdir/sourceCode/ngs_te_mapper.R \
sample=$projectdir/example/sample1_1.fastq\;$projectdir/example/sample1_2.fastq \
genome=$projectdir/dm3.fasta \
teFile=$projectdir/D_mel_transposon_sequence_set.fa \
tsd=20 \
output=$projectdir/analysis \
sourceCodeFolder=$projectdir/sourceCode 

Rscript --vanilla $projectdir/sourceCode/ngs_te_logo.R \
genome=$projectdir/dm3.fasta \
output=$projectdir/analysis/logo \
inputFolder=$projectdir/analysis/bed_tsd \
outputFile=$projectdir/analysis/allSamples.bed \
window=25 \
sourceCodeFolder=$projectdir/sourceCode 

nDiffentLines="`diff -y --suppress-common-lines $projectdir/expectedBedFile.bed $projectdir/analysis/bed_tsd/sample1_1_sample1_2insertions.bed | grep -c '[[:space:]]|[[:space:]]'`"
nNotPresentInNewFile="`diff -y --suppress-common-lines $projectdir/expectedBedFile.bed $projectdir/analysis/bed_tsd/sample1_1_sample1_2insertions.bed | grep -c '[[:space:]]>[[:space:]]'`"
nPresentOnlyNewFile="`diff -y --suppress-common-lines $projectdir/expectedBedFile.bed $projectdir/analysis/bed_tsd/sample1_1_sample1_2insertions.bed | grep -c '[[:space:]]<'`"

if [ $nDiffentLines -eq 0 ]; then
	echo "There was no insertion sites that were different between the test .bed file and the new .bed file from ~analysis/bed_tsd "
else 
	echo "There was $nDiffentLines insertion site/s that was/were different between the test .bed file and the new .bed file from ~analysis/bed_tsd "
fi

if [ $nNotPresentInNewFile -eq 0 ]; then
	echo "All insertion sites in the test .bed file were also present in the new .bed file from ~analysis/bed_tsd"
else
	echo "There was $nNotPresentInNewFile insertion site/s that was/were present in the test .bed file and not in the new .bed file from ~analysis/bed_tsd "
fi

if [ $nPresentOnlyNewFile -eq 0 ]; then
	echo "All insertion sites in the new .bed file from ~analysis/bed_tsd were also present in the test .bed"
else
	echo "There was $nPresentOnlyNewFile insertion site/s that was/were present in the new .bed file from ~analysis/bed_tsd and not in the test .bed file"
fi

