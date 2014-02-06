#!/usr/bin/Rscript --vanilla --slave
#chmod u+x ngs_te_mapper.R
#./ngs_te_mapper.R sample.fasta /Users/user_name/ngs_te_mapper 1 20
#align to the TEs first
#then get the reads that match to a TE
#align the selected reads to the genome of interest
# Author: raquel
###############################################################################


Args <- commandArgs();
toAddArgs<-2
print(sessionInfo())
cat("\n")
print(Args)
sample<-Args[1+toAddArgs]
directory<-Args[2+toAddArgs]
if(length(Args) == (2+toAddArgs))
{
	Args[3+toAddArgs]<-1
	Args[4+toAddArgs]<-20
	Args[5+toAddArgs]<-20
}
if(length(Args) == (3+toAddArgs))
{
	Args[4+toAddArgs]<-20
	Args[5+toAddArgs]<-20
}
if(length(Args) == (4+toAddArgs))
{
	Args[5+toAddArgs]<-20
}
repeated<-as.integer(Args[3+toAddArgs])
tolerance<-as.integer(Args[4+toAddArgs])
tsd<-as.integer(Args[5+toAddArgs])

cat(paste("going to analysze: ",sample, "\n", "in to: ",directory,"\nwith \n\t",
	"number of matches per read in the genome: ",repeated,"\n", "\t", 
	"proportion of missmatches: ",tolerance,"\n", "\t",
	"maximum size of TSD: ", tsd, "\n", sep = ""))

seqtekCo<-"seqtk fq2fa"
blatCommand<-"blat "
bwaCommand<-"bwa mem "
bwaIndex<-"bwa index -p "



#q(save = "no")
#get the folders right
analysisFolder<-paste(directory, "/analysis/", sep = "")
fastaFolder<-paste(directory, "/samples/fasta/", sep = "")		
firstAlignFolder<-paste(analysisFolder, "psl_te/", sep = "")
secondFastaFolder<-paste(analysisFolder, "fasta_aligned_te/", sep = "")
lastBlatFolder<-paste(analysisFolder, "psl_genome/", sep = "")
dataFolder<-paste(analysisFolder, "r_data_files/", sep = "")
bedFolder<-paste(analysisFolder, "bed_tsd/", sep = "")
outputFolder<-paste(analysisFolder, "metadata/", sep = "")

source("sourceCode/ngs_te_mapper_functions.R")

#now for the files
files<-strsplit(sample, split =";")[[1]]
if(length(files) == 1)
{
	if(basename(sample) == sample)
	{
		#only the name of the file is given look for the file in the fasta folder
		fastaFile<-paste(fastaFolder, sample, sep = "")
		sample<-strsplit(sample, split = "\\.")[[1]][1]
	}
	else
	{
		fastaFile<-sample
		sample<-unlist(strsplit(basename(sample), split = "\\."))[1]
	}
}

if(length(files) >1)
{
	
	fastaFile<-paste(fastaFolder, "temp", RandomString(), sep = "")
	myOutput<-file(fastaFile, "w")
	fileName<-NULL
	for ( i in 1:length(files))
	{
		if(basename(files[i]) == files[i])
		{
			#only the name of the file is given look for the file in the fasta folder
			fileName[i]<-strsplit(files[i], split = "\\.")[[1]][1]
			cat(system(paste("sed 's/ .*/_", fileName[i], "/g' ",fastaFolder, files[i],
				sep = ""), intern = TRUE), file = myOutput, sep = "\n")
		}
		else
		{
			fileName[i]<-strsplit(basename(files[i]), split = "\\.")[[1]][1]
			cat(system(paste("sed 's/ .*/_", fileName[i], "/g' ", files[i],	sep = ""), 
				intern = TRUE), file = myOutput, sep = "\n")
		}
	}
	close(myOutput)	
	sample<-paste(fileName, collapse = "_")
}
aligned<-paste(firstAlignFolder, sample, sep= "" )
secondFastaFile<-paste(secondFastaFolder,sample,".fasta", sep= "" )
lastBlatFile<-paste(lastBlatFolder,sample, sep= "" )
dataFile<-paste(dataFolder,sample, ".Rdata", sep= "" )
bedFileReads<-paste(bedFolder, sample, "reads",sep = "")
bedFileInsertions<-paste(bedFolder, sample, "insertions.bed",sep = "")
outputFile<-paste(outputFolder, sample, "insertions.tsv",sep = "")

teFile<-system(paste("ls -rt ", directory, "/reference/te", sep = ""), intern = TRUE)[1]
teFile<-paste(directory, "/reference/te/", teFile, sep = "")
organism<-system(paste("ls ", directory, "/reference/genome", sep = ""), intern = TRUE)[1]
organism<-paste(directory, "/reference/genome/", organism, sep = "")

referenceTE<-unlist(strsplit(teFile, split = "\\."))[1]
if( teFile != paste(referenceTE, ".fasta", sep =""))
{
	system(paste("mv ", teFile, " ",referenceTE, ".fasta", sep = ""))
	teFile<-paste(referenceTE, ".fasta", sep = "")
}
aList<-system(paste("ls -rt ", directory, "/reference/te", sep = ""), intern = TRUE)
if (length(grep(paste(directory, "/reference/te/",referenceTE, ".bwt", sep = ""), aList)) == 0 | 
		length(grep(paste(directory, "/reference/te/",referenceTE, ".bwt", sep = ""), aList)) == 0 )
{
	print("getting the bwa index for the TE set")
	system(paste( bwaIndex, referenceTE, teFile, sep = " "))	
}

referenceGenome<-unlist(strsplit(organism, split = "\\."))[1]
if(organism!= paste(referenceGenome, ".fasta", sep =""))
{
	system(paste("mv ", organism, " ",referenceGenome, ".fasta", sep = ""))
	organism<-paste(referenceGenome, ".fasta", sep = "")
}
aList<-system(paste("ls -rt ", directory, "/reference/genome/", sep = ""), intern = TRUE)
if (length(grep(paste(directory, "/reference/genome/",referenceGenome, ".bwt", sep = ""), aList)) == 0 | 
		length(grep(paste(directory, "/reference/genome/",referenceGenome, ".bwt", sep = ""), aList)) == 0 )
{
	print("getting the bwa index for the Genome")
	system(paste( bwaIndex, referenceGenome, organism, sep = " "))	
}


#align to the TE dataset
system(paste(blatCommand, teFile, " ", fastaFile, " ", aligned,".psl", sep = ""))

system(paste(bwaCommand, referenceTE, " ", fastaFile, " >", aligned,".sam" ,sep = ""))

#select the reads that map uniquely to the end and start of the TE
#aPslFile<-GetPSL(paste(aligned, ".psl", sep = ""))
#selectedReads<-SelectFirstReads(aPslFile, tolerated =tolerance )

aSamFile<-GetSamFile(paste(aligned,".sam", sep = ""))
selectedReads<-SelectFirstReadsSam(aSamFile, tolerated =tolerance ,secondFastaFile = secondFastaFile )

#right the reads to a new fasta file
RightNewFasta(selectedReads, fastaFile, secondFastaFile )
system(paste(blatCommand, organism, " ", secondFastaFile, " ", lastBlatFile, ".psl"))

system(paste(bwaCommand, referenceGenome, " ", secondFastaFile, " >", lastBlatFile,".sam" ,sep = ""))


if(length(files) >1)
{
	system(paste("rm ", fastaFile))
}
#select for the reads that are mapped to both the genome and the TE
#allow for reads to be mapped to different sites by the number given from repeated
otherPslFile<-GetPSL(paste(lastBlatFile, ".psl", sep = ""))
secondReads<-SelectSecondReads(otherPslFile, bedFileReads, withRep = repeated, tolerated =tolerance )
myLocations<-FinalProcessing(secondReads$toKeep,sample, tsd = tsd)


otherSamFile<-GetSamFile(paste(lastBlatFile, ".sam", sep = ""))
secondReads<-SelectSecondReadsSam(otherSamFile, bedFileReads, withRep = repeated, tolerated =tolerance )
myLocations<-FinalProcessing(secondReads$toKeep,sample, tsd = tsd)




if (length(myLocations) == 0)
{
	cat(paste("there were no new insertions found with the following criteria: ", "\n","\t",
			"number of matches per read in the genome: ", repeated,"\n", "\t", 
			"proportion of missmatches: ",tolerance,"\n", "\t",
			"maximum size of TSD: ", tsd, "\n", sep = ""))
	q(save = "no")
}
myLocations2<-matrix(data = unlist(strsplit(myLocations, split = ";")), nrow= length(myLocations), byrow = TRUE)
save(list = ls(), file =dataFile)

myOutput<-file(bedFileInsertions, "w")
cat(paste(myLocations2[,1],myLocations2[,2], myLocations2[,3], myLocations, sep = "\t"), sep = "\n", file = myOutput)
close(myOutput)

myOutput<-file(outputFile, "w")
cat(paste("chrom", "start", "end", "tsd", "strand", "teName", "strain", "nReads", "repeatedStart", 
				sep = "\t"), sep = "\n", file = myOutput)
cat(paste(myLocations2[,1],myLocations2[,2], myLocations2[,3], myLocations2[,4],myLocations2[,5], myLocations2[,6], 
				myLocations2[,7],myLocations2[,8], myLocations2[,9], sep = "\t"), sep = "\n", file = myOutput)
close(myOutput)

cat("finished teh job and found ", length(myLocations), " new insertions run: \nR --no-save < sourceCode/ngs_te_logo.R ",
		directory, " 25\nto get the logos centred at the TSD with +/- 25 bp to both sides\n")

q(save = "no")
