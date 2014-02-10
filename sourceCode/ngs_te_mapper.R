#!/usr/bin/Rscript --vanilla 
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

bwaCommand<-("bwa mem ")
bwaIndex<-"bwa index -p "

#q(save = "no")
#get the folders right
analysisFolder<-paste(directory, "/analysis/", sep = "")
fastaFolder<-paste(directory, "/samples/fasta/", sep = "")		
firstAlignFolder<-paste(analysisFolder, "align_te/", sep = "")
secondFastaFolder<-paste(analysisFolder, "fasta_aligned_te/", sep = "")
lastAlignFolder<-paste(analysisFolder, "align_genome/", sep = "")
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
	fileName<-NULL
	for ( i in 1:length(files))
	{
		if(basename(files[i]) == files[i])
		{
			#only the name of the file is given look for the file in the fasta folder
			fileName[i]<-strsplit(files[i], split = "\\.")[[1]][1]
		}
		else
		{
			fileName[i]<-strsplit(basename(files[i]), split = "\\.")[[1]][1]
		}
	}
	sample<-paste(fileName, collapse = "_")
}
aligned<-paste(firstAlignFolder, sample, sep= "" )
secondFastaFile<-paste(secondFastaFolder,sample,".fasta", sep= "" )
lastFile<-paste(lastAlignFolder,sample, sep= "" )
dataFile<-paste(dataFolder,sample, ".Rdata", sep= "" )
bedFileReads<-paste(bedFolder, sample, "reads",sep = "")
bedFileInsertions<-paste(bedFolder, sample, "insertions.bed",sep = "")
outputFile<-paste(outputFolder, sample, "insertions.tsv",sep = "")

teFile<-system(paste("ls -rt ", directory, "/reference/te", sep = ""), intern = TRUE)[1]
teFile<-paste(directory, "/reference/te/", teFile, sep = "")
organism<-system(paste("ls ", directory, "/reference/genome", sep = ""), intern = TRUE)[grep(
	"fasta", system(paste("ls ", directory, "/reference/genome", sep = ""), intern = TRUE))]
organism<-paste(directory, "/reference/genome/", organism, sep = "")

referenceTE<-unlist(strsplit(teFile, split = "\\."))[1]
if( teFile != paste(referenceTE, ".fasta", sep =""))
{
	system(paste("mv ", teFile, " ",referenceTE, ".fasta", sep = ""))
	teFile<-paste(referenceTE, ".fasta", sep = "")
}
aList<-system(paste("ls -rt ", directory, "/reference/te", sep = ""), intern = TRUE)
if (length(grep(paste(referenceTE, ".bwt", sep = ""), paste(directory, "/reference/te/", aList, sep = ""))) == 0 | 
		length(grep(paste(referenceTE, ".bwt", sep = ""), paste(directory, "/reference/te/", aList, sep = ""))) == 0 )
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
if (length(grep(paste(referenceGenome, ".bwt", sep = ""), paste(directory, "/reference/genome/",aList, sep = ""))) == 0 | 
		length(grep(paste(referenceGenome, ".bwt", sep = ""), paste(directory, "/reference/genome/",aList, sep = ""))) == 0 )
{
	print("getting the bwa index for the Genome")
	system(paste( bwaIndex, referenceGenome, organism  ,sep = " "))	
}


#start the analysis
#align to the TE dataset

if(length(files) == 1)
{
	system(paste(bwaCommand, referenceTE, " ", fastaFile, " >", aligned,".sam" ,sep = ""))
	aSamFile<-GetSamFile(paste(aligned,".sam", sep = ""))
}

if(length(files) >1)
{
	for ( i in 1:length(files))
	{
		if(basename(files[i]) == files[i])
		{
			#only the name of the file is given look for the file in the fasta folder
			fileName<-strsplit(files[i], split = "\\.")[[1]][1]
			system(paste(bwaCommand, referenceTE, " ",fastaFolder, files[i]," >", firstAlignFolder , fileName, ".sam" ,sep = ""))
		}
		else
		{
			fileName<-strsplit(basename(files[i]), split = "\\.")[[1]][1]
			system(paste(bwaCommand, referenceTE, " ", files[i]," >", firstAlignFolder , fileName,".sam" ,sep = ""))
		}
	}
	i<-1
	if(basename(files[i]) == files[i])
	{
		fileName<-strsplit(files[i], split = "\\.")[[1]][1]
		aSamFile<-GetSamFile(paste(firstAlignFolder , fileName,".sam", sep = ""))	
		aSamFile$QNAME<-paste(aSamFile$QNAME, fileName, sep = "_")
	}
	else
	{
		fileName<-strsplit(basename(files[i]), split = "\\.")[[1]][1]
		aSamFile<-GetSamFile(paste(firstAlignFolder , fileName,".sam", sep = ""))
		aSamFile$QNAME<-paste(aSamFile$QNAME, fileName, sep = "_")
	}
	for ( i in 2:length(files))
	{
		if(basename(files[i]) == files[i])
		{
			fileName<-strsplit(files[i], split = "\\.")[[1]][1]
			tempSamFile<-GetSamFile(paste(firstAlignFolder , fileName,".sam", sep = ""))	
			tempSamFile$QNAME<-paste(tempSamFile$QNAME, fileName, sep = "_")
		}
		else
		{
			
			fileName<-strsplit(basename(files[i]), split = "\\.")[[1]][1]
			tempSamFile<-GetSamFile(paste(firstAlignFolder , fileName,".sam", sep = ""))	
			tempSamFile$QNAME<-paste(tempSamFile$QNAME, fileName, sep = "_")
		}
		aSamFile<-mapply(c, aSamFile, tempSamFile, SIMPLIFY=FALSE)
	}
	aSamFile$names<-tempSamFile$names
	aSamFile$length<-tempSamFile$length
}



#select the reads that map uniquely to the end and start of the TE
selectedReads<-SelectFirstReadsSam(aSamFile, tolerated =tolerance ,secondFastaFile = secondFastaFile )
system(paste(bwaCommand, referenceGenome, " ", secondFastaFile, " >", lastFile,".sam" ,sep = ""))

#select for the reads that are mapped to both the genome and the TE
#allow for reads to be mapped to different sites by the number given from repeated
otherSamFile<-GetSamFile(paste(lastFile, ".sam", sep = ""))
secondReads<-SelectSecondReadsSam(otherSamFile, paste(bedFileReads, "test", sep = ""), tolerated =tolerance )
myLocations<-FinalProcessingSam(secondReads$toKeep,sample, tsd = tsd)

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
cat(paste("chrom", "start", "end", "tsd", "strand", "teName", "strain", "nReads", sep = "\t"), sep = "\n", file = myOutput)
cat(paste(myLocations2[,1],myLocations2[,2], myLocations2[,3], myLocations2[,4],myLocations2[,5], myLocations2[,6], 
				myLocations2[,7],myLocations2[,8], sep = "\t"), sep = "\n", file = myOutput)
close(myOutput)

cat("finished the job and found ", length(myLocations), " new insertions run: \nR --no-save < sourceCode/ngs_te_logo.R ",
		directory, " 25\nto get the logos centred at the TSD with +/- 25 bp to both sides\n")

q(save = "no")
