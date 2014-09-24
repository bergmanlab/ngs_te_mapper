#!/usr/bin/Rscript --vanilla 
#chmod u+x ngs_te_mapper.R
#./ngs_te_mapper.R sample.fasta /Users/user_name/ngs_te_mapper 1 20
#align to the TEs first
#then get the reads that match to a TE
#align the selected reads to the genome of interest
# Author: raquel
###############################################################################

args <- commandArgs(trailingOnly = TRUE);
print(sessionInfo())
cat("\n")

if (length(args) ==0 )
{
	cat(paste("need to have arguments passed in ex:", "\n", "sourceCode/ngs_te_mapper.R sample=~/ngs_te_mapper/example/sample1.fastq\\;~/ngs_te_mapper/example/sample2.fastq genome=~/ngs_te_mapper/reference/genome/dm3.fasta teFile=~/ngs_te_mapper/reference/genome/dm3.fasta output=~/ngs_te_mapper/analysis sourceCodeFolder=~/ngs_te_mapper/sourceCode\n", sep = ""))
	q(save = "no")
}

print(args)
sample<-NA
genome<-NA
teFile<-NA
output<-NA
fastaFolder<-NA
tsd<-NA
sourceCodeFolder<-NA
for( i in 1:length(args))
{
	test<-strsplit(args[[i]], split = "=")
	assign(test[[1]][1], test[[1]][2], envir=.GlobalEnv);
#	eval(parse(text=args[[i]]))
}

if (is.na(sample) == TRUE)
{
	print("need a sample to process ex: sample='sample1.fastq;sample2.fastq'")
	q(save = "no")
}
if(is.na(genome) == TRUE)
{
	print("need the full path to the genome file ex: genome='~/ngs_te_mapper/reference/genome/dm3.fasta'")
	q(save = "no")
}
if(is.na(teFile) == TRUE)
{
	print("need the full path to the TE file ex: teFile='~/ngs_te_mapper/reference/genome/dm3.fasta'")
	q(save = "no")
}
if(is.na(output) == TRUE)
{
	print("need the full path to the output folder ex: output='~/ngs_te_mapper/analysis'")
	q(save = "no")
}
if(is.na(fastaFolder) == TRUE)
{
	if(dirname(strsplit(sample, split = ";")[[1]][1]) == "." )
	{
		print("need the path to the samples since only the sample filename was supplied ex: fastaFolder='~/ngs_te_mapper/samples/'")
		q(save = "no")	
	}
}
if(is.na(sourceCodeFolder) == TRUE)
{
	print("need the full path to the sourceCodeFolder folder ex: output='~/ngs_te_mapper/sourceCode'")
	q(save = "no")
}
if(is.na(tsd) == FALSE)
{
	tsd<-as.numeric(tsd)
}
if(is.na(tsd) == TRUE)
{
	tsd<-20
}
cat(paste("going to analyse: ",sample, "\n", "in to: ",output,"\nwith \n\t",
	"maximum size of TSD: ", tsd, "\n", sep = ""))

bwaCommand<-("bwa mem ")
bwaIndex<-"bwa index -p "
source(paste(sourceCodeFolder, "/ngs_te_mapper_functions.R", sep = ""))
#In here it will look if there is an index in the given Te and organism folders and if it does not exist it will create a new one
referenceTE<-unlist(strsplit(teFile, split = "\\.fa"))[1]
aList<-system(paste("ls -rt ", dirname(referenceTE), sep = ""), intern = TRUE)
if (length(grep(paste(referenceTE, ".bwt", sep = ""), paste(dirname(referenceTE),"/", aList, sep = ""))) == 0)
{
	print("getting the bwa index for the TE set")
	system(paste( bwaIndex, referenceTE, teFile, sep = " "))	
}

referenceGenome<-unlist(strsplit(genome, split = "\\.fa"))[1]
aList<-system(paste("ls -rt ", dirname(referenceGenome), sep = ""), intern = TRUE)
if (length(grep(paste(referenceGenome, ".bwt", sep = ""), paste(dirname(referenceGenome),"/", aList, sep = ""))) == 0 )
{
	print("getting the bwa index for the Genome")
	system(paste( bwaIndex, referenceGenome, genome  ,sep = " "))	
}


#q(save = "no")
#get the folders right
firstAlignFolder<-paste(output, "/aligned_te/", sep = "")
secondFastaFolder<-paste(output, "/selected_reads/", sep = "")
lastAlignFolder<-paste(output, "/aligned_genome/", sep = "")
#dataFolder<-paste(output, "/r_data_files/", sep = "")
bedFolder<-paste(output, "/bed_tsd/", sep = "")
#outputFolder<-paste(output, "/metadata/", sep = "")

#make the folders here. if they already exist it will through an error 
system(paste("mkdir ", output, sep = ""))
system(paste("mkdir ", firstAlignFolder, sep = ""))
system(paste("mkdir ", secondFastaFolder, sep = ""))
system(paste("mkdir ", lastAlignFolder, sep = ""))
#system(paste("mkdir ", dataFolder, sep = ""))
system(paste("mkdir ", bedFolder, sep = ""))
#system(paste("mkdir ", outputFolder, sep = ""))


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
			#only the name of the file is given look for the file in the sample folder
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
#dataFile<-paste(dataFolder,sample, ".Rdata", sep= "" )
bedFileReads<-paste(bedFolder, sample, "reads",sep = "")
bedFileInsertions<-paste(bedFolder, sample, "insertions.bed",sep = "")
#outputFile<-paste(outputFolder, sample, "insertions.tsv",sep = "")



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
selectedReads<-SelectFirstReadsSam(aSamFile, secondFastaFile = secondFastaFile )
system(paste(bwaCommand, referenceGenome, " ", secondFastaFile, " >", lastFile,".sam" ,sep = ""))

#select for the reads that are mapped to both the genome and the TE
#allow for reads to be mapped to different sites by the number given from repeated
otherSamFile<-GetSamFile(paste(lastFile, ".sam", sep = ""))
#secondReads<-SelectSecondReadsSam(otherSamFile, paste(bedFileReads, "new.bed", sep = ""), tolerated =tolerance )
secondReads<-SelectSecondReadsSam(otherSamFile)
myLocations<-FinalProcessingSam(secondReads$toKeep,sample, tsd = tsd)
if (length(myLocations) == 0)
{
	cat(paste("there were no new insertions found with the following criteria: ", "\n","\t",
			"maximum size of TSD: ", tsd, "\n", sep = ""))
	q(save = "no")
}

#now for the old TEs
#otherSamFile<-GetSamFile(paste(lastFile, ".sam", sep = ""))
secondReadsOld<-SelectSecondReadsSamOld(otherSamFile,aSamFile)
minDist<-5*tsd
myLocationsOld<-FinalProcessingSamOldTes(secondReadsOld$toKeep,sample, aSamFile, minDist)
if (length(myLocationsOld) == 0)
{
	cat(paste("there were no old insertions found with the following criteria: ", "\n","\t",
					"minimum size of TE: ", minDist, "\n", "\t",
					"maximum size of TE: 1.5 * Te length",  "\n", sep = ""))
	q(save = "no")
}

#######
myLocations<-paste(myLocations, "new", sep = ";")
myLocationsOld<-paste(myLocationsOld, "old", sep = ";")
myLocations2<-matrix(data = c(unlist(strsplit(myLocations, split = ";")), unlist(strsplit(myLocationsOld, split = ";"))), nrow= length(myLocations)+length(myLocationsOld), byrow = TRUE)
myLocations2<-as.data.frame(myLocations2)
myLocations2$V2<-as.numeric(as.character(myLocations2$V2))
myLocations2$V3<-as.numeric(as.character(myLocations2$V3))
myLocations2<-myLocations2[order(myLocations2$V1,myLocations2$V2,myLocations2$V3),]
myOutput<-file(bedFileInsertions, "w")
temp<-paste(myLocations2[,4],myLocations2[,5], myLocations2[,6],myLocations2[,7],myLocations2[,8],myLocations2[,9], sep = ";")
cat(paste(myLocations2[,1],myLocations2[,2], myLocations2[,3], temp, sep = "\t"), sep = "\n", file = myOutput)
close(myOutput)

#myOutput<-file(outputFile, "w")
#cat(paste("chrom", "start", "end", "tsd_length", "strand", "teName", "strain", "nReads","insertion", sep = "\t"), sep = "\n", file = myOutput)
#cat(paste(myLocations2[,1],myLocations2[,2], myLocations2[,3], myLocations2[,4],myLocations2[,5], myLocations2[,6], myLocations2[,7],myLocations2[,8],myLocations2[,9], sep = "\t"), sep = "\n", file = myOutput)
#close(myOutput)

cat(paste("finished the job and found ", length(myLocations), " new insertions run: \n\nsourceCode/ngs_te_logo.R  genome=",genome, " output=",output, "/logo inputFolder=", output, "/metadata outputFile=", output, "/allSamples.bed window=25", "\nto get the logos centred at the TSD with +/- 25 bp to both sides\n\n", sep = ""))
cat("also found ", length(myLocationsOld), " previously known insertions\n")

#save(list = ls(), file =dataFile)
q(save = "no")
