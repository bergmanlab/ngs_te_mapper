#R --no-save < ngs_te_mapper.R sample.fasta /Users/user_name/ngs_te_mapper 1 20
#align to the TEs first
#then get the reads that match to a TE
#align the selected reads to the genome of interest
# Author: raquel
###############################################################################


Args <- commandArgs();

print(Args)
sample<-Args[[3]]
directory<-Args[[4]]
repeated<-as.integer(Args[[5]])
tsd<-as.integer(Args[[6]])

print(sample)
print(directory)

seqtekCo<-"seqtk fq2fa"
blatCommand<-"blat "

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
	fastaFile<-paste(fastaFolder, sample, sep = "")
	sample<-strsplit(sample, split = "\\.")[[1]][1]
}
if(length(files) >1)
{
	
	fastaFile<-paste(fastaFolder, "temp", RandomString(), sep = "")
	myOutput<-file(fastaFile, "w")
	fileName<-NULL
	for ( i in 1:length(files))
	{
		fileName[i]<-strsplit(files[i], split = "\\.")[[1]][1]
		cat(system(paste("sed 's/ .*/_", fileName[i], "/g' ",fastaFolder, files[i],
						sep = ""), intern = TRUE), file = myOutput, sep = "\n")
	}
	close(myOutput)	
	sample<-paste(fileName, collapse = "_")
}
aligned<-paste(firstAlignFolder, sample,".psl", sep= "" )
secondFastaFile<-paste(secondFastaFolder,sample,".fasta", sep= "" )
lastBlatFile<-paste(lastBlatFolder,sample,".psl", sep= "" )
dataFile<-paste(dataFolder,sample, ".Rdata", sep= "" )
bedFileReads<-paste(bedFolder, sample, "reads",sep = "")
bedFileInsertions<-paste(bedFolder, sample, "insertions.bed",sep = "")
outputFile<-paste(outputFolder, sample, "insertions.tsv",sep = "")

teFile<-system(paste("ls ", directory, "/reference/te", sep = ""), intern = TRUE)[1]
teFile<-paste(directory, "/reference/te/", teFile, sep = "")
organism<-system(paste("ls ", directory, "/reference/genome", sep = ""), intern = TRUE)[1]
organism<-paste(directory, "/reference/genome/", organism, sep = "")


#align to the TE dataset
system(paste(blatCommand, teFile, " ", fastaFile, " ", aligned, sep = ""))

#select the reads that map uniquely to the end and start of the TE
aPslFile<-GetPSL(aligned)
selectedReads<-SelectFirstReads(aPslFile)

#right the reads to a new fasta file
RightNewFasta(selectedReads, fastaFile, secondFastaFile )
system(paste(blatCommand, organism, " ", secondFastaFile, " ", lastBlatFile))
system(paste("rm ", fastaFile))
#select for the reads that are mapped to both the genome and the TE
#allow for reads to be mapped to different sites by the number given from repeated
otherPslFile<-GetPSL(lastBlatFile)
secondReads<-SelectSecondReads(otherPslFile, bedFileReads, withRep = repeated)
myLocations<-FinalProcessing(secondReads$toKeep,sample, tsd = tsd)
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

q(save = "no")
