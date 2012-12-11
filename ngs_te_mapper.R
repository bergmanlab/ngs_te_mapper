#align to the TEs first
#then get the reads that match to a TE
#then get those reads to a lig with the genome of interest
#do something about the choice of having repeated (redundant alignements) or not
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

analysisFolder<-paste(directory, "/analysis/", sep = "")
fastaFolder<-paste(directory, "/samples/fasta/", sep = "")		
firstAlignFolder<-paste(analysisFolder, "psl_te/", sep = "")
secondFastaFolder<-paste(analysisFolder, "fasta_aligned_te/", sep = "")
lastBlatFolder<-paste(analysisFolder, "psl_genome/", sep = "")
dataFolder<-paste(analysisFolder, "r_data_files/", sep = "")
bedFolder<-paste(analysisFolder, "bed_tsd/", sep = "")
outputFolder<-paste(analysisFolder, "metadata/", sep = "")

fastaFile<-paste(fastaFolder, sample, sep = "")
aligned<-paste(firstAlignFolder, sample,".psl", sep= "" )
secondFastaFile<-paste(secondFastaFolder,sample,".fasta", sep= "" )
lastBlatFile<-paste(lastBlatFolder,sample,".psl", sep= "" )
dataFile<-paste(dataFolder,sample, sep= "" )
bedFileReads<-paste(bedFolder, sample, "reads",sep = "")
bedFileInsertions<-paste(bedFolder, sample, "insertions",sep = "")
outputFile<-paste(outputFolder, sample, "insertions",sep = "")

teFile<-system(paste("ls ", directory, "/reference/te", sep = ""), intern = TRUE)[1]
teFile<-paste(directory, "/reference/te/", teFile, sep = "")
organism<-system(paste("ls ", directory, "/reference/genome", sep = ""), intern = TRUE)[1]
organism<-paste(directory, "/reference/genome/", organism, sep = "")

source("functions.R")
save(list = ls(), file ="datafile")

#align to the TE dataset
system(paste(blatCommand, teFile, " ", fastaFile, " ", aligned, sep = ""))

#select the reads that map uniquely to the end and start of the TE
aPslFile<-GetPSL(aligned)
selectedReads<-SelectFirstReads(aPslFile)
save(list = ls(), file =dataFile)
#right the reads to a new fasta file
RightNewFasta(selectedReads, fastaFile, secondFastaFile )
system(paste(blatCommand, organism, " ", secondFastaFile, " ", lastBlatFile))

#select for the reads that are mapped to both the genome and the TE
#allow for reads to be mapped to different sites by the number given from repeated
otherPslFile<-GetPSL(lastBlatFile)
secondReads<-SelectSecondReads(otherPslFile, bedFileReads, withRep = repeated)
myLocations<-FinalProcessing(secondReads$toKeep,sample, tsd = tsd)

myLocations2<-matrix(data = unlist(strsplit(myLocations, split = ";")), nrow= length(myLocations), byrow = TRUE)
myOutput<-file(bedFileInsertions, "w")
cat(paste(myLocations2[,1],myLocations2[,2], myLocations2[,3], myLocations, sep = "\t"), sep = "\n", file = myOutput)
close(myOutput)



myOutput<-file(outputFile, "w")
cat(paste("chrom", "start", "end", "tsd", "strand", "teName", "strain", "nReads", "repeatedStart", 
				sep = "\t"), sep = "\n", file = myOutput)
cat(paste(myLocations2[,1],myLocations2[,2], myLocations2[,3], myLocations2[,4],myLocations2[,5], myLocations2[,6], 
				myLocations2[,7],myLocations2[,8], myLocations2[,9], sep = "\t"), sep = "\n", file = myOutput)
close(myOutput)

save(list = ls(), file =dataFile)
q(save = "no")



allTogether<-paste(analysisFolder, "allSamples", sep = "")
system(paste("grep SRS ", outputFolder, '* | cut -d ":" -f2 |sort > ' ,allTogether, sep = "")) 

aGraphName<-paste(analysisFolder, "allSamples.pdf", sep = "")
postscript(aGraphName, paper = 'a4', horizontal = TRUE);

seqWindow<-25

aFastaFile<-GetFasta(organism, sizeLocation = NA)
myData<-read.table(allTogether, header = TRUE)
teNames<-names(table(myLocations2[,6]))
for ( i in 1:length(teNames))
{
	selected<-which(myData[,6] == teNames[i])
	mySequences<- GetSequences(aFastaFile, myData[selected,1], as.numeric(myData[selected,2]),seqWindow, seqWindow, 
			myData[selected,5])
	Logo(mySequences$matrix, title =teNames[i], start = -25,ytitle = c("score"), 
			xtitle = c("Position relative to insertion site"))
}

dev.off()	



