# R --no-save < ngs_te_mapper.R /Users/user_name/ngs_te_mapper 25
# Author: raquel
###############################################################################

Args <- commandArgs();
print(Args)
directory<-Args[[3]]
seqWindow<-Args[[4]]

logoFolder<-paste(directory, "/logo/", sep = "")
inputFolder<-paste(directory, "/analysis/metadata/", sep = "")
organism<-system(paste("ls ", directory, "/reference/genome", sep = ""), intern = TRUE)[1]
organism<-paste(directory, "/reference/genome/", organism, sep = "")

source("sourceCode/ngs_te_mapper_functions.R")

aFastaFile<-GetFasta(organism, sizeLocation = NA)		

allTogether<-paste(analysisFolder, "allSamples", sep = "")
myOutput<-file(allTogether, "w")
teInsertionData<-system(paste("grep -v repeatedStart ", inputFolder, '* | cut -d ":" -f2 |sort ', sep = ""), 
		intern = TRUE) 
cat(paste("chrom", "start", "end", "tsd", "strand", "teName", "strain", "nReads", "repeatedStart", 
				sep = "\t"), sep = "\n", file = myOutput)
cat(paste(teInsertionData[,1],teInsertionData[,2], teInsertionData[,3], teInsertionData[,4],teInsertionData[,5], 
			teInsertionData[,6],teInsertionData[,7],teInsertionData[,8], teInsertionData[,9], sep = "\t"), 
	sep = "\n", file = myOutput)
close(myOutput)

aGraphName<-paste(analysisFolder, "allSamples.ps", sep = "")

teNames<-names(table(teInsertionData[,6]))
postscript(aGraphName, paper = 'a4', horizontal = TRUE);
for ( i in 1:length(teNames))
{
	tempo<-which(teInsertionData[,6] == teNames[i])
	mySequences<- GetSequences(aFastaFile, teInsertionData[tempo,1], as.numeric(teInsertionData[tempo,2]),seqWindow,
			seqWindow, teInsertionData[tempo,5])
	Logo(mySequences$matrix, title =paste(teNames[i], length(tempo), sep = "\t"), 
			start = -seqWindow, ytitle = c("score"), xtitle = c("Position relative to insertion site"))
}
dev.off()

