#!/usr/bin/env Rscript
# R --no-save < ngs_te_logo.R /Users/user_name/ngs_te_mapper 25
# Author: raquel
###############################################################################


args <- commandArgs(trailingOnly = TRUE);
print(sessionInfo())

source("sourceCode/ngs_te_mapper_functions.R")
cat("\n")
if (length(args) ==0 )
{
	cat(paste("need to have arguments passed in ex:", "\n", "sourceCode/ngs_te_logo.R genome=~/ngs_te_mapper/reference/genome/dm3.fasta output=~/ngs_te_mapper/analysis/logo inputFolder=~/analysis/metadata outputFile=~/analysis/allSamples.bed window=25 sourceCodeFolder=~/ngs_te_mapper/sourceCode\n", sep = ""))
	q(save = "no")
}
print(args)

genome<-NA
output<-NA
outputFile<-NA
inputFolder<-NA
window<-NA
for( i in 1:length(args))
{
	test<-strsplit(args[[i]], split = "=")
	assign(test[[1]][1], test[[1]][2], envir=.GlobalEnv);
#	eval(parse(text=args[[i]]))
}
rm(test)
if(is.na(genome) == TRUE)
{
	print("need the full path to the genome file ex: genome='~/ngs_te_mapper/reference/genome/dm3.fasta'")
	q(save = "no")
}
if(is.na(output) == TRUE)
{
	print("need the full path to the output folder ex: output='~/ngs_te_mapper/analysis/logo'")
	q(save = "no")
}
if(is.na(inputFolder) == TRUE)
{
	print("need the path to input folder that contains the results from the analysis from ngs_te_mapper.R ex: inputFolder='$projectdir/analysis/metadata'")
	q(save = "no")	
}
if(is.na(outputFile) == TRUE)
{	
	print("need the full path to the output file that will have all the samples together ex: output='~/analysis/allSamples'")
	q(save = "no")
}
if(is.na(sourceCodeFolder) == TRUE)
{
	print("need the full path to the sourceCodeFolder folder ex: output='~/ngs_te_mapper/sourceCode'")
	q(save = "no")
}
if(is.na(window) == FALSE)
{
	window<-as.numeric(window)
}
if(is.na(window) == TRUE)
{
	window<-25
}

source(paste(sourceCodeFolder, "/ngs_te_mapper_functions.R", sep = ""))
aFastaFile<-GetFasta(genome, sizeLocation = NA)		
myOutput<-file(outputFile, "w")

teInsertionData<-strsplit(system(paste("grep --color=never . ", inputFolder, '/*insertions.bed | cut -d ":" -f2 | sort | sed "s/;/\t/g"', sep = ""), intern = TRUE), split = "\t")
teInsertionData<-matrix(data = unlist(teInsertionData), ncol = length(teInsertionData[[1]]), nrow = length(teInsertionData),byrow = TRUE)
cat(paste("chrom", "start", "end", "tsd", "strand", "teName", "strain", "nReads", sep = "\t"), sep = "\n", file = myOutput)
teInsertionData<-teInsertionData[order(teInsertionData[,1],teInsertionData[,2], teInsertionData[,3]),]
cat(paste(teInsertionData[,1],teInsertionData[,2], teInsertionData[,3], teInsertionData[,4],teInsertionData[,5], teInsertionData[,6],teInsertionData[,7],teInsertionData[,8], sep = "\t"), sep = "\n", file = myOutput)
close(myOutput)

#if it exists it will just through an error
system(paste("mkdir ", output, sep = ""))
aGraphName<-paste(output, "/allSamples.pdf", sep = "")

teNames<-names(table(teInsertionData[which(teInsertionData[,9] == "new"),6]))
pdf(aGraphName, height=8.3,  width= 11.7);
for ( i in 1:length(teNames))
{
    tempo<-which(teInsertionData[,6] == teNames[i] & teInsertionData[,9] == "new")
    if(names(table( teInsertionData[tempo,5]))[1] == "NA")
    {
        next;
    }
	mySequences<- GetSequences(aFastaFile, teInsertionData[tempo,1], as.numeric(teInsertionData[tempo,2]),window, window, teInsertionData[tempo,5])
	Logo(mySequences$matrix, title =paste(teNames[i], length(mySequences$id), sep = " "), start = -window, ytitle = c("score"), xtitle = c("Position relative to insertion site"))
}
dev.off()

