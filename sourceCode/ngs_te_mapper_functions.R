#!/usr/bin/Rscript --vanilla 
# chmod u+x functions.R
# 1:00:56 PM
# TODO: Add comment 
# Author: raquel
###############################################################################

GetSamFile<-function(file, fieldNames=c("QNAME","FLAG",
				"RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL","TAGS"),fieldTypesNoTags = c("character()", 
				"integer()","character()","integer()","integer()","character()","character()","integer()",
				"integer()","character()", "character()", "character()") , toSep=";")
{
	afile<-system(paste("grep -v '^.*[[:space:]]4[[:space:]]\\*[[:space:]]0' " , file,  sep = ""),intern = TRUE)
#	afile<-scan(file,what=character(), sep = "\n", quote = "",comment.char="")
	select<-grep("^@SQ", afile)
	allLength<-vector("list",length(select) )
	temp<-afile[select]
	temp<-strsplit(gsub(":", "\t" , temp), split = "\t")
	startPos<-ListStartPositions(temp)
	temp<-unlist(temp)
	allLength$names<-temp[startPos+2]
	allLength$length<-temp[startPos+4]
	afile<-afile[-select] #to take the comments about the sequences on the top
	afile<-(strsplit(afile, split= "\t"))
	aLength<-length(fieldNames)
	aList<-vector("list", aLength)
	unlistedFile<-unlist(afile)
	for ( i in 1:(aLength-1))
	{
		aList[[i]]<-unlistedFile[(ListStartPositions(afile) +i -1)]
	}
	lengths<-unlist(lapply(afile, length))
	extraTags<-aLength:max(lengths)
	i<-1
	tempo<-unlistedFile[(ListStartPositions(afile)+ extraTags[i] -1)]
	tempo[which(lengths<extraTags[i])]<-"NA"
	aList[[aLength]]<-tempo
	for ( i in 2:length(extraTags))
	{
		tempo<-unlistedFile[(ListStartPositions(afile)+ extraTags[i] -1)]
		tempo[which(lengths<extraTags[i])]<-"NA"
		aList[[aLength]]<-paste(aList[[aLength]],tempo, sep = "_" )
	}
	intOnes<-which(fieldTypesNoTags == "integer()")
	charOnes<-which(fieldTypesNoTags == "character()")
	#to leave the tags out of the system
	for ( i in intOnes)
	{
		aList[[i]]<-as.integer(aList[[i]])
	}
	names(aList)<-fieldNames
	aList$names<-allLength$names
	aList$length<-allLength$length
	return(aList);	
#	Col Field Type
#	1 QNAME String
#	2 FLAG Int
#	3 RNAME String
#	4 POS Int
#	5 MAPQ Int
#	6 CIGAR String
#	7 RNEXT String
#	8 PNEXT Int
#	9 TLEN Int
#	10 SEQ String
#	11 QUAL String  ASCII of Phred-scaled base QUALity+33
#these are the CIGAR string meanings	
#M 0 alignment match (can be a sequence match or mismatch) 
#I 1 insertion to the reference
#D 2 deletion from the reference
#N 3 skipped region from the reference
#S 4 soft clipping (clipped sequences present in SEQ)
#H 5 hard clipping (clipped sequences NOT present in SEQ)
#P 6 padding (silent deletion from padded reference)
#= 7 sequence match
#X 8 sequence mismatch
}

SelectFirstReadsSam<-function(aSamFile, tolerated = 20, st="start", en="end", secondFastaFile)
{
	asc <- function(x) 
	{ 
		strtoi(charToRaw(x),16L) 
	}
	readQualitScore<-unlist(lapply(lapply(aSamFile$QUAL, asc), sum))
	if(length(names(table(readQualitScore)))> 1)
	{
		#for latter	
	}
	let<-unlist(lapply(strsplit(aSamFile$CIGAR, split ="\\d+", perl = TRUE), paste, collapse = ""))
	toKeep<-which(let == "MS" | let =="SM") #only selct for the 
	let<-let[toKeep]
	for ( i in 1:(length(aSamFile)-2))
	{
		aSamFile[[i]]<-aSamFile[[i]][toKeep]
	}
#	tolerance = 3*tolerated
#	toTake<-c(which(aSamFile$MAPQ > tolerance))
#	if(length(toTake) == 0)
#	{
#	}
#	else
#	{	
#		let<-let[-toTake]
#		for ( i in 1:(length(aSamFile)-2))
#		{
#			aSamFile[[i]]<-aSamFile[[i]][-toTake]
#		}
#	}
	num<-strsplit(aSamFile$CIGAR, split ="\\D+", perl = TRUE)
	matchsize<-NULL
	difference<-NULL
	sel<-which(let == "MS")
	matchsize[sel]<-as.numeric(unlist(num[sel])[ListStartPositions(num[sel])])
	difference[sel]<-as.numeric(unlist(num[sel])[ListStartPositions(num[sel])+1])
	sel<-which(let == "SM")
	matchsize[sel]<-as.numeric(unlist(num[sel])[ListStartPositions(num[sel])+1])
	difference[sel]<-as.numeric(unlist(num[sel])[ListStartPositions(num[sel])])
	TELength<-as.numeric(aSamFile$length[match(aSamFile$RNAME , aSamFile$names)])
	location<-rep("NA", length(aSamFile[[1]]))
	location[which(aSamFile$POS == 1)]<-st
	location[which((aSamFile$POS +matchsize -1) == TELength )]<-en
	toTake<-which(location =="NA")
	if(length(toTake) == 0)
	{
	}
	else
	{	
		let<-let[-toTake]
		location<-location[-toTake]
		matchsize<-matchsize[-toTake]
		difference<-difference[-toTake]
		for ( i in 1:(length(aSamFile)-2))
		{
			aSamFile[[i]]<-aSamFile[[i]][-toTake]
		}
	}
	id<-paste(aSamFile$QNAME,aSamFile$RNAME ,matchsize, difference, let, location, sep =",")
	myOut<-file(secondFastaFile, "w")
	cat(paste(">", id, "\n", aSamFile$SEQ, sep = ""), sep = "\n", file = myOut)
	close(myOut)
	firstReads<-vector("list")
	firstReads$id<-id
	firstReads$sequence<-aSamFile$SEQ
	return(firstReads)	
}

SelectSecondReadsSam<-function(otherSamFile, bedFile, tolerated = 20, st="start", en="end", strands=c("+", "-"), idSplit=",", pasteChar = ";")
{
	asc <- function(x) 
	{ 
		strtoi(charToRaw(x),16L) 
	}
	readQualitScore<-unlist(lapply(lapply(otherSamFile$QUAL, asc), sum))
	if(length(names(table(readQualitScore)))> 1)
	{
		#for latter	
	}
	let<-unlist(lapply(strsplit(otherSamFile$CIGAR, split ="\\d+", perl = TRUE), paste, collapse = ""))
	let<-gsub("H", "S",let )#it does not matter now the reads are there 
	toKeep<-which(let == "MS" | let =="SM") #only selct for the 
	if(length(toKeep) >=1)
	{
		let<-let[toKeep]
		for ( i in 1:(length(otherSamFile)-2))
		{
			otherSamFile[[i]]<-otherSamFile[[i]][toKeep]
		}
	}
	else
	{
		return(list(toKeep = NULL,readId=NULL ))
	}
#	tolerance = 3*tolerated
#	toTake<-c(which(otherSamFile$MAPQ > tolerance))
#	if(length(toTake) == 0)
#	{
#	}
#	else
#	{	
#		let<-let[-toTake]
#		for ( i in 1:(length(otherSamFile)-2))
#		{
#			otherSamFile[[i]]<-otherSamFile[[i]][-toTake]
#		}
#	}
	num<-strsplit(otherSamFile$CIGAR, split ="\\D+", perl = TRUE)
	matchsize<-NULL
	difference<-NULL
	sel<-which(let == "MS")
	matchsize[sel]<-as.numeric(unlist(num[sel])[ListStartPositions(num[sel])])
	difference[sel]<-as.numeric(unlist(num[sel])[ListStartPositions(num[sel])+1])
	sel<-which(let == "SM")
	matchsize[sel]<-as.numeric(unlist(num[sel])[ListStartPositions(num[sel])+1])
	difference[sel]<-as.numeric(unlist(num[sel])[ListStartPositions(num[sel])])
	tempId<-strsplit(otherSamFile$QNAME, split = idSplit)
	tempUnlist<-unlist(tempId)
	oldId<-vector("list", length=6)
	for( i in 0:5)
	{
		oldId[[i+1]]<-tempUnlist[ListStartPositions(tempId)+i]
	}
	names(oldId)<-c("realId","te", "match", "diff", "oldLet", "location")
	oldId$match<-as.numeric(oldId$match) 
	oldId$diff<-as.numeric(oldId$diff)
	toKeep<-which(oldId$match == difference & oldId$diff == matchsize)
	if(length(toKeep) >=1)
	{
		let<-let[toKeep]
		matchsize<-matchsize[toKeep]
		difference<-difference[toKeep]
		for ( i in 1:(length(otherSamFile)-2))
		{
			otherSamFile[[i]]<-otherSamFile[[i]][toKeep]
		}
		for( i in 1:6)
		{
			oldId[[i]]<-oldId[[i]][toKeep]
		}
	}
	else
	{
		return(list(toKeep = NULL,readId=NULL ))
	}
	strand<-rep(NA, length(matchsize))
	strand[which(oldId$oldLet == let)]<-strands[2]
	strand[which(oldId$oldLet != let)]<-strands[1]
	start<-otherSamFile$POS-1
	end<-otherSamFile$POS+matchsize-1
	toKeep<-paste(otherSamFile$RNAME,  start, end,strand, oldId$te, oldId$location, let, sep = pasteChar)
	if(length(toKeep) >=1)
	{
		if(missing(bedFile))
		{
			internalVar<-new.env();
			internalVar$toKeep<-toKeep
			internalVar$readId<-oldId$realId
			return(as.list(internalVar));	
		}
		else
		{
			myOutput<-file(bedFile, "w")
			aMat<-matrix(unlist(strsplit(toKeep, split = pasteChar)), nrow = length(toKeep), byrow = TRUE)
			chrom<-1
			start<-2
			end<-3
			bedID<-paste(oldId$realId, toKeep, sep = ";")
			cat(paste(aMat[,1],aMat[,2], aMat[,3],bedID,sep = "\t"), sep = "\n", file = myOutput)
			close(myOutput)
			internalVar<-new.env();
			internalVar$toKeep<-toKeep
			internalVar$readId<-oldId$realId
			return(as.list(internalVar));	
		}
	}
	else
	{
		return(list(toKeep = NULL,readId=NULL ))
	}
}

SelectSecondReadsSamOld<-function(otherSamFile,aSamFile, bedFile, tolerated = 20, st="start", en="end", strands=c("+", "-"), idSplit=",", pasteChar = ";")
{
	asc <- function(x) 
	{ 
		strtoi(charToRaw(x),16L) 
	}
	readQualitScore<-unlist(lapply(lapply(otherSamFile$QUAL, asc), sum))
	if(length(names(table(readQualitScore)))> 1)
	{
		#for latter	
	}
	let<-unlist(lapply(strsplit(otherSamFile$CIGAR, split ="\\d+", perl = TRUE), paste, collapse = ""))
	toKeep<-which(let == "M") #only selct for the 
	if(length(toKeep) >=1)
	{
		let<-let[toKeep]
		for ( i in 1:(length(otherSamFile)-2))
		{
			otherSamFile[[i]]<-otherSamFile[[i]][toKeep]
		}	
	}
	else
	{
		return(list(toKeep = NULL,readId=NULL ))
	}
#	tolerance = 3*tolerated
#	toTake<-c(which(otherSamFile$MAPQ > tolerance))
#	if(length(toTake) == 0)
#	{
#	}
#	else
#	{	
#		let<-let[-toTake]
#		for ( i in 1:(length(otherSamFile)-2))
#		{
#			otherSamFile[[i]]<-otherSamFile[[i]][-toTake]
#		}
#	}
	matchsize<-as.numeric(unlist(strsplit(otherSamFile$CIGAR, split ="\\D+", perl = TRUE)))
	tempId<-strsplit(otherSamFile$QNAME, split = idSplit)
	tempUnlist<-unlist(tempId)
	oldId<-vector("list", length=6)
	for( i in 0:5)
	{
		oldId[[i+1]]<-tempUnlist[ListStartPositions(tempId)+i]
	}
	names(oldId)<-c("realId","te", "match", "diff", "oldLet", "location")
	oldId$match<-as.numeric(oldId$match) 
	oldId$diff<-as.numeric(oldId$diff)
	aMatch<-match(oldId$realId, aSamFile$QNAME)
	strandAgreement<-rep("NA", length(matchsize))
	strandAgreement[which(otherSamFile$SEQ == aSamFile$SEQ[aMatch])]<-"+"
	strandAgreement[which(otherSamFile$SEQ != aSamFile$SEQ[aMatch])]<-"-"
	realLocation<-rep("NA", length(matchsize))
	sel<-which(strandAgreement == "+" & oldId$oldLet =="SM" )
	realLocation[sel]<-otherSamFile$POS[sel]+ oldId$diff[sel] 
	sel<-which(strandAgreement == "+" & oldId$oldLet =="MS" )
	realLocation[sel]<-otherSamFile$POS[sel]+ oldId$match[sel]
	sel<-which(strandAgreement == "-" & oldId$oldLet =="SM" )
	realLocation[sel]<-otherSamFile$POS[sel]+ oldId$match[sel] 
	sel<-which(strandAgreement == "-" & oldId$oldLet =="MS" )
	realLocation[sel]<-otherSamFile$POS[sel]+ oldId$diff[sel] 
	realLocation[which(oldId$location == "end")]<-as.numeric(realLocation[which(oldId$location == "end")])-1
	toKeep<-paste(otherSamFile$RNAME,  realLocation, strandAgreement, oldId$te, oldId$location, oldId$oldLet,  sep = pasteChar)
	if(length(toKeep) >=1)
	{
		if(missing(bedFile))
		{
			internalVar<-new.env();
			internalVar$toKeep<-toKeep
			internalVar$readId<-oldId$realId
			return(as.list(internalVar));	
		}
		else
		{
			myOutput<-file(bedFile, "w")
			aMat<-matrix(unlist(strsplit(toKeep, split = pasteChar)), nrow = length(toKeep), byrow = TRUE)
			chrom<-1
			start<-as.numeric(aMat[,2])
			end<-as.numeric(aMat[,2])
			sel<-which(strandAgreement == "+" & oldId$oldLet =="SM" )
			end[sel]<-as.numeric(aMat[,2])[sel]+ oldId$diff[sel] +1
			sel<-which(strandAgreement == "+" & oldId$oldLet =="MS" )
			start[sel]<-as.numeric(aMat[,2])[sel]- oldId$diff[sel] 
			sel<-which(strandAgreement == "-" & oldId$oldLet =="SM" )
			start[sel]<-as.numeric(aMat[,2])[sel]- oldId$diff[sel] 
			sel<-which(strandAgreement == "-" & oldId$oldLet =="MS" )
			end[sel]<-as.numeric(aMat[,2])[sel]+ oldId$diff[sel]
			bedID<-paste(oldId$realId, toKeep, sep = ";")
			cat(paste(aMat[,1],start,end, bedID,sep = "\t"), sep = "\n", file = myOutput)
			close(myOutput)
			internalVar<-new.env();
			internalVar$toKeep<-toKeep
			internalVar$readId<-oldId$realId
			return(as.list(internalVar));	
		}
	}
	else
	{
		return(list(toKeep = NULL,readId=NULL ))
	}
}

FinalProcessingSam<-function(secondReads, sample, tsd = 20, splitChar = ";", toPaste = ";")
{
	if(length(secondReads) == 0)
	{
		return(x=NULL)
	}
	aList<-strsplit(secondReads, split = splitChar)
	chrom<-1
	start<-2
	end<-3
	firstStrand<-4
	firstTeName<-5
	location<-6
	type<-7
	bigMatrix<-matrix(unlist(aList), byrow = TRUE, nrow = length(aList))
	aTeSet<-names(table(bigMatrix[,firstTeName]))
	readPairs<-NULL
	m<-1
	GetReadStrand<-function(tm, extraInfo, pos ="+_end", neg = "-_start")
	{
		aStrand<-NA
		if(length(tm) == 1)
		{
			if(extraInfo[tm] == pos)
			{
				aStrand<-"+"
			}
			else if(extraInfo[tm] == neg)
			{
				aStrand<-"-"
			}
			nReads<-1
		}
		else
		{
			tempTable<-table(extraInfo[tm])
			namesT<-names(tempTable)
			countTable<-c(tempTable, use.names=FALSE)
			tempName<-"NA"
			if(length(names(tempTable)) == 1)
			{
				tempName<-names(tempTable)
				nReads<-countTable
			}
			else
			{
				select<-which(countTable == max(countTable, na.rm = TRUE))
				if(length(select) == 1)
				{
					tempName<-names(tempTable)[select]
					nReads<-countTable[select]
				}
				else#unknown
				{
					nReads<-sum(countTable)
				}
			}
			if(tempName == pos)
			{
				aStrand<-"+"
			}
			else if(tempName == neg)
			{
				aStrand<-"-"
			}
		}
		return(list(strand = aStrand, reads = nReads))
	}
	for ( i in 1:length(aTeSet))
	{
		print( aTeSet[i])
		locations<-which(bigMatrix[,firstTeName] == aTeSet[i])
		tempo<-matrix(data = bigMatrix[locations,], ncol = type)
		matSize<-5
		smMattrix<-matrix(data = tempo[which(tempo[,type] == "SM"),c(chrom,start,firstStrand,firstTeName,location)], 
				ncol = matSize)
		msMattrix<-matrix(data = tempo[which(tempo[,type] == "MS"),c(chrom,end,firstStrand,firstTeName,location)],
				ncol = matSize)
		aChNameSet<-names(table(tempo[,chrom]))
		for ( j in 1:length(aChNameSet))
		{		
			tempSMMattrix<-matrix(data = smMattrix[which(smMattrix[,chrom] == aChNameSet[j]),], ncol = matSize)
			tempMSMattrix<-matrix(data = msMattrix[which(msMattrix[,chrom] == aChNameSet[j]),], ncol = matSize)
			starts<-as.numeric(tempSMMattrix[,2])
			ends<-as.numeric(tempMSMattrix[,2])
			extraInfoSt<-paste(tempSMMattrix[,3], tempSMMattrix[,5], sep = "_")
			extraInfoEn<-paste(tempMSMattrix[,3], tempMSMattrix[,5], sep = "_")
			tempstarts<-as.numeric(names(table(starts)))
			tempends<-as.numeric(names(table(ends)))		
			aMatrix<-matrix(data = NA, ncol = length(tempstarts), nrow = length(tempends), dimnames= list(tempends,tempstarts))
			if( length(tempstarts) == 0 |length(tempends) ==0 )
			{
				next
			}
			for ( s in 1:ncol(aMatrix))
			{
				aMatrix[,s]<-tempends - tempstarts[s]
			}
			goodStarts<-as.numeric(rep(colnames(aMatrix), each =nrow(aMatrix))[which(aMatrix >0 & aMatrix<=tsd)])
			goodEnds<-as.numeric(rep(rownames(aMatrix), ncol(aMatrix))[which(aMatrix >0 & aMatrix<=tsd)])
			matchStart<-match(starts, goodStarts)
			matchEnd<-match(ends, goodEnds)
			if(length(goodStarts) == 0)
			{
				next
			}
			for( s in 1:length(goodStarts))
			{
				tms<-which(matchStart == s)
				tme<-which(matchEnd == s)				
				st<-GetReadStrand(tms,extraInfoSt )
				en<-GetReadStrand(tme,extraInfoEn ,pos ="+_start", neg = "-_end")
				if(is.na(st$strand) == TRUE)
				{
					aStrand<-en$strand
				}
				else if (is.na(en$strand) == TRUE)
				{
					aStrand<-st$strand
				}
				else{
					if(en$strand == st$strand)
					{
						aStrand<-st$strand
					}
					else if (st$reads>en$reads )
					{
						aStrand<-st$strand
					}
					else if (st$reads<en$reads )
					{
						aStrand<-en$strand
					}
					else
					{
						aStrand<-NA
					}
				}
				readPairs[m]<-paste(aChNameSet[j],goodStarts[s], goodEnds[s], goodEnds[s]-goodStarts[s], aStrand, aTeSet[i], 
						sample, st$reads+en$reads, sep = toPaste)
				m<-m+1
			}		
		}
	}
	return(readPairs)
}

FinalProcessingSamOldTes<-function(secondReads, sample, aSamfile,minDist=100, maxDist = 1.5,splitChar = ";",toPaste = ";")
{
	if(length(secondReads) == 0)
	{
		return(x=NULL)
	}
	aList<-strsplit(secondReads, split = splitChar)
	chrom<-1
	start<-2
	strand<-3
	firstTeName<-4
	location<-5
	type<-6
	bigMatrix<-matrix(unlist(aList), byrow = TRUE, nrow = length(aList))
	aTeSet<-names(table(bigMatrix[,firstTeName]))
	readPairs<-NULL
	m<-1
	GetReadStrand<-function(tm, extraInfo, pos ="+_end", neg = "-_start")
	{
		aStrand<-NA
		if(length(tm) == 1)
		{
			if(extraInfo[tm] == pos)
			{
				aStrand<-"+"
			}
			else if(extraInfo[tm] == neg)
			{
				aStrand<-"-"
			}
			nReads<-1
		}
		else
		{
			tempTable<-table(extraInfo[tm])
			namesT<-names(tempTable)
			countTable<-c(tempTable, use.names=FALSE)
			tempName<-"NA"
			if(length(names(tempTable)) == 1)
			{
				tempName<-names(tempTable)
				nReads<-countTable
			}
			else
			{
				print(countTable)
				select<-which(countTable == max(countTable, na.rm = TRUE))
				if(length(select) == 1)
				{
					tempName<-names(tempTable)[select]
					nReads<-countTable[select]
				}
				else#unknown
				{
					nReads<-sum(countTable)
				}
			}
			if(tempName == pos)
			{
				aStrand<-"+"
			}
			else if(tempName == neg)
			{
				aStrand<-"-"
			}
		}
		return(list(strand = aStrand, reads = nReads))
	}
	TestDistances<-function(test, against,aDist )
	{
		tempo<-as.numeric(names(table(test)))
		for(s in 1:length(tempo))
		{
			aMatch<-which(is.na(match(test, tempo[s])) == FALSE)
			if(length(aMatch)>1)
			{
				tmin<-min(abs(test[aMatch]-against[aMatch]))
				if(tmin >= aDist)
				{
					tot<-which(abs(test[aMatch]-against[aMatch]) != tmin)
					if(length(tot)>0)
					{
						test<-test[-aMatch[tot]]
						against<-against[-aMatch[tot]]
					}	
					tempo<-as.numeric(names(table(test)))
				}
				else
				{
					while((tmin < aDist) & length(aMatch)>1)
					{
						aMatch<-which(is.na(match(test, tempo[s])) == FALSE)
						tot<-which(abs(test[aMatch]-against[aMatch]) == tmin)
						test<-test[-aMatch[tot]]
						against<-against[-aMatch[tot]]
						tempo<-as.numeric(names(table(test)))
						aMatch<-which(is.na(match(test, tempo[s])) == FALSE)
						tmin<-min(abs(test[aMatch]-against[aMatch]))
					}
				}
			}					
		}
		return(list("test" = test, "against" = against))
	}
	for ( i in 1:length(aTeSet))
	{
		print( aTeSet[i])
		tempMax<-maxDist*as.numeric(aSamFile$length[match(aTeSet[i], aSamFile$names)])
		locations<-which(bigMatrix[,firstTeName] == aTeSet[i])
		tempo<-matrix(data = bigMatrix[locations,], ncol = type)
		matSize<-5
		sMattrix<-matrix(data = tempo[which(tempo[,location] == "start"),c(chrom,start,strand,firstTeName,type)],	ncol = matSize)
		eMattrix<-matrix(data = tempo[which(tempo[,location] == "end"),c(chrom,start,strand,firstTeName,type)],	ncol = matSize)
		aChNameSet<-names(table(tempo[,chrom]))
		for ( j in 1:length(aChNameSet))
		{		
			tempSMattrix<-matrix(data = sMattrix[which(sMattrix[,chrom] == aChNameSet[j]),], ncol = matSize)
			tempEMattrix<-matrix(data = eMattrix[which(eMattrix[,chrom] == aChNameSet[j]),], ncol = matSize)
			starts<-as.numeric(tempSMattrix[,2])
			ends<-as.numeric(tempEMattrix[,2])
			extraInfoSt<-tempSMattrix[,3]
			extraInfoEn<-tempEMattrix[,3]
			tempstarts<-as.numeric(names(table(starts)))
			tempends<-as.numeric(names(table(ends)))		
			aMatrix<-matrix(data = NA, ncol = length(tempstarts), nrow = length(tempends), dimnames= list(tempends,tempstarts))
			if( length(tempstarts) == 0 |length(tempends) ==0 )
			{
				next;
			}
			for ( s in 1:ncol(aMatrix))
			{
				aMatrix[,s]<-tempends - tempstarts[s]
			}
			goodStarts<-as.numeric(rep(colnames(aMatrix), each =nrow(aMatrix))[which(aMatrix >=minDist & aMatrix<=tempMax)])
			goodEnds<-as.numeric(rep(rownames(aMatrix), ncol(aMatrix))[which(aMatrix >=minDist & aMatrix<=tempMax)])
			if(length(goodStarts) == 0)
			{
				next;
			}
			if(length(goodEnds) != length(names(table(goodEnds))))
			{
				t<-TestDistances(goodEnds, goodStarts, minDist*8)
				goodEnds<-t$test
				goodStarts<-t$against
			}
			if(length(goodStarts) != length(names(table(goodStarts))))
			{
				t<-TestDistances(goodStarts,goodEnds, minDist*8)
				goodStarts<-t$test
				goodEnds<-t$against
			}
			matchStart<-match(starts, goodStarts)
			matchEnd<-match(ends, goodEnds)
			for( s in 1:length(goodStarts))
			{
				tms<-which(matchStart == s)
				tme<-which(matchEnd == s)				
				st<-GetReadStrand(tms,extraInfoSt,pos ="+", neg = "-" )
				en<-GetReadStrand(tme,extraInfoEn ,pos ="+", neg = "-")
				if(is.na(st$strand) == TRUE)
				{
					aStrand<-en$strand
				}
				else if (is.na(en$strand) == TRUE)
				{
					aStrand<-st$strand
				}
				else{
					if(en$strand == st$strand)
					{
						aStrand<-st$strand
					}
					else if (st$reads>en$reads )
					{
						aStrand<-st$strand
					}
					else if (st$reads<en$reads )
					{
						aStrand<-en$strand
					}
					else
					{
						aStrand<-NA
					}
				}
				readPairs[m]<-paste(aChNameSet[j],goodStarts[s], goodEnds[s], goodEnds[s]-goodStarts[s], aStrand, aTeSet[i], 
						sample, st$reads+en$reads, sep = toPaste)
				m<-m+1
			}		
		}
	}
	return(readPairs)
}

ListStartPositions<-function(aList)
{
	return(as.numeric(which(sequence(unlist(lapply(aList,length))) == 1)))
}

GetFasta<-function(file, aNameLocation = 1, charSplit = " ", sizeLocation = 7, sizeSplit = "=", allChroms = "yes")
{
	aFile<-scan(file, what = "character", sep = "\n")
	locations<-grep(">" , aFile)
	nChar<-length(aFile)
	tempo<-strsplit(aFile[locations], charSplit)
	aNameLocation<-aNameLocation-1
	seqNames<-unlist(tempo)[ListStartPositions(tempo)+aNameLocation]
	nChrom<-length(seqNames)
	seqNames<-unlist(strsplit(seqNames, split= ">"))[seq(2, nChrom*2, 2)]
	sequences<-NULL
	if(allChroms == "yes")
	{
		for( i in 1:(nChrom-1))
		{
			sequences[i]<-paste(aFile[seq(locations[i]+1, locations[i+1]-1, 1)], sep = "", collapse = "")
		}
		sequences[nChrom]<-paste(aFile[seq(locations[nChrom]+1, nChar, 1)], sep = "", collapse = "")
	}
	else
	{
		sequences[1]<-paste(aFile[seq(2, length(aFile), 1)], sep = "", collapse = "")
	}
	internalVar<-new.env();
	if(is.numeric(sizeLocation) == TRUE)
	{
		sizeLocation<-sizeLocation-1
		seqSize<-unlist(strsplit(unlist(tempo)[ListStartPositions(tempo)+sizeLocation], 
						split = sizeSplit))[seq(2, nChrom*2, 2)]
		if(is.numeric(seqSize))
		{
			internalVar$seqSize = seqSize;
		}
		else
		{
			internalVar$seqSize = nchar(sequences);
		}
		internalVar$sequences = sequences;
		internalVar$seqNames = seqNames;
		
		return(as.list(internalVar))	
	}
	else{
		internalVar$sequences = sequences;
		internalVar$seqNames = seqNames;
		internalVar$seqSize = nchar(sequences);
		return(as.list(internalVar))				
	}	
}

GetSequences<-function (fastaFile, chrom, start,winLeft, winRigth, strand,  complement = 1, posPairs  = c("A","C","G","T"), negPairs = c("T","G", "C", "A"))
{
	start<-as.vector(start[which (strand != "NA")])
	chrom<-as.vector(chrom[which(strand != "NA")])
	strand<-strand[which(strand != "NA")]
	end<-start
	start [strand == 1 | strand == "+"]<-start[strand == 1 | strand == "+"] - winLeft 
	end[strand == 1 | strand == "+"]<-end[strand == 1 | strand == "+"] + winRigth 
	start [strand == 0 | strand == "-"]<-start[strand == 0 | strand == "-"] - winRigth
	end[strand == 0 | strand == "-"]<-end[strand == 0 | strand == "-"] + winLeft 
	myLength<- end[1]-start[1] + 1 
	print (length (start))
	toTake<-which(start <= 0)
	if(length(toTake) > 0)
	{
		start<-start[- toTake]
		end<-end[- toTake]
		chrom<-chrom[- toTake]
		strand<-strand[- toTake]	
		print (paste("left ",length(toTake), " out because they were before the start of the chrommosome"))
	}
	mySequences<-NULL
	myStrands<-NULL
	myLocations<-NULL
	n<-1
	chrNames<-names(table(chrom))
	for ( i in 1:length (chrNames)){
		selected<-which(chrom == chrNames[i])
		startSelected<-start[selected]
		endSelected<-end[selected]
		strandSelected<-strand[selected]
		chrLocation<-which(fastaFile$seqNames == chrNames[i])
		toTake<-which(endSelected > fastaFile$seqSize[chrLocation])
		if(length(toTake) > 0)
		{
			print (paste("left ",length(toTake), " out because they were before after the end of chrommosome " ,
							chrNames[i]))
			startSelected<-startSelected[- toTake]
			endSelected<-endSelected[- toTake]
			strandSelected<- strandSelected[- toTake]
		}
		tempLocations<-paste(chrNames[i],startSelected,endSelected,strandSelected, sep = ";")
		for(j in 1:length(startSelected))
		{
			mySequences[n]<-substr(fastaFile$sequences[chrLocation],startSelected[j], endSelected[j])
			myStrands[n]<-strandSelected[j]
			myLocations[n]<-tempLocations[j]
			n<-n+1
		}
	}
	mySequences<-strsplit(mySequences, split = "")
	mySequences<-matrix(data = unlist(mySequences),nrow = length(mySequences), ncol = length(mySequences[[1]]),
			byrow = TRUE)
	if(complement == 1)
	{
		print("getting the reverse complement\n")
		selected<-which(myStrands == 0 | myStrands == "-" )
		tempSequences<-matrix(data = mySequences[selected,], ncol = ncol(mySequences))
		reversedSequences<-tempSequences
		myLength<-ncol(tempSequences)
		for ( i in 1:length(posPairs)){
			reversedSequences[which(tempSequences == posPairs[i])]<-negPairs[i]
		}
		reversedSequences<-reversedSequences[,seq(myLength, 1,-1)]
		mySequences[selected,]<-reversedSequences
	}
	myLength<-ncol(mySequences)
	freqMatrix<-matrix(data = NA, ncol = myLength, nrow = length(posPairs))
	for ( i in 1:myLength)
	{
		for (j in 1:length(posPairs))
		{
			freqMatrix[j,i]<-length(which(mySequences[,i] == posPairs[j]))
		}
	}
	percentMatrix<-freqMatrix/nrow(mySequences)
	internalVar<-new.env();
	internalVar$sequences = mySequences;
	internalVar$strands = myStrands;
	internalVar$id = myLocations;
	internalVar$matrix = percentMatrix
	return(as.list(internalVar))	
}


myA<-data.frame(V1 = c(0.715058211295326,0.9,0.597367982119625,0.398015144944458,0.086976682122961,0.27071755012176,
0.329562664709611,0.65741401741335,0.608176268472496,0.38000133435634,0.494088801414418,
0.608176268472496),
V2 =  c(0,0,1,1,0,0,0.201646666800031,0.201646666800031,0.37311492088169,0.37311492088169,0.764059682383637,
0.37311492088169))

myC<-data.frame(V1 = c(0.895085462853275, 0.893580132002031, 0.891808258588593, 0.889757150110002, 0.887406498561517, 
0.884735995938399, 0.881730411237096, 0.87837197495346, 0.874642917583347, 0.870522931122017,
0.865996784565916, 0.861046708410899, 0.855652394652225, 0.849801150786935, 0.843467591809105,
0.836639025215772, 0.830120155694703, 0.823365205618548, 0.816369097986123, 0.80914452530039,
0.801688949060755, 0.794007446268404, 0.78610255542393, 0.777979353528516, 0.769640379082755,
0.761090709087832, 0.752327805043155, 0.743361820951092, 0.734190218311051, 0.724823151125402, 
0.715258080893552, 0.705500084616686, 0.695554239295989, 0.68542054493146, 0.675106617024877,
0.66461245557624, 0.653943137586732, 0.643098663056355, 0.632089185987477, 0.620909629378913,
0.609570147233034, 0.598073278050432, 0.586419021831105, 0.574609917075647, 0.562656117786427,
0.550552546962261, 0.538309358605517, 0.525926552716195, 0.502673887290574, 0.479992384498223,
0.457887121340328, 0.436373328820443, 0.415461160940938, 0.39516077170418, 0.375484853613133,
0.356443560670164, 0.338049585378236, 0.320310543239127, 0.303241665256389, 0.286853105432391,
0.271152479268912, 0.256155017769504, 0.241870874936537, 0.228312743272973, 0.215488238280589,
0.203410052462346, 0.192090878321205, 0.181538331358944, 0.171767642579117, 0.162788965984092, 
0.154609917075647, 0.147245726857336, 0.140706549331528, 0.135, 0.130143848366898, 0.126145709933999,
0.123015738703672, 0.120766627178880, 0.119408529361990, 0.118954137755965, 0.119405990861398,
0.120751396175326, 0.122982738195972, 0.126084785919783, 0.130049923844982, 0.134865459468607,
0.140516161787104, 0.146996953799289, 0.154292604501608, 0.162392959891691, 0.171285327466576, 
0.180959553223896, 0.191405483161279, 0.202607886275173, 0.214559147063801, 0.227246573024200, 
0.24065747165341, 0.254781688949061, 0.269609070908783, 0.285126925029616, 0.301322558808597,
0.31818835674395, 0.335706549331528, 0.353872059570147, 0.372672194956845, 0.392094262988661,
0.412125571162633, 0.432758503976984, 0.453977830428161, 0.475775935014385, 0.498137586732104,
0.521052631578947, 0.541512946353021, 0.561531562024031, 0.581100863090201, 0.600205618547978,
0.618835674394991, 0.636980876628871, 0.654625994246065, 0.671763411744796, 0.68837789812151, 
0.704461837874429, 0.72, 0.734982230495854, 0.749395836859029, 0.763230665087155, 0.776476561177864, 
0.789118294127602, 0.801145709933999, 0.812548654594686, 0.823314435606702, 0.833430360467084, 
0.842886275173464, 0.851672025723473, 0.859772381113555, 0.867177187341344, 0.873876290404468, 
0.879854459299374, 0.885106617024877, 0.889614994076832, 0.893369436452868, 0.896359790150618,
0.898573362667118, 0.9, 0.722101878490438, 0.720896090709088, 0.719337451345405, 0.717428498899983,
0.715171771873413, 0.712572347266881, 0.70963530208157, 0.706363174818074, 0.702758503976984, 
0.698823828058893, 0.694566762565578, 0.689987307497038, 0.68509053985446, 0.679878998138433,
0.674357759350144, 0.668529361990185, 0.662396344559147, 0.655963784058216, 0.649236757488577,
0.642215264850229, 0.634904383144356, 0.627309189372144, 0.619429683533593, 0.611273481130479,
0.602843120663395, 0.594141140632933, 0.585170079539685, 0.575935014384837, 0.566441022169572, 
0.556688102893891, 0.54668387205957, 0.53642832966661, 0.525926552716195, 0.51284819766458, 
0.500115078693518, 0.487734811304789, 0.475712472499577, 0.464048062277881, 0.452749196141479, 
0.441820951091555, 0.431265865628702, 0.421089016754104, 0.411295481468946, 0.40188779827382, 
0.392873582670503, 0.384255373159587, 0.376038246742257, 0.368224741919106, 0.360822474191911,
0.353833982061262, 0.347264342528347, 0.341116094093755, 0.335396852259266, 0.330106617024877,
0.32525554239296, 0.320846166864106, 0.316881028938907, 0.313365205618548, 0.310303773904214, 
0.307701810797089, 0.305561854797766, 0.303888982907429, 0.302690810627856, 0.301967337959046,
0.301726180402775, 0.301977491961415, 0.302731426637333, 0.303982907429345, 0.305724318835674,
0.307950583855136, 0.310656625486546, 0.313837366728719, 0.317490269081063, 0.321605178541209,
0.326179556608563, 0.331205787781351, 0.336681333558978, 0.342598578439668, 0.348954983922830,
0.355740396006092, 0.362954814689457, 0.370588085970553, 0.378640209849382, 0.387101032323574,
0.395968014892537, 0.405233542054493, 0.414895075308851, 0.424944999153833, 0.435378236588255,
0.446189710610932, 0.45737434422068, 0.468927060416314, 0.480842782196649, 0.493113894059909, 
0.5057378575055, 0.518707057031647, 0.532018954137756, 0.54376967337959, 0.55525131155864, 
0.566453714672534, 0.577369267219496, 0.587982738195972, 0.598289050600778, 0.608270434929768,
0.617921814181757, 0.627230495853783, 0.636183787442884, 0.644776611947876, 0.652991199864613,
0.66082247419191, 0.668255203926214, 0.675281773565747, 0.68189202910814, 0.685435775935015, 
0.688766288712134, 0.691901336943645, 0.694858690133695, 0.69764850228465, 0.700291081401252,
0.702801658487054, 0.705195464545609, 0.707490269081063, 0.709701303096971, 0.711843797596886, 
0.713932983584363, 0.71598916906414, 0.718025046539178, 0.720055847013031, 0.722101878490438, 
0.896344559147064),
V2 =  c(0.676111682241716, 0.690569618699877, 0.70441127948102, 0.717688475372151, 0.730442109626172, 
0.742718539263032, 0.75456957506974, 0.766033393415667, 0.777164531971346, 0.787998440222623, 
0.79858965584003, 0.808981808959994, 0.819221256602467, 0.829354355787401, 0.839422009767697,
0.84947602933036, 0.858431114831792, 0.867129873281722, 0.87557230468015, 0.883755682143548, 
0.891677278788391, 0.899337094614678, 0.906729675855355, 0.91385774939395, 0.920718588346936,
0.92730673894726, 0.933624928078447, 0.939667701973446, 0.945435060632255, 0.95092427717135, 
0.95613262470720, 0.96106010323981, 0.96570398588565, 0.970061545761196, 0.974132782866445, 
0.977912243434346, 0.981402654348425, 0.984598561841628, 0.987499965913956, 0.990101412798355,
0.992405629378352, 0.99440988877042, 0.996108737207508, 0.99750490157314, 0.998592928100262, 
0.999375543672403, 0.999841840755456, 1, 0.999438261993515, 0.997763955508168, 0.994987988078065,
0.991129447887892, 0.98619651558823, 0.980208279363764, 0.973172919865074, 0.965109525276847,
0.95602900313319, 0.94594771473526, 0.934876567617167, 0.92283192308007, 0.909824688658073, 
0.895873952535865, 0.880987895364025, 0.865182878443713, 0.848475263076088, 0.830873229911731, 
0.812395867135327, 0.79305680916451, 0.772866963533387, 0.75183996465959, 0.729992173844279,
0.707337225505087, 0.683888754059648, 0.659660393925594, 0.634665779520559, 0.608918545262176,
0.582435052451605, 0.555226208622951, 0.527308375077375, 0.498689732465457, 0.470313782487409,
0.442616826507489, 0.415617952710386, 0.389328068630205, 0.363760808684578, 0.338935261058194,
0.314862333285158, 0.29155838666663, 0.269037055620243, 0.247309247680104, 0.226396777914425, 
0.206305100090260, 0.187056029275821, 0.168660473005216, 0.151134792579604, 0.134489895533092,
0.118742143166839, 0.103907896782005, 0.0900007907962225, 0.0770317327435993, 0.0650198108088209,
0.0539732056424672, 0.0439137323127517, 0.0348495714702536, 0.0267998112996599, 0.0197753593350766,
0.0137925768776638, 0.0088678252285809, 0.00500928503840805, 0.00223604449183151, 0.000561738006484404,
0, 0.000400851878413571, 0.00159795374660152, 0.00357767118693057, 0.00633182354882064, 
0.00985223018169208, 0.0141252566679118, 0.0191399954733733, 0.0248909928310232, 0.0313618874397016,
0.0385472255323557, 0.0464333726913522, 0.055012148266111, 0.064272644722526, 0.0742039545264902,
0.084797897027424, 0.0960408378076946, 0.107924596216722, 0.120438264720399, 0.133573662668146, 
0.147317155642331, 0.161657836108846, 0.176590250300639, 0.192098036916549, 0.208175742189524, 
0.224812458585456, 0.241997278570240, 0.259716567726243, 0.27796487228641, 0.296731284716636,
0.316002170599287, 0.335769349283784, 0.35602191323602, 0.35602191323602, 0.344762611154590, 
0.333814173795195, 0.323179328041361, 0.312860800776616, 0.302861318884486, 0.293183609248498,
0.283830398752178, 0.274809868046106, 0.266116563363229, 0.257758665354127, 0.249741627785852, 
0.242062723774879, 0.234727407088261, 0.227741131493051, 0.221106623872774, 0.214823884227433,
0.208898366324079, 0.203330070162713, 0.198127176393915, 0.193289685017684, 0.188823049801074,
0.184727270744085, 0.181005074730243, 0.177661915526602, 0.174703246900215, 0.172126341967555, 
0.169939381379203, 0.168139638251631, 0.166738020118947, 0.165731800097622, 0.165123705071185, 
0.164921915690215, 0.165281864315729, 0.166353529541692, 0.168134184484578, 0.170615648493806,
0.173789740918796, 0.177648281108969, 0.182185815297271, 0.187396889716649, 0.193276050600051,
0.199809663529842, 0.206995001622496, 0.214826611110960, 0.223293584461127, 0.232393194789471,
0.242114534561885, 0.252454876894843, 0.263403314254238, 0.274954392873017, 0.287102658984127,
0.299837205053460, 0.313155304197492, 0.327048775765641, 0.341509439107327, 0.356531840455499, 
0.372110526043101, 0.388231861452502, 0.404895846683700, 0.422094301086118, 0.439816317125646, 
0.458059167918761, 0.47681467281488, 0.496074651163425, 0.515580049029366, 0.534578246559355,
0.553066516869865, 0.571031225543263, 0.588469645696023, 0.605370869794038, 0.621729444070255,
0.637537187874094, 0.652788647438502, 0.667472915229372, 0.68158453747965, 0.695115333538758,
0.708059849639642, 0.720409905131722, 0.732154592480891, 0.743291184803623, 0.753811501449339,
0.76370463488393, 0.772967858223872, 0.781590263935056, 0.789566398250431, 0.79688535363589,
0.803547130091432, 0.809535366315898, 0.814847335425762, 0.819474856770443, 0.823412476582888,
0.82664928732899, 0.829182562125224, 0.830998666553955, 0.832094873731658, 0.832460276124226, 
0.832184860888037, 0.83135861517947, 0.829984265882051, 0.828064539879308, 0.82559943717124, 
0.822591684641374, 0.81904128228971, 0.8149536838833, 0.810328889422146, 0.8051723526733,
0.79947861986971, 0.793255871661954, 0.786504108050033, 0.779223329033947, 0.771418988380749,
0.763088359206913, 0.758190876393097, 0.753279759161647, 0.748322284910245, 0.743272096618937,
0.738091017918352, 0.732740872439115, 0.72718075692833, 0.721369768133094, 0.715264275916983,
0.70882883079415, 0.702022529511697, 0.694801741933197, 0.687128291689277, 0.678964002410565,
0.670265243960635, 0.660993839970113, 0.660993839970113))

myG<-data.frame(V1 = c( 0.543588526479962, 0.543588526479962, 0.746206008001627, 0.745201566420289, 0.744016579643317,
0.742638333220316, 0.741056655590968, 0.739258832304876, 0.737237234691802, 0.73497914830135, 
0.732474401573201, 0.72971282294704, 0.726681697972469, 0.72337085508917, 0.719772665626907, 
0.715871872245202, 0.711663389163898, 0.707129416152438, 0.702264867430664, 0.693758900115278,
0.684739269003865, 0.67523140299722, 0.665265816776293, 0.654870482131959, 0.644078456635248,
0.632910083406795, 0.62140096290771, 0.609576524038787, 0.5974647385909, 0.585098664135078,
0.572501186682037, 0.559702820912728, 0.546734081508103, 0.533622940259036, 0.520397368956398,
0.506675934088289, 0.493206245338035, 0.48000610293619, 0.467083135552994, 0.454452600528921,
0.44212466942429, 0.430109513799417, 0.418419848104699, 0.407070929680613, 0.396067844307317,
0.385428392215366, 0.375162744965078, 0.365283617006849, 0.355798637010917, 0.346725605207839, 
0.338072150267851, 0.329850986641351, 0.322074828778735, 0.31475384824032, 0.307903302366583,
0.301530819827762, 0.295651657964332, 0.290275988336611, 0.285416525394996, 0.281083440699803,
0.277286905811352, 0.274044720960195, 0.271367057706652, 0.26926154472096, 0.267743439343595,
0.266825456024954, 0.266515223435275, 0.266789855563843, 0.26760866616939, 0.268966569471757, 
0.270855936800705, 0.273274225266156, 0.276211263307791, 0.279661965145453, 0.283623787889062,
0.288084017088221, 0.29304519563301, 0.298494609073032, 0.304427171628128, 0.310840340408219, 
0.317723943852987, 0.325072896182274, 0.332884654506001, 0.34114904726385, 0.349860988675663,
0.359017935851360, 0.368607174340544, 0.378628704143215, 0.389074896589137, 0.399938123008069,
0.411210754729776, 0.422892791754255, 0.434974062521191, 0.447446938360344, 0.460308876381637,
0.47355479080491, 0.487171967179765, 0.501162948396284, 0.515515020004069, 0.527298772631722, 
0.538851122262155, 0.550161897335051, 0.561228554960331, 0.572046009357835, 0.582606631857327,
0.592905336678646, 0.602939580931715, 0.612701735946294, 0.622186715942226, 0.631386892249271,
0.640299721977352, 0.648922662236387, 0.657242998575981, 0.665260730996135, 0.672970773716688,
0.677636977012274, 0.682038719739608, 0.686198887909406, 0.690135281752221, 0.693873330168848,
0.697430833389842, 0.700830677425917, 0.704095748287787, 0.707246389096088, 0.710302942971452,
0.713288295924595, 0.71622533396623, 0.719131857326914, 0.722033294907439, 0.72494744693836, 
0.727897199430393, 0.9, 0.897568997084153, 0.894408184715536, 0.890525191564386, 0.885930189191022, 
0.88063589204584, 0.874655014579237, 0.867997728351529, 0.860674204923035, 0.852697158744151, 
0.844076761375195, 0.834820641486404, 0.824949142198413, 0.814464806401302, 0.803382891435546,
0.791713568861463, 0.779467010239371, 0.766655930019665, 0.753290499762664, 0.739383433918763, 
0.724944904048281, 0.709987624601614, 0.694519224249, 0.678552417440835, 0.662099918627517, 
0.645171899369363, 0.62777853122669, 0.609932528649895, 0.591644063199295, 0.572925849325286,
0.553785515698108, 0.534238319658236, 0.514294432765986, 0.502337763612938, 0.490505696073778,
0.478798230148505, 0.467215365837119, 0.455759646029701, 0.444436156506408, 0.433244897267241,
0.422188411202278, 0.411266698311521, 0.400487387265206, 0.389845392283176, 0.379348342035668,
0.368996236522683, 0.358791618634299, 0.348737031260595, 0.338835017291653, 0.32908557672747, 
0.319491252458127, 0.310054587373703, 0.300780667254357, 0.291666949210009, 0.282718519020818,
0.273937919576863, 0.265322607988065, 0.256880212924663, 0.248610734386655, 0.240514172374042,
0.232595612666983, 0.224860141045636, 0.217300128839764, 0.209928290499763, 0.202739540245474,
0.195738963857056, 0.188929104224588, 0.182309961348071, 0.175886621007663, 0.169656540313284,
0.163627347935173, 0.157796500983251, 0.152171628127755, 0.146747643588527, 0.141534718925883,
0.136527768359666, 0.131731877670035, 0.127149589747067, 0.122783447480843, 0.118633450871364,
0.114704685698786, 0.110997151963111, 0.107513392554418, 0.104253407472706, 0.101224825388215,
0.0984251034108633, 0.0958567844307318, 0.0935224113378993, 0.0914270699125245, 0.0895707601546077,
0.0879534820641487, 0.0865777785312266, 0.0854487353360006, 0.0845688953685496, 0.0839357157387942,
0.0835542822268937, 0.0834271377229266, 0.083556825116973, 0.0839382586288737, 0.0845739811487082,
0.085458906896318, 0.0865930358717028, 0.0879687394046247, 0.0895911032752424, 0.0914474130331592,
0.0935452973486133, 0.0958771275513665, 0.0984403607513392, 0.101232454058452, 0.104253407472706,
0.107495592323862, 0.110961551502000, 0.114648742117041, 0.118549535498745, 0.122666474537194, 
0.126994473452228, 0.131530989353767, 0.136276022241812, 0.141221943446125, 0.146371295856784,
0.151718993693633, 0.157265036956669, 0.163004339865735, 0.168934359530752, 0.175052553061640,
0.181358920458398, 0.18784583305079, 0.194515833728894, 0.201366379602631, 0.208389842001763, 
0.215588763816369, 0.222960602156371, 0.230497728351529, 0.238200142401844, 0.246070387197396,
0.254098291177867, 0.262283854343256, 0.270627076693565, 0.279122872448634, 0.287768698718383,
0.296562012612735, 0.305502814131688, 0.314586017495084, 0.323809079812843, 0.333172001084966,
0.342667152641215, 0.352297077371669, 0.362056689496169, 0.371945989014715, 0.381959890147148,
0.392095850003391, 0.402351325693361, 0.412726317217061, 0.423215738794331, 0.433817047535092,
0.444527700549264, 0.455347697836848, 0.466271953617685, 0.477297925001695, 0.4884230690988,
0.499647385908998, 0.51074710110531, 0.521623041974639, 0.532282837187224, 0.542723943852987,
0.552959076422323, 0.562990777785312, 0.572821590832034, 0.582456601342646, 0.591903437987387,
0.601164643656337, 0.610247847019733, 0.619155590967654, 0.62789296128026, 0.63646504373771, 
0.64487946701024, 0.65313623109785, 0.66124296467078, 0.669207296399268, 0.677026683393232, 
0.684713840103072, 0.692268766528786, 0.699699091340611, 0.707009900318709, 0.714203736353156,
0.721285685224114, 0.728263375601817, 0.735139350376348, 0.741918695327863, 0.748609039126602,
0.755210381772564, 0.761730351935987, 0.768176578287109, 0.790147148572591, 0.9, 0.9),
V2 =  c(0.54319247161996, 0.379582186906051, 0.379582186906051, 0.369759952443151, 0.360395834412725, 
0.351457110212451, 0.342911057240012, 0.334727679776614, 0.32687152833641, 0.319309880317082,
0.312010013116310, 0.304941931015300, 0.298070184528208, 0.291362051052713, 0.284787534870023,
0.278313913377818, 0.271905737090252, 0.265533010288532, 0.259163010370338, 0.248642693724623,
0.238689568852446, 0.229319997054966, 0.220550339633343, 0.212391504121684, 0.204865305588202,
0.197979924683477, 0.191754449592195, 0.186205241615515, 0.181343208287544, 0.177190164676496,
0.173757018316477, 0.171060130508645, 0.169115862554163, 0.167935121987134, 0.167539723875774, 
0.167943302637714, 0.169145858272956, 0.171131029480338, 0.173885181842228, 0.177391954057466, 
0.181637711708420, 0.186608820377455, 0.192286191879886, 0.198658918681606, 0.205710639481456, 
0.213427719861802, 0.221793798521484, 0.230792514159343, 0.240412959241272, 0.250638772466112, 
0.261453592532702, 0.27284378502341, 0.284795715520603, 0.297290295839594, 0.310316618446276, 
0.323861048923017, 0.337904499085131, 0.352433334514983, 0.367433920794941, 0.382889896623845,
0.39879035446759, 0.415113479257960, 0.431851090344378, 0.448984099542156, 0.466501599317188,
0.484384501484788, 0.502619171627322, 0.521644637992577, 0.54018744597362, 0.55823941491987,
0.575792364180749, 0.592838113105675, 0.609368481044069, 0.625375287345352, 0.640850351358942, 
0.655788219317788, 0.670177983687783, 0.684016917585399, 0.697288659709478, 0.709993210060019,
0.722116934219388, 0.733657105304061, 0.74460281577993, 0.754943158112887, 0.764678132302935,
0.773791377048912, 0.782280165467292, 0.790136316907496, 0.797348923835416, 0.803915259367527,
0.809821689086194, 0.815062759224365, 0.81963028913146, 0.823518825040426, 0.826717459417156,
0.82921801161107, 0.831015027855115, 0.832097600615185, 0.832460276124226, 0.83223121790799, 
0.83154404325928, 0.830398752178098, 0.828800798431496, 0.826750182019475, 0.824246902942035, 
0.821296414966227, 0.817898718092054, 0.814056539203041, 0.809772605182715, 0.805046916031075, 
0.79988219863165, 0.794278452984437, 0.788241132856492, 0.781770238247814, 0.774868496041928, 
0.770459125379378, 0.76603884718272, 0.76156403131553, 0.756993774524909, 0.752287173557956, 
0.747397871394719, 0.742287691665826, 0.736910277351324, 0.731230178965366, 0.725198312604474, 
0.718773775015748, 0.711918389829815, 0.70458852691025, 0.696740556120626, 0.688330847324518,
0.679318497269026, 0.679318497269026, 0.697514991042188, 0.715275183451089, 0.732599074495731,
0.749467575991427, 0.765880687938176, 0.781819322151293, 0.797283478630777, 0.812254069191942, 
0.826731093834789, 0.84069546437463, 0.85414445392794, 0.867064428077084, 0.87944993305501, 
0.891287334444084, 0.902571178477254, 0.91328510385336, 0.923426383688873, 0.93298411044969, 
0.941944649718177, 0.950302547727279, 0.95804689694289, 0.965169516714432, 0.97165677262427,
0.977503210905353, 0.982697924023571, 0.98723000444482, 0.99109127151852, 0.99427354459409, 
0.996765916137424, 0.998557478614416, 0.99963732449096, 1, 0.999852748289562, 0.999408266274722, 
0.99867473460606, 0.99764942640005, 0.996337795423744, 0.994739841677142, 0.992861018927298, 
0.99070405405774, 0.988266220184937, 0.985555697959473, 0.982572487381347, 0.979319315334084,
0.975798908701213, 0.972013994366259, 0.967964572329222, 0.963658823240683, 0.959094020217114,
0.954272890142043, 0.949200886782523, 0.943878010138553, 0.93830698709366, 0.932493271414898, 
0.926436863102266, 0.920140489039292, 0.913604149225974, 0.906836024312893, 0.899833387416523,
0.892601692303917, 0.885140938975074, 0.877456581197047, 0.869548618969838, 0.861419779176972, 
0.853075515585503, 0.84451582819543, 0.835743443890281, 0.826758362670055, 0.817568765185333,
0.808171924552586, 0.79857329453887, 0.788772875144184, 0.77877612013558, 0.768585756396587,
0.758201783927203, 0.747626929610955, 0.736863920331371, 0.725915482971976, 0.714784344416297,
0.703473231547861, 0.691984871250194, 0.680319263523297, 0.668481862134223, 0.656475393966497,
0.644299859020122, 0.631957984178622, 0.619455223209051, 0.60679157611141, 0.593969769769224, 
0.58099253106602, 0.5678598600018, 0.554579937227141, 0.541150035858518, 0.527575609662984, 
0.513859385524066, 0.500001363441763, 0.486497836217921, 0.473117018752778, 0.459853457279279,
0.446715332448005, 0.433705371142482, 0.420823573362711, 0.408075392875744, 0.395463556565108,
0.382993518197857, 0.370662550890464, 0.358478835293508, 0.34644237140699, 0.334555886114436, 
0.322827560066427, 0.311254666379435, 0.299842658820514, 0.288591537389663, 0.277509482737464,
0.266596494863915, 0.255855300652543, 0.245288626986875, 0.234904654517491, 0.224697929477338,
0.214679359400522, 0.204846217403516, 0.195203957253374, 0.185755305833622, 0.176505716911314,
0.167452463602922, 0.158603726559027, 0.149959505779629, 0.141525255031782, 0.133303701199011,
0.125294844281316, 0.117506864929278, 0.109937036259370, 0.102593538922172, 0.0954736460341569,
0.0885855382459048, 0.0819319424409426, 0.0755128586192698, 0.0693337405479398, 0.0633973151104797,
0.0577035823068888, 0.0522607227877474, 0.0470687365530556, 0.042127623602813, 0.0374482914711263,
0.0330252863909423, 0.0288695158963675, 0.0249782531038750, 0.0213542248969919, 0.0180056119262977,
0.0149324141917924, 0.0121346316934765, 0.00962044508192914, 0.00738985435715067, 0.00544831328619442,
0.00379582186906052, 0.00243783387280174, 0.00137707618094482, 0.000613548793490332,
0.000155432361017495, 0, 0.000136344176331169, 0.000542649821798078, 0.00122164381992747,
0.00217877993777239, 0.00341678505885973, 0.0049356591831892, 0.00674085607781414, 0.00883510262626144,
0.0112211257120574, 0.0139016522187286, 0.0168794090298020, 0.0201571230288041, 0.0237347942157347, 
0.0276206032411738, 0.0318172769886478, 0.0363220885746306, 0.0411432186497016, 0.0462806672138611,
0.0517398880341624, 0.057520881110605, 0.0636291002102425, 0.0700645453330751, 0.0768326702461555, 
0.0839334749494843, 0.0913724132101146, 0.0991522119115726, 0.107275597937385, 0.115745298171079, 
0.124564039496181, 0.133734548796217, 0.143259552954715, 0.153141778855200, 0.0274869859483690, 
0.0274869859483690, 0.54319247161996))


myT<-data.frame(V1 = c(0.855, 0.855, 0.124670337305034, 0.124670337305034, 0.401045089797478, 0.401045089797478, 
				0.588629763434884, 0.588629763434884),
		V2 =  c(0.828531745918341, 1, 1, 0.828531745918341, 0.828531745918341, 0, 0, 0.828531745918341))

Logo <-function(matrix, maxScore=2, letterType = 2,title = "", start=1, weights = c(0.25,0.25,0.25,0.25), 
		ytitle = "", xtitle = "", myAxis = TRUE, starting, ending, sites =0, bases, 
		myColors= c("green", "blue", "orange", "red"), myLetters = c("A", "C", "G","T")){
	#sites is the total number of sequences and bases is the number od bases (DNA =4)
	if ( missing(bases)){
		bases = nrow(matrix)
	}
	if (missing(starting)){
		starting<-1
		ending<-length(matrix[1,])
	}
	else{
	}
	plot(c(starting,ending+1), c(0,maxScore),  type = "n",axes = FALSE, main = c(title), ylab = ytitle, xlab =xtitle)
	if(is.integer(length(matrix[1,])/2) == TRUE){
		if (length(matrix[1,])<20){
			breaks <-2
		}
		else{
			breaks <-length(matrix[1,])/10 
		}
	}
	else{
		if (length(matrix[1,])<20){
			breaks <-2
		}
		else{
			breaks <-(length(matrix[1,])-1)/10 
		}
	}
	if (myAxis == TRUE){
		if(start<0){
			axis(1, seq(1, length(matrix[1,])+1, breaks),labels = FALSE ) 
			end<-length(matrix[1,])+start-1
			axis(1, seq(1+0.5, length(matrix[1,])+0.5, breaks),seq(start,end,breaks), tick = FALSE)
			par(mgp=c(c(3, 1, 0)))
			axis(2, seq(0, maxScore,maxScore/5)) #seq(0, 4, 0.8),
		}
		else{
			axis(1, seq(1, length(matrix[1,])+1, breaks),labels = FALSE ) 
			end<-length(matrix[1,])+start-1
			axis(1, seq(1+0.5, length(matrix[1,])+1.5, breaks),seq(start,end+1,breaks), tick = FALSE)
			par(mgp=c(c(3, 1, 0)))
			axis(2, seq(0, maxScore,maxScore/5)) #seq(0, 4, 0.8),
		}	
	}
	else{
	}
	if (letterType ==1){
	}
	else{
		myA<-myA
		myC<-myC
		myG<-myG
		myT<-myT
	}
	if(missing(weights)){
	}
	else{
		weights <-(0.25-weights)+c(0.25, 0.25,0.25, 0.25)#to make the weight the rigth number otherwise it would be inverted
		matrix<-matrix*(weights/0.25)	
	}
	#freqMatrix<-matrix*sites #just if using variance
	for ( i in seq(starting,ending,1)){
		AEHnb = -sum(matrix[,i]*(log2(matrix[,i])), na.rm = TRUE) - ((bases-1)/(2*log(2)*sites)) 
		#in here uses the natural log instead of the log of two
		mySum<-sum(matrix[,i]*(log2(matrix[,i])), na.rm = TRUE)#just to compensate for those that are equal to 0}
		myEntropy<-log2(bases)+mySum
		inner<-0.00
		myTemp <- sort(matrix[,i])#to start ploting from the less frequent one
		myTempOr <-order(matrix[,i])#to get the right position of the letters 
		for ( j in 1:length(myTemp)){
			if (myTemp[j] == 0){
				next;
			}
			if(myLetters[myTempOr[j]] =="A"){
				multy<-myTemp[j]*myEntropy
				polygon(myA[[1]]+i,(myA[[2]]*multy)+inner,  lty = 0,col  = myColors[[1]] )
				inner<-inner+multy
			}
			else if(myLetters[myTempOr[j]] =="C"){
				multy<-myTemp[j]*myEntropy
				polygon(myC[[1]]+i,(myC[[2]]*multy)+inner, lty = 0, col  = myColors[[2]] )
				inner<-inner+multy
			}
			else if(myLetters[myTempOr[j]] =="G"){
				multy<-myTemp[j]*myEntropy
				polygon(myG[[1]]+i,(myG[[2]]*multy)+inner, lty = 0, col  = myColors[[3]] )
				inner<-inner+multy
			}
			else if(myLetters[myTempOr[j]] =="T"){
				multy<-myTemp[j]*myEntropy
				polygon(myT[[1]]+i,(myT[[2]]*multy)+inner, lty = 0, col  = myColors[[4]] )
				inner<-inner+multy
			}
			else{
				multy<-myTemp[j]*myEntropy
				print (myLetters[myTempOr[j]])
				inner<-inner+multy
			}
		}
	}	
}

RandomString <- function(n=1, lenght=10)
{
	randomString <- NULL
	aBucket<-c(0:9, letters, LETTERS)
	for (i in 1:n)
	{
		randomString[i]<-paste(sample(aBucket,lenght, replace=TRUE),collapse="")
	}
	return(randomString)
}
