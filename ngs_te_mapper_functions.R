#!/usr/bin/Rscript --vanilla --slave
# chmod u+x functions.R
# 1:00:56 PM
# TODO: Add comment 
# Author: raquel
###############################################################################



GetPSL<-function(file, fields="all", fieldNames = c("match", "mismatch", "repmatch", "NCount", 
				"QgapCount", "QgapBases", "TgapCount", "TgapBases", "strand", "Qname", "Qsize", 
				"Qstart", "Qend", "Tname", "Tsize", "Tstart", "Tend", "blockCount", "blockSizes", 
				"qStarts", "tStarts"), fieldTypes = list(integer(),integer(),integer(),integer(),
				integer(),integer(),integer(),integer(),character(), character(), integer(),integer(),
				integer(),character(),integer(),integer(),integer(),integer(),character(), character(),
				character()))
{
	internalVar<-new.env();
	if (fields[1] == "all"){
		afile<-scan(file, fieldTypes, sep = "\t", quote = "", fill = TRUE, skip = 5);
		for ( i in 1:length(fieldNames)){
			assign(fieldNames[i], afile[[i]], envir=internalVar);
		}
	}
	else{
		fields<-sort(fields);
		afile<-scan(file, fieldTypes, sep = "\t", quote = "", fill = TRUE,skip = 5);
		for ( i in fields){
			assign(fieldNames[i], afile[[i]], envir=internalVar);
		}
	}
	return(as.list(internalVar));	
}

#only compare the quality matchsize and blockcount
#take 3 since the minimumis 4 and add 1 since position 3 was taken
GetBest<-function(test1, test2, notComparable, toCompare, add=3, take = 3, aTest = sample(c(1,2),1))
{
	different<-which(test1[-notComparable] != test2[-notComparable])
	if(length(different) == 0)
	{
		#select a random one
		if(aTest == 2)
		{
			return(test2)
		}
		if(aTest == 1)
		{
			return(test1)
		}
		else
		{
			return("NA")
		}	
	}
	else
	{
		#the choices are difference in the matchsize, quality and blockcount
		#take only the maximum of the two since I satrt chosing by the highest of them
		different<-max(different) + add 
		comp<-different-take
		result<-eval(call(toCompare[comp], test1[different],test2[different]))
		if(result == TRUE)
		{
			return(test2)
		}
		else
		{
			return(test1)
		}
		
	}
}

Overlap<-function(vector1, vector2, gap = 5)
{
	#they have to be the same length so that there is only a type of overlap
	if(vector1[1] <= vector2[1] & vector1[2]<=vector2[2])
	{
		return(TRUE)
	}
	else if(vector1[1] >= vector2[1] & vector1[2] >=vector2[2])
	{
		return(TRUE)
	}
	else
	{
		if (vector1[1]< vector2[1])
		{
			if((vector2[1]-vector1[2])<=gap)
			{
				return(TRUE)
			}
			else{
				return(FALSE)
			}
		}
		else
		{
			if((vector1[1]-vector2[2]) <=gap)
			{
				return(TRUE)
			}
			else{
				return(FALSE)
			}
		}
	}
}

SelectFirstReads<-function(aPslFile, tolerated = 20, mid="middle", st="start", en="end")
{
	aPslFile<-data.frame(aPslFile)
	totalGap = aPslFile$QgapBases + aPslFile$TgapBases;
	totalMiss = aPslFile$mismatch + aPslFile$repmatch + aPslFile$NCount;
	quality = totalGap +  totalMiss;
	matchsize = aPslFile$Qend - aPslFile$Qstart;
	uniqueBlockStarts = split (aPslFile$tStarts, "," );	
	tolerance = matchsize/tolerated
	toTake<-c(which(quality > tolerance))
	if(length(toTake) == 0)
	{
	}
	else
	{			
		quality<-quality[- toTake ]
		matchsize<-matchsize[- toTake ]
		tolerance<-tolerance[-toTake]
		aPslFile<-aPslFile[- toTake, ]
	}
	toTake<-which(aPslFile$blockCount > tolerance)
	if(length(toTake) == 0)
	{
	}
	else
	{
		quality<-quality[- toTake ]
		matchsize<-matchsize[- toTake ]
		tolerance<-tolerance[-toTake]
		aPslFile<-aPslFile[- toTake, ]
	}
	location<-rep(mid, length(aPslFile$Qstart))
	location[which(aPslFile$Tstart == 0)]<-st
	location[which(aPslFile$Tend == aPslFile$Tsize )]<-en
	queryStarts<-paste(aPslFile$Qstart, aPslFile$Qend, sep ="-")
	difference<-aPslFile$Qsize - aPslFile$Qend + aPslFile$Qstart
	aPslFile$Qname<-as.character(aPslFile$Qname)
#	to keep
	id<-strsplit(paste(aPslFile$Tname,location,aPslFile$strand,aPslFile$blockCount,quality,matchsize,queryStarts,
					aPslFile$Qsize,difference, sep = ","), split = ",")
	tName<-1
	location<-2
	strand<-3
	blockCount<-4
	quality<-5
	match<-6
	qStart<-7
	qSize<-8
	diff<-9
	test<-id[order(aPslFile$Qname)]
	readTable<-table(aPslFile$Qname)
	notComparable<-c(tName, location, strand, qSize, qStart, diff)
	toCompare<-c(">", ">", "<")
	toAddMid<-7
	toKeep<-NULL
	position<-1
	aNameSet<-names(readTable)
	for ( i in 1:length(aNameSet))
	{
		if(readTable[[i]] == 1)
		{
			toKeep[i]<-test[position]
		}
		else
		{
			toKeep[i]<-test[position]
			for(j in 1:(readTable[[i]]-1))
			{
				if(toKeep[i]=="NA")
				{
					#if the previous one was NA only keep if the aligment is better
					toKeep[[i]]<-GetBest(test[[position+j-1]], test[[position+j]], notComparable, toCompare, aTest = "NA")
				}
				else if(identical(toKeep[[i]] ,test[[position+j]]))
				{
					#they are the same so it does not make a difference
				}
				else
				{
					if(toKeep[[i]][tName] == test[[position+j]][tName])
					#they are the same so we don't need to care about overlaps
					{
						if(toKeep[[i]][location] == test[[position+j]][location])
						#they are the same so we don't need to care about overlaps
						{
							toKeep[[i]]<-GetBest(toKeep[[i]], test[[position+j]], notComparable, toCompare)
						}
						else
						{
							midPos<-grep(mid, c(toKeep[[i]],test[[position+j]]))
							if( length(midPos) > 0)
							{
								if(c(toKeep[[i]],test[[position+j]])[midPos+toAddMid]  == 0)
								{
									#keep the middle one since the read is all in there
									if(midPos > diff)
									{
										toKeep[[i]]<-test[[position+j]]
									}
								}
								else
								{
									#keep either the start or end one if they are equal since they are more important
									if(midPos > diff)
									{
										aTest<-1
									}
									else
									{
										aTest<-2
									}
									toKeep[[i]]<-GetBest(toKeep[[i]], test[[position+j]], notComparable, toCompare, aTest = aTest)
								}
							}
							else
							{
								#it is a start and an end so selct a random one is they are equal
								toKeep[[i]]<-GetBest(toKeep[[i]], test[[position+j]], notComparable, toCompare)
							}
						}
					}
					else
					{
						#test overlap
						temp<-as.numeric(unlist(strsplit(toKeep[[i]][qStart], split = "-")))
						temp2<-as.numeric(unlist(strsplit(test[[position+j]][qStart], split = "-")))
						if(Overlap(temp, temp2) == TRUE)
						{
							toKeep[[i]]<-GetBest(toKeep[[i]], test[[position+j]], notComparable, toCompare, aTest = "NA")
						}
						else
						{
							toKeep[[i]]<-GetBest(toKeep[[i]], test[[position+j]], notComparable, toCompare)
						}
					}	
				}
			}
		}
		toKeep[[i]]<-paste(toKeep[[i]], collapse = ",")
		position<-position+readTable[[i]]
	}
	internalVar<-new.env();
	toTake<-which(toKeep == "NA")#those that are umbiguoes with the matching to a specific TE
	print (paste(round(length(toTake)/length(toKeep), 5)*100  , 
					"% of the reads were ambigous in the match to the TE dataset", sep = ""))
	aNameSet<-aNameSet[-toTake]
	toKeep<-toKeep[-toTake]
	internalVar$toKeep = toKeep
	internalVar$aNameSet= aNameSet
	return(as.list(internalVar))	
}

RightNewFasta<-function(selectedReads, fastaFile, outputFile, mid="middle", st="start", en="end", noMid=TRUE)
{
	aFile<-file(fastaFile, open = "r")
	if (noMid == TRUE)
	{
		toTake<-grep(mid, selectedReads$toKeep)
		selectedReads$aNameSet<-selectedReads$aNameSet[-toTake]
		selectedReads$toKeep<-selectedReads$toKeep[-toTake]
		aNameSet<-paste(">", selectedReads$aNameSet, sep = "")
		newID<-paste(aNameSet,selectedReads$toKeep, sep = ",")
	}
	else
	{
		aNameSet<-paste(">", selectedReads, sep = "")
		newID<-aNameSet
	}
	myOutput<-file(outputFile, "w")
	toWrite<-TRUE
	while (length(x<-readLines(aFile, 1))>0) 
	{
		if(substr(x, 1, 1) == ">")
		{	
			aMatch<-match(x,aNameSet)
			if(is.na(aMatch))
			{
				toWrite = FALSE
				next;
			}
			else
			{
				cat(newID[aMatch], file = myOutput)
				cat("\n", file = myOutput)
				toWrite =TRUE
			}
		}
		else
		{
			if (toWrite == TRUE)
			{
				cat(x, file = myOutput)
				cat("\n", file = myOutput)
			}
			else
			{
				next;
			}
		}
	}
	close(aFile)
	close(myOutput)
}


SelectSecondReads<-function(aPslFile, bedFile, tolerated = 20, sizeTolerance = 20, mid="middle", 
		st="start", en="end", strands=c("+", "-"), idSplit=",", withRep  = 1, pasteChar = ";")
{
	aPslFile<-data.frame(aPslFile)
	newId<-strsplit(as.character(aPslFile$Qname), split = idSplit)
	newId<-matrix(unlist(newId), byrow = TRUE, nrow = length(newId))
	#select only for those that the full length of the TE non mapped sequence
	readId<-1
	tName<-2
	location<-3
	strand<-4
	blockCount<-5
	qual<-6
	match<-7
	qStart<-8
	qSize<-9
	diff<-10
	#1
	totalGap = aPslFile$QgapBases + aPslFile$TgapBases;
	totalMiss = aPslFile$mismatch + aPslFile$repmatch + aPslFile$NCount;
	quality = totalGap +  totalMiss;
	matchsize = aPslFile$Qend - aPslFile$Qstart;
	tolerance = matchsize/tolerated
	toTake<-which((newId[,qSize] != aPslFile$Qsize) | (quality > tolerance))
	if(length(toTake) == 0)
	{
	}
	else
	{	
		newId<-newId[- toTake,]
		aPslFile<-aPslFile[- toTake, ]
		tolerance<-tolerance[- toTake]
		quality<-quality[- toTake]
	}
	#2
	temp<-unlist(strsplit(newId[,qStart], split = "-"))
	aStart<-temp[seq(1, 2*length(newId[,qStart]), 2)]
	anEnd<-temp[seq(2, 2*length(newId[,qStart]), 2)]
	toTake<-which((aPslFile$Qend != aStart) & (aPslFile$Qstart != anEnd))
	if(length(toTake) == 0)
	{
	}
	else
	{	
		newId<-newId[- toTake,]
		aPslFile<-aPslFile[- toTake, ]
		aStart<-aStart[- toTake]
		anEnd<-anEnd[- toTake]
		tolerance<-tolerance[- toTake]
		quality<-quality[- toTake]
	}
	#3
	sizeTolerated<-aPslFile$Qsize/sizeTolerance
	aStart<-as.numeric(aStart)
	anEnd<-as.numeric(anEnd)
	difference<-aPslFile$Qsize - anEnd + aStart + aPslFile$Qstart - aPslFile$Qend
	toTake<-which((difference-1) > sizeTolerated)
	if(length(toTake) == 0)
	{
	}
	else
	{	
		newId<-newId[- toTake,]
		aPslFile<-aPslFile[- toTake, ]
		aStart<-aStart[- toTake]
		anEnd<-anEnd[- toTake]
		tolerance<-tolerance[- toTake]
		difference<-difference[- toTake]
		quality<-quality[- toTake]
	}
	#4
	toTake<-which(difference <0)
	if(length(toTake) == 0)
	{
	}
	else
	{	
		print(aPslFile[toTake, ])
		newId<-newId[-toTake,]
		aPslFile<-aPslFile[- toTake, ]
		aStart[-toTake]
		anEnd[-toTake]
		tolerance<-tolerance[-toTake]
		quality<-quality[-toTake]
	}
	#5
	toTake<-which(aPslFile$blockCount >1 & aPslFile$blockCount >tolerance)
	if(length(toTake) == 0)
	{
	}
	else
	{	
		newId<-newId[-toTake,]
		aPslFile<-aPslFile[- toTake, ]
		aStart<-aStart[- toTake]
		anEnd<-anEnd[- toTake]
		quality<-quality[-toTake]
	}
	#6
	toTake<-which(aPslFile$Qsize  != anEnd & aStart !=0 & aPslFile$Qend != aPslFile$Qsize & aPslFile$Qstart !=0 )
	if(length(toTake) == 0)
	{
	}
	else
	{	
		newId<-newId[-toTake,]
		aPslFile<-aPslFile[- toTake, ]
		aStart<-aStart[- toTake]
		anEnd<-anEnd[- toTake]
		quality<-quality[-toTake]
	}
	readIds<-table(newId[,readId])
	anOrder<-order(newId[,readId])
	predictedStrand<-rep(NA, length(aPslFile$strand))
	predictedStrand[which(aPslFile$strand == newId[,strand])]<-strands[1]
	predictedStrand[which(aPslFile$strand != newId[,strand])]<-strands[2]
	finalId<-paste(aPslFile$Tname, aPslFile$Tstart,  aPslFile$Tend, predictedStrand, newId[,tName],newId[,location], 
			sep = pasteChar)
	finalId<-finalId[anOrder]
	quality<-quality[anOrder]
	readName<-names(readIds)
	toKeep<-NULL
	position<-1
	for ( i in 1:length(readName))
	{
		if(readIds[[i]] == 1)
		{
			#print(finalId[position])
			toKeep[i]<-finalId[position]
		}
		else
		{
			tempo<-min(quality[seq(position,(position + readIds[[i]] -1), 1)])
			if(length(tempo) >1)
			{
				if((tempo[1] == 0) & length(tempo <= withRep))
				{
					tempo<-tempo-1 +position
					toKeep[i]<-paste(finalId[tempo], sep = idSplit)
				}
				else{
					toKeep[i]<-"NA"
				}
			}
			else
			{
				toKeep[i]<-finalId[position+tempo-1]
			}
		}
		position<-position+readIds[[i]]
	}
	readNameToKeep<-readName[- toTake]
	toTake<-which(is.na(toKeep) == TRUE)
	readNameToKeep<-readNameToKeep
	toKeep<-unlist(strsplit(toKeep, split = idSplit))
	#those that are umbiguoes with the matching to a specific TE
	print (paste(round(length(toTake)/length(toKeep), 5)*100  , 
					"% of the ambigous matching to the genome", sep = ""))
	if(length(toTake)>0)
	{
		toKeep<-toKeep[- toTake]
	}
	if(missing(bedFile)){
		
		internalVar<-new.env();
		internalVar$toKeep<-toKeep
		internalVar$readId<-readNameToKeep
		return(as.list(internalVar));	
	}
	else
	{
		myOutput<-file(bedFile, "w")
		aMat<-matrix(unlist(strsplit(toKeep, split = pasteChar)), nrow = length(toKeep), byrow = TRUE)
		chrom<-1
		start<-2
		end<-3
		bedID<-paste(readNameToKeep, toKeep, sep = ";")
		cat(paste(aMat[,1],aMat[,2], aMat[,3],bedID,sep = "\t"), sep = "\n", file = myOutput)
		close(myOutput)
		internalVar<-new.env();
		internalVar$toKeep<-toKeep
		internalVar$readId<-readNameToKeep
		return(as.list(internalVar));	
	}
	
}

PredictStrand<-function(nstarts,nends,strandStart,strandEnd, tePosStart,tePosEnd, decision = c("start", "end"))
{	
	aStartStrandTable<-table(strandStart)
	anEndtrandTable<-table(strandEnd)
	if(length(aStartStrandTable) == 1)
	{
		stStrand<-names(aStartStrandTable)
	}
	else
	{
		stStrand<-names(which(aStartStrandTable == max(aStartStrandTable)))
		if(length(stStrand) == 2)
		{
			aPosTable<-table(tePosStart)
			temp<-names(which(aPosTable == max(aPosTable)))
			if(length(temp) == 1)
			{
				if(temp == decision[2])
				{
					stStrand<-"+"
				}
				else
				{
					stStrand<-"-"
				}
			}
			else
			{
				stStrand<-NA
			}
		}
	}
	if(length(anEndtrandTable) == 1)
	{
		enStrand<-names(anEndtrandTable)
	}
	else
	{
		enStrand<-names(which(anEndtrandTable == max(anEndtrandTable)))
		if(length(enStrand) == 2)
		{
			aPosTable<-table(tePosEnd)
			temp<-names(which(aPosTable == max(aPosTable)))
			if(length(temp) == 1)
			{
				if(temp == decision[2])
				{
					enStrand<-"+"
				}
				else
				{
					enStrand<-"-"
				}
			}
			else
			{
				enStrand<-NA
			}
		}
	}
	bothStrands<-c(enStrand,stStrand)
	bothTable<-table(bothStrands)
	if(length(bothTable) == 1)
	{
		return(names(bothTable))
	}
	else
	{
		if(nstarts > nends)
		{
			return(stStrand)
		}
		else if (nstarts < nends)
		{
			return(enStrand)
		}
		else
		{
			return(NA)
		}
	}
}

GetDistMatrix<-function(positions, tsd)
{
	test<-as.matrix(dist(positions,  diag = TRUE, upper = TRUE))
	test[which(col(test) == row(test))]<-NA
	test[which(col(test) > row(test))]<-NA
	toKeep<-which(test<=tsd)
	aSize<-length(test[1,])
	theColumns<-as.integer(toKeep/aSize)
	theRows<-round((toKeep/aSize-theColumns)*aSize)
	theColumns<-theColumns+1
	select<-which(theRows == 0)
	theRows[select]<-aSize
	theColumns[select]<-theColumns[select]-1
	internalVar<-new.env();
	internalVar$theRows = theRows
	internalVar$theColumns= theColumns
	internalVar$realEnds = positions[theRows]
	internalVar$realStarts= positions[theColumns]
	return(as.list(internalVar))	
}

FinalProcessing<-function(secondReads, sample, tsd = 20, starts = "begin", ends = "stop", 
		splitChar = ";", lastSelection = "min" , distBetStart = tsd/2 , toPaste = ";")
{
	aList<-strsplit(secondReads, split = splitChar)
	chrom<-1
	start<-2
	end<-3
	firstStrand<-4
	firstTeName<-5
	location<-6
	aMatrix<-matrix(unlist(aList), byrow = TRUE, nrow = length(aList))
	aTeSet<-names(table(aMatrix[,firstTeName]))
	readPairs<-NULL
	m<-1
	for ( i in 1:length(aTeSet))
	{
		locations<-which(aMatrix[,firstTeName] == aTeSet[i])
		tempo<-matrix(data = aMatrix[locations,], ncol = location)
		newMatrix<-matrix(c(rep (tempo[,chrom],2,),tempo[,start],tempo[,end],rep(tempo[,firstStrand],2), 
						rep(tempo[,firstTeName],2),rep(tempo[, location], 2),rep(starts, length(tempo[,1])), 
						rep(ends, length(tempo[,1]))), byrow = FALSE, nrow = length(tempo[,1])*2) 
		genomeLocation<-2
		strand<-3
		teName<-4
		tePlace<-5#start or an end
		readLocation<-6 #begin or stops 
		newMatrix[,genomeLocation]<-as.numeric(newMatrix[,genomeLocation])
		anOrder<-order(newMatrix[,chrom], newMatrix[,genomeLocation])
		newMatrix<-newMatrix[anOrder,]
		aChNameSet<-names(table(newMatrix[,chrom]))
		for ( j in 1:length(aChNameSet))
		{
			tempo<-newMatrix[which(newMatrix[,chrom] == aChNameSet[j]),]
			aTable<-table(tempo[,genomeLocation])#get the genome locations
			genomeLocations<-as.numeric(names(aTable))
			firstDist<-GetDistMatrix(genomeLocations, tsd)
			tempEnd<-matrix(data = tempo[which(is.na(match(tempo[,genomeLocation],firstDist$realEnds)) == FALSE),], 
					ncol = readLocation)
			tempEnd<-matrix(data = tempEnd[which(tempEnd[,readLocation] == ends),],	ncol = readLocation)
			tempStart<-matrix(data = tempo[which(is.na(match(tempo[,genomeLocation],firstDist$realStarts)) == FALSE),],
					ncol = readLocation)
			tempStart<-matrix(data = tempStart[which(tempStart[,readLocation] == starts),],	ncol = readLocation)
			tablestarts<-table(tempStart[,genomeLocation])
			lastStarts<-names(tablestarts)	
			#first select the best starts the most abundant ones and those that agree 
			#more among each other in terms of strand
			if(length(lastStarts) == 0)
			{
				next;
			}
			startDist<-GetDistMatrix(lastStarts, distBetStart)
			#there is no start that has a small distance between with it
			#select the good start and end pairs
			toTake<-NULL
			if(length(startDist$realStarts)  == 0)
			{
				aFlag<-rep("0", length(tempStart[,genomeLocation]))
			}
			else
			{
				aFlag<-rep("0", length(tempStart[,genomeLocation]))
				for ( s in 1:length(startDist$realStarts))
				{
					startSet<-which(tempStart[,genomeLocation] == startDist$realStarts[s])
					startSet2<-which(tempStart[,genomeLocation] == startDist$realEnds[s])
					if(length(startSet) == length(startSet2))
					{	
						aTable1<-table(tempStart[startSet,strand])
						aTable2<-table(tempStart[startSet2,strand])
						if(length(names(aTable1)) > length(names(aTable2)))
						{
							#flag startSet2 has close (1)
							aFlag[startSet2]<-1
							toTake<-c(toTake, startSet)
						}
						else if(length(names(aTable1)) < length(names(aTable2)))
						{
							aFlag[startSet]<-1
							toTake<-c(toTake, startSet2)
						}
						else
						{
							#keep both
						}
					}
					else if(length(startSet) > length(startSet2))
					{
						aFlag[startSet]<-1
						toTake<-c(toTake, startSet2)
					}
					else
					{
						aFlag[startSet2]<-1
						toTake<-c(toTake, startSet)
					}
				}
			}			
			#get the best pair of reads for the starts
			#and keep the best pair or the pair with the most number of reads
			if(length(toTake)>0)
			{
				tempStart<-tempStart[- as.numeric(names(table(toTake))),]
				tablestarts<-table(tempStart[,genomeLocation])
				lastStarts<-names(tablestarts)			
				aFlag<-aFlag[- as.numeric(names(table(toTake)))]
			}
			position<-1
			if(length(lastStarts) == 0)
			{
				next;
			}
			for(s in 1:length(lastStarts))
			{
				distances<-as.numeric(tempEnd[,genomeLocation])-as.numeric(lastStarts[s])
				selected<-which(distances <= tsd & distances>0)
				if(length(selected) == 0)
				{
				}
				else
				{
					#print(paste("hello", s, i,j ,sep = " " ))
					nStarts<-tablestarts[[s]]
					startSelect<-seq(position, (position+tablestarts[[s]]-1),1)
					strandStart<-tempStart[startSelect,strand]
					tableEnds<-table(tempEnd[selected,genomeLocation])
					tePosStart<-tempStart[startSelect,tePlace]
					thisFlag<-names(table(aFlag[startSelect]))
					if(length(names(tableEnds)) == 1)
					{						
						anEnd<-names(tableEnds)[1]
						distance<-distances[selected[1]] #there is only one
						nEnds<-tableEnds[1]
						strandEnd<-tempEnd[selected,strand]
						tePosEnd<-tempEnd[selected,tePlace]
						aStrand<-PredictStrand(nStarts,nEnds,strandStart,strandEnd, tePosStart,tePosEnd)
						nReads<-nStarts+nEnds
					}
					else
					{
						anEnd<-names(tableEnds)[1]
						distance<-distances[1]
						nEnds<-tableEnds[1]
						selected<-which(tempEnd[,genomeLocation] == anEnd)
						strandEnd<-tempEnd[selected,strand]
						tePosEnd<-tempEnd[selected,tePlace]
						aStrand<-PredictStrand(nStarts,nEnds,strandStart,strandEnd, tePosStart,tePosEnd)
						for(t in 2:length(names(tableEnds)))
						{
							if(nEnds > tableEnds[t])
							{
								next
							}
							else if (nEnds < tableEnds[t])
							{
								anEnd<-names(tableEnds)[t]
								distance<-distances[t]
								nEnds<-tableEnds[t]
								selected<-which(tempEnd[,genomeLocation] == anEnd)
								strandEnd<-tempEnd[selected,strand]
								tePosEnd<-tempEnd[selected,tePlace]
								aStrand<-PredictStrand(nStarts,nEnds,strandStart,strandEnd, tePosStart,tePosEnd)	
							}
							else
							{
								tempEnd<-names(tableEnds)[t]
								tempdistance<-distances[t]
								tempEnds<-tableEnds[t]
								selected<-which(tempEnd[,genomeLocation] == anEnd)
								strandEnd<-tempEnd[selected,strand]
								tePosEnd<-tempEnd[selected,tePlace]
								tempStrand<-PredictStrand(nStarts,nEnds,strandStart,strandEnd, tePosStart,tePosEnd)		
								test<-which(is.na(c(tempStrand,aStrand)) == TRUE)
								if ( length(test) == 1)
								{
									if(test == 2)
									{
										anEnd<-tempEnd
										distance<-tempdistance
										nEnds<-tempEnds
										aStrand<-tempStrand
									}
								}
								else
								{
									if(lastSelection == "min" )
									{
										
									}
									else if(lastSelection == "max")
									{
										anEnd<-tempEnd
										distance<-tempdistance
										nEnds<-tempEnds
										aStrand<-tempStrand
									}
								}
							}
						}
					}
					readPairs[m]<-paste(aChNameSet[j],lastStarts[s], anEnd, distance, aStrand, aTeSet[i], 
							sample, nReads, thisFlag, sep = toPaste)
					m<-m+1
				}
				position<-position+tablestarts[[s]]
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

GetSequences<-function (fastaFile, chrom, start,winLeft, winRigth, strand,  complement = 1, 
		posPairs  = c("A","C","G","T"), negPairs = c("T","G", "C", "A"))
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

myA<-read.table("LogoA")
myC<-read.table("LogoC")
myG<-read.table("LogoG")
myT<-read.table("LogoT")

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
		AEHnb = -sum(matrix[,i]*(log2(matrix[,i])), na.rm = TRUE) - ((bases-1)/(2*log(2)*sites)) #in here uses the natural log instead of the log of two
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

