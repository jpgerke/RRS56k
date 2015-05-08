
#function to get genotypes from a table of lines

get.genos <- function(linetable, sample_codes, genotypes){
	markernames <- genotypes[,1]
	justgenotypes <- genotypes[,-1]
	linenames <- linetable$Line
	calls <- matrix(ncol=length(linetable$SampleID), nrow=dim(genotypes)[[1]])
	for(i in 1:length(linetable$SampleID)){
		calls[,i] <- as.character(justgenotypes[,which(sample_codes == linetable$SampleID[i])])
	}
	dimnames(calls) <- list(markernames,linenames)
	return(calls)
}


get.scores <- function(linetable, sample_codes, genotypes){
	markernames <- genotypes[,1]
	justgenotypes <- genotypes[,-1]
	linenames <- linetable$Line
	calls <- matrix(ncol=length(linetable$SampleID), nrow=dim(genotypes)[[1]])
	for(i in 1:length(linetable$SampleID)){
		calls[,i] <- justgenotypes[,which(sample_codes == linetable$SampleID[i])]
	}
	dimnames(calls) <- list(markernames,linenames)
	return(calls)
}


#function to change hets to missing data
is.het <- function(call){
	return(substr(call,1,1) != substr(call,2,2))
	}

remove.hets <- function(calls){
	hets <- is.het(calls)
	calls[which(hets=="TRUE")] <- NA
	return(calls)
}

#cut out markers with a certain missing data threshold
na.cut.markers <- function(dataset, threshold) {
	nacounts <- apply(dataset$calls, 1, function(x){sum(is.na(x))})
	torun <- dataset
	keep <- nacounts <= threshold
	newset <- row.prune(dataset=dataset, keep=keep)
	return(newset)
	}

#remove a particular line from the dataset	
remove.line <- function(line, dataset){
	column <- dimnames(dataset$calls)[[2]]!=line
	newset <- column.prune(dataset=dataset, keep=column)
	return(newset)
	}	

#remove monomorphic calls	
remove.monos <- function(dataset){
	alleles <- apply(dataset$calls, 1, function(x){
		counts <- rle(sort(as.vector(x)))
		counts <- length(counts$values)
		return(counts)
		})
	alleles <- alleles!=1
	newset <- row.prune(dataset=dataset, keep=alleles)
	return(newset)
	}		
	
	get.alleles <- function(genotypes){
	alleles <- names(table(unlist(strsplit(genotypes, split=""))))
	return(alleles)
	}
	
add.line <- function(dataset, name, genotype, group, status, cycle){
	newname <- dataset
	newname$calls <- cbind(dataset$calls, genotype)
	dimnames(newname$calls)[[2]][length(dimnames(newname$calls)[[2]])] <- name
    newname$scores <- cbind(dataset$scores, rep(NA, length(genotype)))
    dimnames(newname$scores)[[2]] <- dimnames(newname$calls)[[2]]
	newname$group <- as.factor(c(as.character(dataset$group), group))
	newname$status <- as.factor(c(as.character(dataset$status), status)) 
	newname$cycle <-	as.factor(c(as.character(dataset$cycle), cycle))
	newname$lines <- as.factor(c(as.character(dataset$lines), name))
	return(newname)
}



rF2.to.rRIL <- function(rF2, generations=4){
	num <- 1-2*rF2
	denom <- 1+2*rF2
	factor <- (1-rF2)^generations
	term <- (num/denom)*factor
	rRIL <- 0.5*(1-term)
	return(rRIL)
	}
	
cmF2.to.cmRIL <- function(cmF2){
	rF2 <- 0.5*(1-(exp(-2*(cmF2/100))))
	rRIL <- rF2.to.rRIL(rF2)
	cmRIL <- -0.5*(log(1-2*rRIL))
	return(cmRIL*100)
	}	


remove.line <- function(line, dataset){
	column <- which(dimnames(dataset$calls)[[2]]!=line)
	newset <- column.prune(dataset, column)
	return(newset)
	}	

