#function to get frequency of hets at each population_cycle
get.hetfreq <- function(dataset){
	hets <- is.het(dataset$calls[1:length(dataset$calls)])
	hetfreq <- (sum(hets, na.rm=T))/(sum(hets, na.rm=T) + sum(!hets, na.rm=T))
	return(hetfreq)
}


#get HWE significance of each marker in a set
get.HWE <- function(marker){
	if (length(levels(as.factor(marker))) > 1) {
	marker <- genotype(marker, sep="")
	teststat <- HWE.exact(marker)$p.value
	}
	else{teststat <- NA}
	return(teststat)
}

#get the HWE probabilities for a dataset
dataset.HWE <- function(dataset){
	pvals <- apply(dataset$calls, 1, get.HWE)
	return(pvals)	
}


#count number of switches in phased data
countphase <- function(prephased, phased) {
	prephased[prephased=="?"] <- NA
	hetcode <- apply(prephased, 2, function(x){
			hets <- (x[1] != x[2])
			hets_na <- hets==TRUE & !is.na(hets)
			return(hets_na)
			})	
	prephasedhets <- prephased[,hetcode]
	phasedhets <- phased[,hetcode]

	phasecount <- c()
	for (i in 1:ncol(phasedhets)){
		phasecount[i] <- which(phasedhets[,i]==prephasedhets[1,i])
	}
	#print(paste("# of sites:", ncol(phasedhets)))
	possible <- ncol(phasedhets)-1
	switches <- sum(abs(diff(phasecount)))
	#print(paste("# of switches", switches ))
	switchrate <- switches/(ncol(phasedhets)-1)
	#ratename<-"rate="
	#switchrateout <- sprintf("%s\t%.7f\n", ratename,switchrate) 
	#cat(switchrateout)
	results <- list(possible, switches)
	return(results)	

}		

calcH <- function(genotypes){
		alleles <- unlist(strsplit(genotypes, split=""))
		allelecounts <- table(alleles)
		total <- sum(allelecounts)
		if(length(allelecounts)==2){
		freq1 <- allelecounts[1]/total
		freq2 <- allelecounts[2]/total
		H <- 1-(freq1^2 + freq2^2)
		}
		else if (length(allelecounts)==1){
			H<-0
			}
		else if(length(allelecounts)==0){
			H <- NA}	
		return(as.numeric(H))
	}
	
averageH <- function(calls){
	het <- apply(calls, 1, calcH)
	finalH <- mean(het)
	return(finalH)
	}


#fishers exact test on a single marker
fishers <- function(pop1, pop2, alleles){
	pop1alleles <- unlist(strsplit(pop1, split=""))
	pop2alleles <- unlist(strsplit(pop2, split=""))
	pop1counts <- c(sum(pop1alleles==alleles[1], na.rm=T), sum(pop1alleles==alleles[2], na.rm=T))
	pop2counts <- c(sum(pop2alleles==alleles[1], na.rm=T), sum(pop2alleles==alleles[2], na.rm=T))
	totest <- rbind(pop1counts, pop2counts)
	score <- fisher.test(totest)$p.value
	return(score)
	}	

apply.fishers <- function(dataset){
	alleles <- dataset$alleles
	calls <- dataset$calls
	pop1 <- which(dataset$group=="BSSS")
	pop2 <- which(dataset$group=="BSCB")
	scores <- c()
	for (i in 1:nrow(calls)){
		scores[i] <- fishers(calls[i,pop1], calls[i,pop2], alleles[i,])
		
		
		}
	return(scores)
	}

genetic.coverage <- function(dataset, chromosome, window=5, step=5) {
	
	tokeep <- dataset$LG==chromosome
	dataset <- row.prune(dataset,tokeep)
	measures <- dataset$calls
	positions <- dataset$F2map
	checkpts <- seq((min(positions) + window/2), (max(positions)-window/2), by= step)
	

	number <- sapply(checkpts, function(x){
		return(length(positions[(positions > (x-window/2)) & (positions < x+ (window/2)) & (positions != Inf) & (positions != -Inf)]))})
	
	gendata <- sapply(checkpts, function(x) {
		index <- which((positions > (x-window/2)) & (positions < x+ (window/2)) & (positions != Inf) & (positions != -Inf))
		start <- min(index)
		stop <- max(index)
		genstart <- dataset$F2map[start]
		genstop <- dataset$F2map[stop]
		physstart <- dataset$position[start]
		physstop <- dataset$position[stop]
		cmMb <- (genstop-genstart)/(physstop/1000000-physstart/1000000)
		
		return(cmMb)
		})
		
		finaldata <- data.frame(markers=number, bin=1:length(number), chrom=rep(chromosome, length(number)), rate=gendata)
		return(finaldata)
}


simH <- function(simpop){
	values <- apply(simpop, 2, function(line){
		allelecounts <- table(line)
		total <- sum(allelecounts)
		if(length(allelecounts)==2){
			freq1 <- allelecounts[1]/total
			freq2 <- allelecounts[2]/total
			H <- 1-(freq1^2 + freq2^2)
			}
		else if (length(allelecounts)==1){
			H<-0
		}
		return(as.numeric(H))
	})
}
 
het_per_line <- function(dataset){
		hetcounts <- apply(dataset$calls,2,function(x){
			hets <- sapply(x, is.het)
			return(sum(hets, na.rm=T))
			})
			return(hetcounts)
			}

num.poly <- function(calls){
	alleles <- apply(calls, 1, function(x){
		counts <- rle(sort(as.vector(x)))
		counts <- length(counts$values)
		return(counts)
		})
	alleles <- alleles!=1
	poly <- sum(alleles)
	return(poly)
	}

sim.num.poly <- function(calls){
	alleles <- apply(calls, 2, function(x){
		counts <- rle(sort(as.vector(x)))
		counts <- length(counts$values)
		return(counts)
		})
	alleles <- alleles!=1
	poly <- sum(alleles)
	return(poly)
	}


HBK.permute <- function(dataset){
	measures <- dataset$calls
	firstnum = as.numeric(table(dataset$group)[dataset$group[1]])
	secondnum = as.numeric(table(dataset$group)[dataset$group[length(dataset$group)]])
	write.fasta(measures, "tempfile.fa")
	command <- paste(" HBKpermute -i tempfile.fa -c 2 ", 2*firstnum, " ", 2*secondnum, " -s -n 0 ", sep="")
	system(command)
	rmcommand <- "rm tempfile.fa"
	system(rmcommand)
} 

is.poly <- function(dataset){
  alleles <- apply(dataset$calls, 1, function(x){
  	counts <- rle(sort(as.vector(x)))
		counts <- length(counts$values)
		return(counts)
		})
	alleles <- alleles!=1
	return(alleles)  
}

#define X marker windows
define.even.windows <- function(sequence, window){

checkpts <- seq(min(sequence), max(sequence), by=window)
	ends <- checkpts+window-1
	if(max(ends)>max(sequence)){
		ends[which.max(ends)]=max(sequence)
	}
	return(cbind(checkpts,ends))
}


define.genetic.windows <- function(dataset, chromosome, window, step){
	tokeep <- dataset$LG==chromosome
dataset <- row.prune(dataset,tokeep)
positions <- dataset$F2map
checkpts <- seq(min(positions), max(positions), by=step)
bounds=c()
for(i in checkpts){
	currentwindow=which(positions>=i & positions < i+window)
	currentbound=c(min(currentwindow), max(currentwindow))
	bounds= rbind(bounds,currentbound)
	}	
	bounds=bounds[which(bounds[,1] != Inf & bounds[,1] != -Inf),]
	bounds=bounds[which(bounds[,2] != Inf & bounds[,2] != -Inf),]
	return(bounds)
}
 
define.overlapping.windows <- function(sequence, window, step){
 
        checkpts <- seq(min(sequence), max(sequence), by=step)
         ends <- checkpts+window-1
         ending=which(ends>max(sequence))
         ends <- ends[1:(min(ending)-1)]
         checkpts <- checkpts[1:(min(ending)-1)]           
         
        return(cbind(checkpts,ends))
 }
