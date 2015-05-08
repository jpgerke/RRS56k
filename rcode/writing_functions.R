write.flopfile <- function(dataset, outfile){
	#transpose the matrix and replace NA's
	calls <- dataset$calls
	calls <- t(calls)
	missing <- is.na(calls)
	calls[which(missing==TRUE)] <- "--"
	#add the sequence position
	calls <- rbind(as.character(dataset$position), calls)
	#add chromosome
	calls <- rbind(as.character(dataset$LG), calls)
	#add line names
	first_column <- c("Chromosome", "Position", dimnames(dataset$calls)[[2]])
	calls <- cbind(first_column, calls)
	write.table(calls, file=outfile, quote=F, row.names=F, col.names=F, sep=" ")
}

#write a dataset to PLINK outfiles	
write.plink <- function(dataset, pedfile, mapfile, delimiter, assn="FALSE"){
	pheno <- 1
	linetag <- "0 0 0"
	#write the pedfile
	calls <- dataset$calls
	calls[which(is.na(calls))] <- "00"
	if(delimiter=="group"){
		popname <- as.character(dataset$group)
		}
	else if(delimiter=="status"){
		popname <- as.character(dataset$status)
		}
	else if (delimiter=="cycle"){
		popname <- paste(as.character(dataset$cycle), as.character(dataset$status), sep="_")
		}		
	if(assn=="TRUE"){
		pheno <- as.numeric(as.factor(dataset$group))
		linetag <- ""
		}	
	
	lines <- apply(rbind(popname, dimnames(calls)[[2]], pheno, calls), 2, function(line) {
		var <- unlist(strsplit(line[-(1:2)], split=""))
		paste(c(line[1:2], linetag, var), collapse=" ")
		})
	cat(lines, file=pedfile, sep="\n")
			
	#write the mapfile
	mapdata <- cbind(as.character(dataset$LG), as.character(dataset$markernames), rep(0,length(dataset$LG)), as.character(dataset$position))
	write.table(mapdata, file=mapfile, col.names=F, row.names=F, quote=F, sep="\t")
	}

#write out a fasta file
write.fasta <- function(calls, outfile) {
	missing <- is.na(calls)
	calls[which(missing==TRUE)] <- "NN"
	
	lines <- apply(rbind(dimnames(calls)[[2]], calls), 2, function(line){
		allsequence <- unlist(strsplit(line[-1], split=""))
		allsequence <- matrix(allsequence, nrow=2, byrow=FALSE)
		newline1 <- c(paste(">",line[1], sep=""), paste(allsequence[1,], collapse=""))
		newline2 <- c(paste(">",line[1], sep=""), paste(allsequence[2,], collapse=""))
		bothlines <- cbind(newline1, newline2)
		return(bothlines)
		}) 
cat(lines, file=outfile, sep="\n")
	}





write.window <- function(dataset, prefix, window=1000, step=1000){
	measures <- dataset$calls
	positions <- 1:dim(measures)[[1]]	
	checkpts <- seq((window/2), length(positions)-window/2, by=step)
	sapply(checkpts, function(x){
		current <- measures[(positions >= ((x+1)-window/2)) & (positions <= (x+window/2)),]
		filename <- paste(prefix, x, sep="_")
		outfile <- file(filename, "w")
		write.fasta(current, outfile)
		close(outfile)
	return(filename)
		})
	}
	
	
	phaser <- function(dataset, chromosome, snpfile, hapfile, popfile) {
	calls <- dataset$calls[dataset$LG==chromosome,dataset$status=="Outbred"]
	calls[which(is.na(calls))] <- "BB"
	ind <- dim(calls)[[2]]
	loc <- dim(calls)[[1]]	
	cat(ind, loc, sep="\n", file=snpfile)
	pops <- as.numeric(dataset$group[dataset$lines%in%dimnames(calls)[[2]]])
	
	lines <- apply(rbind(dimnames(calls)[[2]], calls), 2, function(line){
		allsequence <- unlist(strsplit(line[-1], split=""))
		allsequence <- matrix(allsequence, nrow=2, byrow=FALSE)
		name <- paste("#", line[1], sep=" ")
		allsequence[which(allsequence=="B")] <- "?"
		newline1 <- paste(as.character(allsequence[1,]), collapse="")
		newline2 <-  paste(allsequence[2,], collapse="")
		alllines <- c(name, newline1, newline2)
		return(alllines)
		})
	cat(lines, file=snpfile, sep="\n", append=T)
	
	haps <- dataset$calls[dataset$LG==chromosome,dataset$status!="Outbred"]
	haps <- remove.hets(haps)	
	hapind <- dim(haps)[[2]]
	haploc <- dim(haps)[[1]]	
	haps[which(is.na(haps))] <- "BB"
	happops <- as.numeric(dataset$group[dataset$lines%in%dimnames(haps)[[2]]])
	cat(hapind, sep="\n", file=hapfile)
	haplines <- apply(rbind(dimnames(haps)[[2]], haps), 2, function(line){
		allsequence <- unlist(strsplit(line[-1], split=""))
		allsequence <- matrix(allsequence, nrow=2, byrow=FALSE)
		name <- paste("#", line[1], sep=" ")
		allsequence[which(allsequence=="B")] <- "?"
		newline1 <- paste(as.character(allsequence[1,]), collapse="")
		alllines <- c(name, newline1)
		return(alllines)
		})
	cat(haplines, file=hapfile, sep="\n", append=T)
	pops <- c(pops, happops)
	cat(pops, file=popfile, sep=" ")
	}

inbredplinker <- function(dataset, pedfile, mapfile) {
	 
	 calls <- dataset$calls
#	 calls[is.na(calls)] <- 0
	 markername <- dataset$markernames
	 physmap <- dataset$position
	 chrom <- dataset$LG
	 genmap <- round(dataset$F2map, digits=7)
	 #genmap[is.na(genmap)] <- 0
	 mapdata <- sprintf("%d\t%s\t%.7f\t%d",  chrom, markername, genmap, physmap)
  cat(mapdata,  file = mapfile,  sep = "\n")
  
  strains <- dimnames(calls)[[2]]
  lines <- apply(rbind(colnames(calls), calls), 2, function(line) 	{
  	var <- paste(rep(line[-1],  each = 2),  collapse = " ")
  	strain <- line[1]
  	printline <- sprintf("%s 1 0 0 0 1 %s", strain, var )
  	return(printline)
  	}
  	)
  	cat(lines, file= pedfile, sep="\n")
}


write.hapalleles <- function(dataset, chromosome, file){
	strains <- dimnames(dataset$calls)[[2]]
	nstrains <- length(strains)
	strains <- paste(dimnames(dataset$calls)[[2]], collapse=" ")
	nmarkers <- length(dataset$markernames[dataset$LG==chromosome])
	out <- file(file, "w")
  	on.exit(close(out))
  	calls <- dataset$calls[dataset$LG==chromosome,]
  	markernames <- dimnames(calls)[[1]]
  	f2map <- dataset$F2map[dataset$LG==chromosome]
  	cat(sprintf("markers %d strains %s\n", nmarkers, nstrains), file=out)
  	cat(sprintf("strain_names %s\n", strains), file=out)
  	for(i in 1:length(markernames)){
  		
  		cat(sprintf("marker %s, 3, %.2f\n", markernames[i], f2map[i]), file=out)
  		tabled <- create.problines(calls[i,])
  		write.table(tabled, file=out, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  		}
  	
  	
  	
	}


create.problines <- function(marker){
	alleles <- unique(marker)
	allele1number <- sum(marker==alleles[1])
	allele2number <- sum(marker==alleles[2])
	allele1prob <- 1/allele1number
	allele2prob <- 1/allele2number
	allele1vector <- marker==alleles[1]
	allele2vector <- marker==alleles[2]
	allele1vector[allele1vector==TRUE] <- allele1prob
	allele1vector[allele1vector==FALSE] <- 0.000
	allele1vector <- sprintf("%.3f", allele1vector)
	allele2vector[allele2vector==TRUE] <- allele2prob
	allele2vector[allele2vector==FALSE] <- 0.000
	allele2vector <- sprintf("%.3f", allele2vector)
	navector <- rep(1/length(marker), length(marker))
	navector <- sprintf("%.3f", navector)
	line1 <- c("allele", "ND", navector)
	line2 <- c("allele", alleles[1], allele1vector)
	line3 <- c("allele", alleles[2], allele2vector)
	alleleframe <- rbind(line1,line2,line3)
	return(alleleframe)
	}
	

write.happygenotypes <- function(dataset, chromosome, file){	calls <- dataset$calls[dataset$LG==chromosome,]
	strains <- dimnames(calls)[[2]]
  lines <- apply(rbind(colnames(calls), calls), 2, function(line) 	{
  	var <- paste(rep(line[-1],  each = 2),  collapse = " ")
  	strain <- line[1]
  	printline <- sprintf("%s 1 %s", strain, var )
  	return(printline)
  	}
  	)
  	cat(lines, file= file, sep="\n")
}

	
	
	
	
