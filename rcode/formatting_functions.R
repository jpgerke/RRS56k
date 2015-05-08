
####get rid of het calls, for phasing calibration using derived inbreds
is.het <- function(call){
	return(substr(call,1,1) != substr(call,2,2))
	}

remove.hets <- function(calls){
	hets <- is.het(calls)
	calls[which(hets=="TRUE")] <- NA
	return(calls)
}






#split dataset by group	
split.lines <- function(dataset, groupname){
	keep <- which(dataset$group==groupname)
	newname <- dataset
	newname$calls <- dataset$calls[,keep]
    newname$scores <- dataset$scores[,keep]
   	newname$status <- dataset$status[keep]
   	newname$group <- dataset$group[keep]
   	newname$lines <- dataset$lines[keep]
	newname$cycle <- dataset$cycle[keep]
	return(newname)
	}	

#split based on derived/founder/outbred
split.status <- function(dataset, statusname){
	keep <- which(dataset$status==as.character(statusname))
	newname <- dataset
	newname$calls <- dataset$calls[,keep]
    newname$scores <- dataset$scores[,keep]
   	newname$status <- dataset$status[keep]
   	newname$group <- dataset$group[keep]
   	newname$lines <- dataset$lines[keep]
	newname$cycle <- dataset$cycle[keep]
	return(newname)
	}
	
	
split.chromosomes <- function(dataset){
	finalset <- lapply(levels(as.factor(dataset$LG)), function(x){
		keep <- which(dataset$LG==x)
		newname <- row.prune(dataset, keep)
		return(newname)
		})
	return(finalset)
}	
	
	
	
	
	
#prune out rows(markers) given a logical or vector of indexes
row.prune <- function(dataset, keep){
	newname <- dataset
	newname$calls <- dataset$calls[keep,]
    newname$scores <- dataset$scores[keep,]
    newname$LG <- dataset$LG[keep]
    newname$position <- dataset$position[keep]
    newname$markernames <- dataset$markernames[keep]
	newname$markercodes <- dataset$markercodes[keep]
	newname$alleles <- dataset$alleles[keep,]
	newname$genmap <- dataset$genmap[keep]
	newname$F2map <- dataset$F2map[keep]
	return(newname)
	}



#prune out columns(lines) given a logical or vector of indexes
column.prune <- function(dataset, keep){
	newname <- dataset
	newname$calls <- dataset$calls[,keep]
    newname$scores <- dataset$scores[,keep]
	newname$group <- dataset$group[keep]
	newname$status <- dataset$status[keep] 
	newname$cycle <- dataset$cycle[keep]
	newname$lines <- dataset$lines[keep]
	return(newname)
}


#split out all the groups
split.all <- function(dataset){
	dummy <- as.factor(paste(as.character(dataset$group), dataset$cycle, sep="_"))
	output <- list()
	for(i in 1:length(levels(dummy))){
		keep <- which(dummy%in%levels(dummy)[i])
		current_set <- column.prune(dataset=dataset, keep=keep)
		output[[i]] <- current_set
		names(output)[i] <- as.character(levels(dummy)[i])
		}
	return(output)
}

	
		

#merge a dataset that has been split

data.merge <- function(dataset1, dataset2){
	newdata <- dataset1
	newdata$calls <- cbind(dataset1$calls, dataset2$calls)
	newdata$scores <- cbind(dataset1$scores, dataset2$scores)
	newdata$status <- c(as.character(dataset1$status), as.character(dataset2$status))
	newdata$group <- c(as.character(dataset1$group), as.character(dataset2$group))
	newdata$lines <- c(dataset1$lines, dataset2$lines)
	newdata$cycle <- c(as.character(dataset1$cycle), as.character(dataset2$cycle))
	return(newdata)
	}	
	



	
		




