# ПОИСК ORF

revcomp <- function(DNA){
	DNA <- strsplit(DNA,split="")[[1]]
	DNA <- rev(DNA)
	nucls <- c("A","C","G","T")
	names(nucls) <- c("T","G","C","A")
	DNA_new <- c()
	for(nucl in DNA){
		DNA_new <- append(DNA_new,nucls[nucl])
	}
	return(paste(DNA_new,collapse=""))
}


# 1. Simulator

sim.nucls <- function(len=100){
	# creates fake nucleotides
	if(is.numeric(len)){
		nucls <- c('A','C','G','T')
		return(paste(sample(nucls,len,TRUE),collapse=""))
	}else{
		print("ERROR: wrong type of input")
	}
}


# СОЗДАЕМ ОТРЕЗОК ОДИНАКОВЫЙ У ВСЕХ
set.seed(12345)
dat <- sim.nucls(10000)


# 2. Triplets

find_tripl <- function(DNAseq,frame=1){
	# find all triplets in 1st frame
	
	if(frame==3){
		frame <- 0
	}
	
	if(is.character(DNAseq)){
		triplets <- c()
		for(i in 1:nchar(DNAseq)){
			if(i %% 3 == frame){
				triplets <- append(triplets,substr(DNAseq,i,i+2))
			}
		}
		return(triplets)
	}else{
		print("ERROR: wrong input data for find_tripl")
	}
}

# 3. START/STOP CODONS

is.start <- function(triplet){
	if(is.character(triplet)){
		if(triplet == "ATG"){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
}

is.stop <- function(triplet){
	if(is.character(triplet)){
		if(triplet %in% c("TAA","TAG","TGA")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
}

 
# 4. ORF
save_orfs <- function(triplets,minlen=200){
	ORFs <- c()
	ORF <- c()
	
	#### STATISTICS
	counter <- 0
	ORF_lens <- c()
	starts <- c()
	stops <- c()
	stops_id <- c()
	
	check <- FALSE

	for(i in 1:length(triplets)){
		if(!(check) && is.start(triplets[i])){
			ORF <- append(ORF,triplets[i])
			check <- TRUE
			start_id <- i
			next
		}
		if(is.stop(triplets[i])){
			ORF <- append(ORF,triplets[i])
			if(length(ORF) * 3 >= minlen){
				ORFs <- append(ORFs,paste(ORF,collapse=""))
				counter <- counter + 1
				ORF_lens <- c(ORF_lens,length(ORF)*3)
				starts <- c(starts, start_id)
				stops <- c(stops, i)
				stops_id <- c(stops_id, triplets[i])
			}
			ORF <- c()
			check <- FALSE
		}
		if(check){
			ORF <- append(ORF,triplets[i])
		}
	}
	if(length(ORF)>=minlen){
		ORFs <- append(ORFs,paste(ORF,collapse=""))
		counter <- counter + 1
		ORF_lens <- c(ORF_lens,length(ORF)*3)
		starts <- c(starts, start_id)
		stops <- c(stops, NA)
		stops_id <- c(stops_id, NA)
	}
	return(list(ORFs,1:counter, ORF_lens, starts, stops, stops_id))
} 

remove_orf  <- function(ORFdata){
	removal_orf <- c()
	forw_start_sort <- order(ORFdata$START[ORFdata$FRAME>0])
	rev_stop_sort <- order(ORFdata$STOP[ORFdata$FRAME<0])+ length(forw_start_sort)

	check <- 0
	checkid <- NA
	for(orf in forw_start_sort){
		if((ORFdata$START[orf] < check) && checkid){
			if((ORFdata$LEN[orf] >= ORFdata$LEN[checkid])|| checkid %in% removal_orf){
				removal_orf <- c(removal_orf,checkid)
			}else{
				removal_orf <- c(removal_orf,orf)
			}
		}
		check <- ORFdata$STOP[orf]
		checkid <- orf
	}
	
	check <- 0
	checkid <- NA
	for(orf in rev_stop_sort){
		if((ORFdata$STOP[orf] < check) && checkid){
			if((ORFdata$LEN[orf] >= ORFdata$LEN[checkid]) || checkid %in% removal_orf){
				removal_orf <- c(removal_orf,checkid)
			}else{
				removal_orf <- c(removal_orf,orf)
			}
		}
		check <- ORFdata$START[orf]
		checkid <- orf
	}

	if(!is.null(removal_orf)){
		removal_orf <- unique(sort(removal_orf))
		for(k in names(ORFdata)){
			ORFdata[[k]] <- ORFdata[[k]][-removal_orf]
			}	
	ORFdata$ID <- 1:length(ORFdata$ID)
	}

	return(ORFdata)
}



#5. Find  ORFs

orf_finder <- function(DNA, minlen=50){
	
	#Statistics
	ALLDATA <- list(c(),c(),c(),c(),c(),c())	

	
	if(is.character(DNA)){
		frames  <- c(1,2,3,-1,-2,-3)
		for(frameid in frames){
			if(frameid > 0){
				triplets <- find_tripl(DNA,frame= frameid)
				ORFs <- save_orfs(triplets,minlen=minlen)
				ALLDATA[[1]]<- c(ALLDATA[[1]], ORFs[[1]])
				ALLDATA[[2]]<- c(ALLDATA[[2]], length(ALLDATA[[2]])+ORFs[[2]])
				ALLDATA[[3]]<- c(ALLDATA[[3]], ORFs[[3]])
				ALLDATA[[4]]<- c(ALLDATA[[4]], ORFs[[4]]*3-3+ frameid)
				ALLDATA[[5]]<- c(ALLDATA[[5]], ORFs[[5]]*3-1+ frameid)
				ALLDATA[[6]]<- c(ALLDATA[[6]], rep(frameid,length(ORFs[[1]])))

			}else{
				triplets <- find_tripl(revcomp(DNA),frame=abs(frameid))
				ORFs <- save_orfs(triplets,minlen=minlen)
				ALLDATA[[1]]<- c(ALLDATA[[1]], ORFs[[1]])
				ALLDATA[[2]]<- c(ALLDATA[[2]], length(ALLDATA[[2]])+ORFs[[2]])
				ALLDATA[[3]]<- c(ALLDATA[[3]], ORFs[[3]])
				ALLDATA[[4]]<- c(ALLDATA[[4]], nchar(DNA)-ORFs[[4]]*3+4-abs(frameid))
				ALLDATA[[5]]<- c(ALLDATA[[5]], nchar(DNA)-ORFs[[5]]*3+2-abs(frameid))
				ALLDATA[[6]]<- c(ALLDATA[[6]], rep(frameid,length(ORFs[[1]])))
			}
		}
				
		names(ALLDATA) <- c("ORF","ID","LEN","START","STOP","FRAME")
		return(remove_orf(ALLDATA))
		
	}else{
		print("ERROR: wrong format for input data")
	}
}








