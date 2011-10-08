
parseGff3 <- function(filegff, fileg="genes.txt", filet="transcripts.txt", filee="exons.txt", nofiles=FALSE)
{
	
	cat ("Bonjour, starting GFF3 parsing for: ", filegff, "\n")
	gff <- read.table(filegff,sep='\t') 
	cat ( filegff, " file read happily! \n")
	

	idx <- which(gff[,3]=="gene")
	
	if (length(idx)==0) 
	{
			idx <- which(gff[,3]=="mRNA")
			gene_id <- apply(gff, 1, function(x) if (x[3]=="mRNA") rnaSeqMap:::.wytnij(x[9],"Parent")  else "---")  
		    gff <- cbind(gff,gene_id)
			genes <- gff[idx,c(1,4,5,7,10)]
			colnames(genes) <- c("chr","start","end","strand","gene_id")    
	}
    else
	{
		gene_id <- apply(gff, 1, function(x) if (x[3]=="gene") rnaSeqMap:::.wytnij(x[9],"ID=")  else "---")
		gff <- cbind(gff,gene_id)
		genes <- gff[idx,c(1,4,5,7,10)]
		colnames(genes) <- c("chr","start","end","strand","gene_id")
	}
	cat ( dim(genes)[1], " genes identified ! \n")
	
	
#transcript
	idx_t <- which(gff[,3]=="transcript")
	
	if (length(idx_t)==0)
	{
		idx_t <- which(gff[,3]=="mRNA")	
		transcript_id <- apply(gff, 1, function(x) if (x[3]=="mRNA") rnaSeqMap:::.wytnij(x[9],"ID=") else "---")
	    gff <- cbind(gff,transcript_id)

		parent <- apply(gff, 1, function(x) if (x[3]!="gene") rnaSeqMap:::.wytnij(x[9],"Parent") else "---")
	    gff <- cbind(gff,parent)

		transcripts <- gff[idx_t,c(1,4,5,7,11,12)]
	    colnames(transcripts) <- c("chr","start","end","strand","transcript_id","gene_id")
	}
	else
	{
		transcript_id <- apply(gff, 1, function(x) if (x[3]=="transcript") rnaSeqMap:::.wytnij(x[9],"ID=") else "---")
	    gff <- cbind(gff,transcript_id)

		parent <- apply(gff, 1, function(x) if (x[3]!="gene") rnaSeqMap:::.wytnij(x[9],"Parent") else "---")
	    gff <- cbind(gff,parent)
	    
		transcripts <- gff[idx_t,c(1,4,5,7,11,12)]
	    colnames(transcripts) <- c("chr","start","end","strand","transcript_id","gene_id")
	}
		
		
	cat ( dim(transcripts)[1], " transcripts identified ! \n")
	
#exons
	idx_e  <- which(gff[,3]=="exon")
	
	if (length(idx_e)==0) 
	{
	idx_e <- which(gff[,3]=="CDS" |  (gff[,3]!="mRNA" & gff[,3]!="gene")) 
	exon_id <- apply(gff, 1, function(x) if (x[3]!="mRNA" & x[3]!="gene" ) rnaSeqMap:::.wytnij(x[9],"ID=") else "---")
	gff <- cbind(gff,exon_id)
	exons <- gff[idx_e,c(1,4,5,7,12,13)]

	}
	else
	{
		exon_id <- apply(gff, 1, function(x) if (x[3]!="mRNA" & x[3]!="gene" ) rnaSeqMap:::.wytnij(x[9],"ID=") else "---")
		gff <- cbind(gff,exon_id)
		exons <- gff[idx_e,c(1,4,5,7,12,13)]
	}

	temptr <- transcripts
	colnames(temptr) <- c("chr","start","end","strand","parent","gene_id")
	exons <- merge(exons,temptr,by="parent")
	
	exons <- exons[,c(1,2,3,4,5,6,11)]
	colnames(exons) <- c("transcript_id","chr","start","end","strand","exon_id","gene_id")
		
	cat ( dim(exons)[1], " exons identified ! \n")
	
# a tera zapis
   if (!nofiles)
     {
	write.table (genes, file=fileg, sep="\t")
	write.table (transcripts, file=filet, sep="\t")
	write.table (exons , file=filee, sep="\t")
	cat (" output files written :", fileg, filet, filee, "\n")
     }	
      cat (" Gia finito, grazie per Sua attenzione!")
	
	out <- list()
	out[["genes"]] <- genes
	out[["transcripts"]] <- transcripts
	out[["exons"]] <- exons
	out
}


##################################################################################################################################################################



