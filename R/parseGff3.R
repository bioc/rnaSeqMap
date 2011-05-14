
parseGff3 <- function(filegff, fileg="genes.txt", filet="transcripts.txt", filee="exons.txt", nofiles=FALSE)
{
	
	cat ("Buongiorno, starting GFF3 parsing for: ", filegff, "\n")
	gff <- read.table(filegff,sep='\t') 
	cat ( filegff, " file read happily! \n")
	
	idx <- which(gff[,3]=="gene")
	gene_id <- apply(gff, 1, function(x) if (length(x[9]>4) & x[3]=="gene") substr(x[9],4,regexpr(";",x[9])-1) else "---")
	gff <- cbind(gff,gene_id)
	genes <- gff[idx,c(1,4,5,7,10)]
	cat ( dim(genes)[1], " genes identified ! \n")
	
	
#transcript
	idx_t <- which(gff[,3]=="mRNA")
	transcript_id <- apply(gff, 1, function(x) if (length(x[9]>4) & x[3]=="mRNA") substr(x[9],4,regexpr(";",x[9])-1) else "---")
	gff <- cbind(gff,transcript_id)
	pacid <- apply(gff, 1, function(x) if (length(x[9]>4) & x[3]!="gene") rnaSeqMap:::.wytnij(x[9],"pacid") else "---")
	gff <- cbind(gff,pacid)
	gene_name <- apply(gff, 1, function(x) if (length(x[9]>4) & x[3]=="mRNA") rnaSeqMap:::.wytnij(x[9],"Parent") else "---")
	gff <- cbind(gff,gene_name)
	transcripts <- gff[idx_t,c(1,4,5,7,11,12,13)]
	cat ( dim(transcripts)[1], " transcripts identified ! \n")
	
#exons
	idx_e <- which(gff[,3]=="CDS" & gff[,3]!="mRNA" & gff[,3]!="gene")
	exon_id <- apply(gff, 1, function(x) if (length(x[9]>4) & x[3]!="mRNA" & x[3]!="gene" ) substr(x[9],4,regexpr(";",x[9])-1) else "---")
	gff <- cbind(gff,exon_id)
	exons <- gff[idx_e,c(1,4,5,7,12,14)]
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



