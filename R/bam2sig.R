bam2sig <- function(annotlib, covdesc="covdesc", species = NULL, level = "gene") 
# annot is for now an annot library *.db
{
	tout <- NULL
	
#library(annotlib)
	
    ConStr <- paste(substr(annotlib,1,(regexpr(".db",annotlib,T)-1)),"_dbconn",sep="")
	
	a_con <- get(ConStr, envir=.GlobalEnv)
	
	if (level=="transcript") query <- "SELECT chr as CHROMOSOME,START,END,STRAND,GENE_ID FROM TRANSCRIPTS" else if (level=="exon")	query <- "SELECT chr as CHROMOSOME,START,END,STRAND,GENE_ID FROM EXONS" else query <- "SELECT chr as CHROMOSOME,START,END,STRAND,GENE_ID FROM GENES"
	
	annot = dbGetQuery(a_con(), query)
	
	numgenes <- dim(annot)[1]
	cvd <- read.table(covdesc)
	numbams <- dim(cvd)[1] 
	
	if(is.null(species))
	{
		cat("Chromosome mapping not used for this run \n")
		
		for (i in 1:numgenes)
		{  
			ch <- annot[i,"CHROMOSOME"];
			st <- as.numeric(annot[i,"START"]);
			en <- as.numeric(annot[i,"END"]);
			str <- as.numeric(annot[i,"STRAND"]);
			rs <- rnaSeqMap:::newSeqReads(ch,st, en, str);
			rs <- rnaSeqMap:::getBamData(rs,1:numbams)
			counts <- lapply(rs@data, rnaSeqMap:::.countz) 
			tout <- rbind(tout, counts)
			
		} 
	}
	else
	{
		setSpecies(species)
		
		for (i in 1:numgenes)
		{  
			ch <- rnaSeqMap:::.chr.convert((annot[i,"CHROMOSOME"]));
			st <- as.numeric(annot[i,"START"]);
			en <- as.numeric(annot[i,"END"]);
			str <- as.numeric(annot[i,"STRAND"]);
			rs <- rnaSeqMap:::newSeqReads(ch,st, en, str);
			rs <- rnaSeqMap:::getBamData(rs,1:numbams)
			counts <- lapply(rs@data, rnaSeqMap:::.countz) 
			tout <- rbind(tout, counts)
		
		} 
	}
	
   if (level=="transcript") rownames(tout) <- annot[,"TRANSCRIPT_ID"] else if (level=="exon") rownames(tout) <- annot[,"EXON_ID"] else rownames(tout) <- annot[,"GENE_ID"]
   tout
 
}

#counts2sig <- function(compAttr="group", compVal=NULL)
#{
# preparing for DESeq
#    cvd <- read.table("covdesc")
#    conds <- factor(cvd[,compAttr])
#    cds <- newCountDataSet( in.deseq.table, conds )
#    cds <- estimateSizeFactors( cds )
#    cds <- estimateVarianceFunctions( cds , pool=T)
#    if (length(compValues)!=2 ) 
#     {
#       v1 <- conds[1]
#       v2 <- conds[2]
#       warning("Using two first values of the decisive attribute for comparisons!")
#     }
#     else
#     {
#       v1 <- compValues[1]
#       v2 <- compValues[2]
#     }
#    res <- nbinomTest( cds, v1, v2)
#    res
#}

.countz <- function(x)
{
  out <- dim(x)[1]
  if (out==1 & x[1,2]==0) out=0  # case of (0,0) read
  out
}

