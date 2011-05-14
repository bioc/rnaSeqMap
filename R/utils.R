#######################################################################
#utils for rna seq data
#AL, MO, ostatnie update 18.09.2010
#######################################################################
.sigRegions <- function(ex, chr,start,end,strand,mi, minsup, db = "FALSE" )
{
	
	if (db=='FALSE')
    {
		for (i in 1:length(ex))
		{
         
			
			cov <- regionCoverage(ex,chr,start,end,strand)

			outp<- .Call("regionmining",cov[,1],cov[,2],mi, minsup)
			out<-cbind(outp[seq(1,length(outp)-1,2)],outp[seq(2,length(outp),2)])
			
		}
		out
		
	}
	else
	{
			query <- paste("call seq_window(",ex,",",chr,",",start,",",end,",",strand,",",mi,");", sep="")
    
			out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
    }
		out
}
 



#funkcja pokrycia
regionCoverage <- function(chr, start, end, strand, ex, db = "FALSE" )
{
	
	if (db=='FALSE')
    {
		for (i in 1:length(ex))
		{
         	query <- paste("call readsInRange(",ex,",",chr,",",start,",",end,",",strand,");", sep="")
  	 		reads <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
			
			Nucleotide<-c(start:end)
			b<-reads[,1]
			c<-reads[,2]
			
			Nr_reads<-.Call("gcoverage",Nucleotide,b,c)
			out<-cbind(Nucleotide,Nr_reads)
			 
		}
		out
	}
	else
	{
		for (i in 1:length(ex))
		{
   			query <- paste("call seq_coverage(",ex[i],",",chr,",",start,",",end,",",strand,");", sep="")
   			cat(query,"\n") 
   			out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)	
		
	     }
		out
	}
	
}

# coverage on a matrix of (start,end) reads
.covFun <- function(m, stR=NULL, endR=NULL)
{
   if (is.null(stR)) stR <- min(m[,1])
   if (is.null(endR)) endR <- max(m[,2])

	# diagnostic cats to see what we actually get - messy, but nice
        #cat(".covFun",stR," ",endR," dim",dim(m),"\n")
	#cat("\nm1:",m[,1],"\n")
	#cat("m2:",m[,2],"\n")
	#cat("zakresCovFun:",stR,"",endR,"\n")
	
	.Call("gcoverage",(stR:endR),(m[,1]),(m[,2]))

}


readsInRange <- function(chr, start, end, strand, ex)
{
	query <- paste("call readsInRange(",ex,",",chr,",",start,",",end,",",strand,");", sep="")
	temp <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
	if  (!is.null(temp)) colnames(temp) <- c("start","end")
	temp
}

#funkcja pokrycia dla genu

.geneCoverage <- function(ex, gene_id, db = "FALSE")
{
	
	gd <- gene.details(gene_id)
	
	chr    <- gd$space
	start  <- gd$ranges@start
	end    <- gd$ranges@start+gd$ranges@width
	strand <- gd$strand 	 
	
	if (db=='FALSE')
    {
		
		for (i in 1:length(ex))
		{
			out <- regionCoverage(ex,chr,start,end,strand)
		}
		out
	}
	else
	{
		for (i in 1:length(ex))
		{
			query <- paste("call seq_coverage(",ex[i],",",chr,",",start,",",end,",",strand,");", sep="")
			out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
		}
		out
	}
}


#funkcja pokrycia dla exonu

.exonCoverage <- function(ex,exon_id,db = "FALSE")
{

	ed <- exon.details(exon_id)
	
	chr    <- ed$space
	start  <- ed$ranges@start
	end    <- ed$ranges@start+ed$ranges@width
	strand <- ed$strand 	 
	
	if (db=='FALSE')
    {
		
		for (i in 1:length(ex))
		{
         	out <- regionCoverage(ex,chr,start,end,strand)
		}
		out	
	}
	else
	{
		for (i in 1:length(ex))
		{
  			query <- paste("call seq_coverage(",ex[i],",",chr,",",start,",",end,",",strand,");", sep="")
  			out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
		}
		out
	}
	
}

.munion <- function(gene_id)
# playing with exons on Iranges 
{
	ed <- exon.details(gene.to.exon(gene_id))

    estart <- ed$ranges@start
 	eranges     <- ed$ranges
     w1<-c(0)
     w2<-c(0)
     
     w <- IRanges(w1,w2)
    
 	 for (e in 1:length(estart))
  	 {
  	 	w <- union(w,eranges)
  	 }
w
}  	

.tunion <- function(transcript_id)
# playing with transcript in IRangges
{
	td <- exon.details(transcript.to.exon(transcript_id))
	
    estart <- td$ranges@start
 	eranges     <- td$ranges
	w1<-c(0)
	w2<-c(0)
	
	w <- IRanges(w1,w2)
    
	for (e in 1:length(estart))
	{
  	 	w <- union(w,eranges)
	}
	w
}  	


geneInChromosome <- function(chr, start, end, strand)
{
	query <- paste("call geneInChromosome('",as.character(chr),"',",start,",",end,",",strand,");", sep="")
	out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con) 
	out
}

spaceInChromosome <- function(chr, start, end, strand)
{
	
	
	query <- paste("call spaceInChromosome('",as.character(chr),"',",start,",",end,",",strand,");", sep="")
	out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con) 
	out
}


.readsForGene <- function(chr,start,end,strand)
{
	query <- paste("call readsForGene(",chr,",",start,",",end,",",strand,");", sep="")
	reads <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
reads
}

.countRSReads <- function(x)
{  
	if (dim(x)[1]==1 & x[1,1]==0) cnt <-0
	else cnt <- dim(x)[1]
	cnt
}

.rsCount <- function(g, exps)
{
	
	gd <- gene.details(g)
	g.st <- as.numeric(start(gd@ranges))
	g.end <- as.numeric(end(gd@ranges))
	ch <- gd@ranges@partitioning@NAMES
        g.ch <- .chromosome.number(ch)
	g.strand <- gd$strand

        		
	rs <- newSeqReads(g.ch,g.st,g.end,g.strand)
	rs <- addExperimentsToReadset(rs,exps)
	out <- unlist(lapply(rs@data, .countRSReads))
	out
}

# sets species name in case we do not connect to database
setSpecies <- function(name=NULL)
{
	if (is.null(name))
	{
		
		choices <- c(as.character("homo_sapiens"), as.character("rattus_norvegius"),as.character("mus_musculus"), as.character("methylobacterium_extorquens"))	
		r = menu(choices, title="Select species:")
		
		if (r == 1) 
			name <- "homo_sapiens"
		else if (r == 2)
		    name <- "rattus_norvegius"
		else if (r == 3)  
			name <- "mus_musculus"
		else if (r == 4)
		    name <- "methylobacterium_extorquens"
		else {
		    name <- "homo_sapiens" 
			cat("As you left the parameter blank, I set the default species to homo_sapiens.\nBut please remember: \n\"Humanism is a superstition.\" J.M.Bochenski \n")
		}}
		
	assign("species",name,envir=xmapcore:::.xmap.internals)
}



.chr.convert <- function(name)
{
	s <- get("species", envir=xmapcore:::.xmap.internals)
	
		switch(s,
			   "mus_musculus" = .mus_musculus(name),
			   "rattus_norvegius" = .rattus_norvegius(name),
			   "methylobacterium_extorquens" = .methylobacterium_extorquens(name),
			   .homo_sapiens(name))	
}
			
.homo_sapiens <- function(name){
			if (name=='1') out='chr1'
			else if (name=='2') out='chr2'
			else if (name=='3') out='chr3'
			else if (name=='4') out='chr4'
			else if (name=='5') out='chr5'
			else if (name=='6') out='chr6'
			else if (name=='7') out='chr7'
			else if (name=='8') out='chr8'
			else if (name=='9') out='chr9'
			else if (name=='10') out='chr10'
			else if (name=='11') out='chr11'
			else if (name=='12') out='chr12'
			else if (name=='13') out='chr13'
			else if (name=='14') out='chr14'
			else if (name=='15') out='chr15'
			else if (name=='16') out='chr16'
			else if (name=='17') out='chr17'
			else if (name=='18') out='chr18'
			else if (name=='19') out='chr19'
			else if (name=='20') out='chr20'
			else if (name=='21') out='chr21'
			else if (name=='22') out='chr22'
			else  if (name=='X') out='chrX'
			else if (name=='Y') out='chrY'
			else if (name=='MT') out='chrM'
			else out=30
			}
		
	 .rattus_norvegius <- function(name)
		{
			if (name=='1') out='chr1'
			else if (name=='2') out='chr2'
			else if (name=='3') out='chr3'
			else if (name=='4') out='chr4'
			else if (name=='5') out='chr5'
			else if (name=='6') out='chr6'
			else if (name=='7') out='chr7'
			else if (name=='8') out='chr8'
			else if (name=='9') out='chr9'
			else if (name=='10') out='chr10'
			else if (name=='11') out='chr11'
			else if (name=='12') out='chr12'
			else if (name=='13') out='chr13'
			else if (name=='14') out='chr14'
			else if (name=='15') out='chr15'
			else if (name=='16') out='chr16'
			else if (name=='17') out='chr17'
			else if (name=='18') out='chr18'
			else if (name=='19') out='chr19'
			else  if (name=='X') out='chrX'
			else if (name=='Y') out='chrY'
			else if (name=='MT') out='chrM'
			else out=30
		
		}
	
		.mus_musculus <- function(name)
		{
			if (name=='1') out='chr1'
			else if (name=='2') out='chr2'
			else if (name=='3') out='chr3'
			else if (name=='4') out='chr4'
			else if (name=='5') out='chr5'
			else if (name=='6') out='chr6'
			else if (name=='7') out='chr7'
			else if (name=='8') out='chr8'
			else if (name=='9') out='chr9'
			else if (name=='10') out='chr10'
			else if (name=='11') out='chr11'
			else if (name=='12') out='chr12'
			else if (name=='13') out='chr13'
			else if (name=='14') out='chr14'
			else if (name=='15') out='chr15'
			else if (name=='16') out='chr16'
			else if (name=='17') out='chr17'
			else if (name=='18') out='chr18'
			else if (name=='19') out='chr19'
			else  if (name=='X') out='chrX'
			else if (name=='Y') out='chrY'
			else if (name=='MT') out='chrM'
			else out=30
		}

		.methylobacterium_extorquens  <- function(name)
	{
		if (name=='META1_META1') out='META1'
		else if (name=='META2_META2') out='META2'      
		else  if (name=='p1META_p1META') out='p1META'
		else if (name=='p2META_p2META') out='p2META'
		else if (name=='p3META_p3META') out='p3META'
		else out=30
	}
###

.chromosome.number <- function(ch)
{
  if (class(ch)=="factor") ch <- as.vector(ch) 
  cho <- NULL
  s <- get("species", envir=xmapcore:::.xmap.internals)
  if (s=="homo_sapiens")
  {
	if (ch=="X") cho=23
	else if (ch=="Y") cho=24
	else if (ch=="MT") cho=25
        else cho <- as.numeric(ch)
  }
  else if (s=="mus_musculus")
  {
	if (ch=="X") cho=20
	else if (ch=="Y") cho=21
	else if (ch=="MT") cho=22
        else cho <- as.numeric(ch)
  }
  else if (s=="rattus_norvegicus")
  {
	if (ch=="X") cho=21
	else if (ch=="Y") cho=22
	else if (ch=="MT") cho=23
        else cho <- as.numeric(ch)
  }
  else cho <- as.numeric(ch)
  cho
}


.chromosome.char <- function(ch)
{
if (class(ch)=="factor") ch <- as.vector(ch) 
chr <- NULL	
s <- get("species", envir=xmapcore:::.xmap.internals)
	
	if (s=="homo_sapiens")
	{
		if (ch==23) chr=("X")
			else if (ch==24) chr=("Y")
			else if (ch==25) chr=("MT")
			else chr<- as.character(ch)
	}
	else 
	if (s=="mus_musculus")
	{
		if (ch==20) chr=("X")
		else if (ch==21) chr=("Y")
		else if (ch==22) chr=("MT")
        else chr <- as.character(ch)
	}
	else if (s=="rattus_norvegicus")
	{
		if (ch==21) chr=("X")
		else if (ch==22) chr=("Y")
		else if (ch==23) chr=("MT")
        else chr <- as.character(ch)
	}
	else chr <- as.character(ch)
chr	
}



# get experiment (all samples) description from xmap database
getExpDescription <- function()
{
        query <- paste("call showDescription();", sep="")
        description <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
        description
}

getSumsExp <- function()
{
   query <- "call seq_count_exps;"
   CountDB <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
   CountDB[,1]
}

xmapConnected <- function()
{
 exists("con",xmapcore:::.xmap.internals)
}
 

.wytnij <- function(text,pathern){
nn <- strsplit(text,";")[[1]] 
n2 <- grep(pathern,nn) 
n3 <- nn[n2]
out <- strsplit(n3,"=")[[1]][2]
out
}


