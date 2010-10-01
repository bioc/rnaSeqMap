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
	query <- paste("call geneInChromosome(",chr,",",start,",",end,",",strand,");", sep="")
	out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con) 
	out
}

spaceInChromosome <- function(chr, start, end, strand)
{
	
	
	query <- paste("call spaceInChromosome(",chr,",",start,",",end,",",strand,");", sep="")
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
     name <- "homo_sapiens" 
     cat("As you left the parameter blank, I set the default species to homo_sapiens.\nBut please remember: \n\"Humanism is a superstition.\" J.M.Bochenski \n")
   }
   assign("species",name,envir=xmapcore:::.xmap.internals)
}


.chromosome.number <- function(ch)
{ 
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
 
