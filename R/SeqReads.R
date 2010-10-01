#####################################################################
# The class SeqReads keep the reads for a given region on a chromosome
#  for a number of experimnets. Thus the genomic coordinates 
#  must be defined on initialization/construction 
# It is not S4 class, it is a half-product for NucleotideDistribution, which is S4
#  AL,MO, 18.08.2010
#####################################################################

setClass( "SeqReads",
    representation(
      data="list", 
      chr ="character",
      start="numeric",
      end="numeric",
      strand="numeric")
)

# create SeqReads object from a matrix with 2 rows or list
newSeqReads <- function( chr, start, end, strand, datain=NULL)
{ 
   if (is.null(datain)) datain <- list()
   else (stopifnot(datain,"list"))

   chr <- as.character(chr)
   start <-  as.numeric(start)
   end <-  as.numeric(end)
   strand <-  as.numeric(strand)
   rs <- new ("SeqReads", 
         data=datain,              
         chr=chr, start=start, end=end, strand=strand)
  #validObject(rs)
   rs
}

# creating empty SR with gene coordinates
newSeqReadsFromGene <- function(g)
{
	gd <- gene.details(g)
        if (is.null(gd)) stop("Gene not found in the database.")
        g.st <- as.numeric(start(gd@ranges))
	g.end <- as.numeric(end(gd@ranges))
	ch <- gd@ranges@partitioning@NAMES
	g.ch <- .chromosome.number(ch)
	g.strand <- gd$strand

	rs <- newSeqReads(g.ch,g.st,g.end,g.strand)
        rs
}


setValidity("SeqReads", function(object)
{
 # czy dane w liscie sa macierzami readow czyli: n x 2 (start, end)
 if (length(data)<(-10)) #temporary -10 to switch off
 {
   for (i in 1:length(data))
      if (dim(data[[i]])[2]!=2 | !is.null(data)) stop("Reads matrices should have 2 columns.")
 }
}
)


addDataToReadset <- function(rs, datain, spl)
{
  stopifnot( is( rs, "SeqReads" ) )
  if (dim(datain)[2]!=2)  stop("Reads matrices should have 2 columns.")
  ll <- rs@data
  ll[[spl]] <- datain
  rs@data <- ll
  rs
}


addExperimentsToReadset <- function(rs,exps)
# e - identifier of the experiment in the database, an int
{
  rs@data <- vector("list",max(exps))
  stopifnot( is( rs, "SeqReads" ) )
  for (e in exps)
  {
    datain <- readsInRange (e, rs@chr, rs@start, rs@end, rs@strand)
    if (is.null(datain)) rs@data[[e]] <- t(as.data.frame(as.numeric(c(0,0)))) 
      #empty readstable is marked as (0,0) 
      # see the joke about programmer looking for elephants in Africa
    else rs@data[[e]] <- datain
  }
rs
}

getCoverageFromRS <- function(rs, exps)
# returns NucleotideDistr object of coverages
{
  stopifnot( is( rs, "SeqReads" ) )
  covVal <- NULL
  for (e in exps)
	{       if (rs@data[[e]][1,1]==0 & rs@data[[e]][1,2]==0) 
		       {
				   # checking if there is an elephant at the end of Africa...
				   vzero <- vector(length=rs@end-rs@start+1)
				   vzero[1:(rs@end-rs@start)] <- 0
				   covVal <- cbind(covVal, vzero, deparse.level=0)
			   } 
		else
		{
		xxx <- .covFun(rs@data[[e]],rs@start, rs@end)
		 covVal <- cbind(covVal,xxx , deparse.level=0)
		}
  }
  if (is.null(xmapcore:::.xmap.internals$initialised)) pd <- NULL
  else 
  {
  pd <- getExpDescription()[exps,]
  pd <- new("AnnotatedDataFrame", data=as.data.frame(pd))
  }
  nd <- newNuctleotideDistr(covVal,rs@chr, rs@start, rs@end, rs@strand, "COV", phenoData=pd)

  nd
}



  	



