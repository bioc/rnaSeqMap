#####################################################################
# The class SeqReads keep the reads for a given region on a chromosome
#  for a number of experimnets. Thus the genomic coordinates 
#  must be defined on initialization/construction 
# NucleotideDistribution, which is S4
#  AL,MO, 18.08.2010 update 5.10.2010, re-design: Mar/Apr 2011
#####################################################################

setClass( "SeqReads",
    contains= "eSet",
    representation(
      data="list", 
      chr ="character",
      start="numeric",
      end="numeric",
      strand="numeric")
)

# create SeqReads object from a matrix with 2 rows or list
newSeqReads <- function( chr, start, end, strand, datain=NULL, phenoData=NULL, featureData=NULL, covdesc=NULL)
{ 
   if (is.null(datain)) datain <- list()
   else (stopifnot(datain,"list"))

   chr <- as.character(chr)
   start <-  as.numeric(start)
   end <-  as.numeric(end)
   strand <-  as.numeric(strand)
   	
   if( !is.null( covdesc ) )
    {
     cvd <- read.table(covdesc)
     if( is.null( phenoData ) )
     phenoData <- as(cvd, "AnnotatedDataFrame")
    }
   else phenoData <- annotatedDataFrameFrom(phenoData, byrow=FALSE)

   rs <- new ("SeqReads", 
         data=datain,              
         chr=chr, start=start, end=end, strand=strand, phenoData=phenoData)
  #validObject(rs)
   rs
}

getBamData <- function( rs, exps=NULL, files=NULL, covdesc="covdesc")
{
  if (!is.null(files))
  {
    bams <- files
    phenoData(rs) <- as(as.data.frame(files), "AnnotatedDataFrame")
  }
  else 
  {
  cvd <- read.table(covdesc)
  bams <-rownames(cvd)
  phenoData(rs) <- as(cvd, "AnnotatedDataFrame")
  }
  
  if (length(bams)<1) stop("No .bam files?")
  if (!is.null(exps)) bams <- bams[exps]
   chr <- paste("chr", as.character(rs@chr), sep="")
   start <-  as.numeric(rs@start)
   end <-  as.numeric(rs@end)
   strand <-  as.numeric(rs@strand)
   if (strand==1)  strand <- "+"
   if (strand==-1)  strand <- "-"
   gr <- GRanges(seqnames =  chr,
                 ranges = IRanges(start,end),
                 strand =  strand)
   attrs <- c("strand", "pos", "qwidth")
   param <- ScanBamParam(which = gr, what=attrs)
   for (i in 1:length(bams))
   {
    # cat (bams[i],"  \n")
    # cat(str(param))
     outbam <- scanBam(bams[i],index=bams[i],param = param)[[1]]
     idx <- which(outbam$strand==strand)
     if (length(idx)>0)
     {
        ttt <- cbind (outbam$pos[idx],outbam$pos[idx] + outbam$qwidth[idx] )
        rs@data[[i]] <- ttt
     #   cat (dim(ttt))
     }
     else rs@data[[i]] <- t(as.data.frame(as.numeric(c(0,0))))
   }
   rs
}


addBamData <- function( rs, file, exp, phenoDesc=NULL)
{
  bams <- file
  if (length(bams)<1) stop("No .bam files?")
   chr <- paste("chr", as.character(rs@chr), sep="")
   start <-  as.numeric(rs@start)
   end <-  as.numeric(rs@end)
   strand <-  as.numeric(rs@strand)
   if (strand==1)  strand <- "+"
   if (strand==-1)  strand <- "-"
   gr <- GRanges(seqnames =  chr,
                 ranges = IRanges(start,end),
                 strand =  strand)
   attrs <- c("strand", "pos", "qwidth")
   param <- ScanBamParam(which = gr, what=attrs)
   #cat(str(param), )

     outbam <- scanBam(bams,index=bams,param = param)[[1]]
     idx <- which(outbam$strand==strand)
     if (length(idx)>0)
     {
        ttt <- cbind (outbam$pos[idx],outbam$pos[idx] + outbam$qwidth[idx] )
        rs@data[[exp]] <- ttt
     }
     else rs@data[[exp]] <- t(as.data.frame(as.numeric(c(0,0))))
   # adds only id size matches!
    if (!is.null(phenoDesc) & dim(phenoData(rs))[2]==length(phenoDesc) )  phenoData(rs) <- rbind(phenoData(rs), phenoDesc)
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
    datain <- readsInRange (rs@chr, rs@start, rs@end, rs@strand,e)
    if (is.null(datain)) rs@data[[e]] <- t(as.data.frame(as.numeric(c(0,0)))) 
      #empty readstable is marked as (0,0) 
      # see the joke about programmer looking for elephants in Africa
    else rs@data[[e]] <- datain
  }
rs
}






