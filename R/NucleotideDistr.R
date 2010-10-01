setClass( "NucleotideDistr",
   contains = "eSet",
   representation = representation(
      chr ="character",
      start="numeric",
      end="numeric",
      strand="numeric", 
      type="character"),
      prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), NucleotideDistr = "1.1.0" ) ) )
)


newNuctleotideDistr <- function(distribs, chr, start, end, strand, type="UNKNOWN", phenoData=NULL, featureData=NULL)
{
   distribs <- as.matrix(distribs)
   chr <- as.character(chr)
   start <-  as.numeric(start)
   end <-  as.numeric(end)
   strand <-  as.numeric(strand)
   type <- as.character(type)
   if( is.null( phenoData ) )
      phenoData <- annotatedDataFrameFrom(distribs, byrow=FALSE )
   if( is.null( featureData ) )
      featureData <- annotatedDataFrameFrom(distribs, byrow=TRUE )
   nd <- new ("NucleotideDistr", 
         assayData=assayDataNew("environment",distribs=distribs),               
         chr=chr, start=start, end=end, strand=strand,
         phenoData=phenoData, featureData=featureData )
#   validObject(nd)
   nd
}


setValidity("NucleotideDistr", function(object)
{
   #if (dim(distribs)[1]<1) return("A proper matrix is needed.")
   #if (dim(distribs)[1]!= (end-start)) return("Data matrix not consistent with genome coordinates.")
  # check if pData has as many rows as assayData
  # check if pData has at least one column
}
)


distribs <- function(nd)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  get("distribs",envir=assayData(nd))
}

# get a single distribution as a vector
getDistr <- function(nd, i, a=NULL, b=NULL )
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  if (is.null(a)| is.null(b))
  {
    a <- nd@start
    b <- nd@end
  }
  if (b<a | a<nd@start | b>nd@end) stop ("Problem with coordinates.")
  di.out <- assayData(nd)[["distribs"]][(a-nd@start+1):(b-nd@start+1),i]
  di.out
}





