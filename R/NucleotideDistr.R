##################
#Last update, MO 18.4.2011
#
#

setClass( "NucleotideDistr",
contains = "eSet",
representation = representation(
chr ="character",
start="numeric",
end="numeric",
strand="numeric", 
type="character", 
data="list"),
)


newNuctleotideDistr <- function(distribs, chr, start, end, strand, type="UNKNOWN", phenoData=NULL, featureData=NULL)
{
	
	chr <- as.character(chr)
	start <-  as.numeric(start)
	end <-  as.numeric(end)
	strand <-  as.numeric(strand)
	type <- as.character(type)
	data <- as.list(distribs)
	
	if( is.null( phenoData ) )
	phenoData <- annotatedDataFrameFrom(distribs, byrow=FALSE )
	if( is.null( featureData ) )
	featureData <- annotatedDataFrameFrom(distribs, byrow=TRUE )
	
	nd <- new ("NucleotideDistr", 
			   data=data,               
			   chr=chr, start=start, end=end, strand=strand,
			   phenoData=phenoData, featureData=featureData )
	
#   validObject(nd)
	nd
}


distribs <- function(nd)
{
	stopifnot( is( nd, "NucleotideDistr" ) )
# na fjuczer trza rozpakowac i zrobic macierz jedna cbindem
	nd@data
}

# get a single distribution as a object Rle

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
	as(di.out,"Rle")
}


#method for getting data from ND withou @ operator/ampliQueso reviewer request/ MW
#16.09.2013
getData <- function(iND){
		
   return(iND@data)	
}

setSAXPYData <-function(iND1,iParam,i){

	iND1@data[[i]] <- iND1@data[[i]]*iParam 
}

setData <-function(iND1,iND2) {

	iND1<-iND2
}
###############################################################################
# new version using GAlignments in the input in RS. Output and ND content - no change
# update: 28 06 2012, MO
getCoverageFromRS <- function (rs, exps) 
{
    stopifnot(is(rs, "SeqReads"))
    covVal <- list()
    for (e in exps) {
    	
    	tmp <- coverage(granges(rs@data[[e]]))[[rs@chr]]
    	tmp <- tmp[rs@start:rs@end]
    	covVal <- c(covVal,tmp)
        }
    pd <- phenoData(rs)
    nd <- newNuctleotideDistr(covVal, rs@chr, rs@start, rs@end, rs@strand, "COV", phenoData = pd)
    nd
}


#### zamiana Rle na macierz
# update: 28 06 2012, MO

RleList2matrix <- function (ll) 
{
    stopifnot(is(ll, "list"))
    for (i in 1:length(ll)) {
        m <- max(length(ll[[i]]))
    }
    macierz <- matrix(0, m, length(ll))
    for (e in 1:length(ll)) {
        macierz[, e] <- as.vector(ll[[e]])
    }
    macierz
}



####################################################################


findRegionsAsND <- function(nd,mi,minsup=5)
{
	stopifnot(is(nd, "NucleotideDistr"))
	
	regVal <- NULL
	regVal <- RleList2matrix(nd@data)
	regVal <- regVal - regVal
	regVall <- list()
	
	exps <- dim(RleList2matrix(nd@data))[2]
	
	for (e in 1:exps) 
	{
		reg <- .Call("regionmining", as.integer(nd@start:nd@end), as.integer(RleList2matrix(nd@data[e])),
					 as.integer(mi), as.integer(minsup))
		
		if (length(reg) > 1) 
		{
			cc <- 1
			while (cc < length(reg))
			{
				regVal[(reg[cc] - nd@start):(reg[cc + 1] - nd@start + 1), e] <- mi
				cc <- cc + 2
			}
		}
	regVall <- c(regVall,as(regVal[,e],"Rle"))
	}
	
	nd <- newNuctleotideDistr(regVall, nd@chr, nd@start, nd@end, 
							  nd@strand, "REG", phenoData = nd@phenoData)
	nd
}



########################################################################

#rs <- newSeqReads(1,14000,30000,-1) 
#rs <- addBamData(rs,bamFile, 1)
#rs <- addExperimentsToReadset(rs,1:6)
#nd <- newgetCoverageFromRS(rs,1:6)
#nd.reg <- newfindRegionsAsND(nd,10)


