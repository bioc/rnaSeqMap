##################
#Last update, MO 4.4.2011
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




###############################################################################

getCoverageFromRS <- function (rs, exps)

{
    stopifnot(is(rs, "SeqReads"))
    covVal <- list()
for (e in exps) {

        if (rs@data[[e]][1, 1] == 0 & rs@data[[e]][1, 2] == 0) {
            vzero <- vector(length = rs@end - rs@start + 1)
            vzero[1:(rs@end - rs@start)] <- 0
            covVal <- c(covVal, as(vzero,"Rle"))
}
        else
{
            xxx <- rnaSeqMap:::.covFun(rs@data[[e]], rs@start, rs@end)
covVal <- c(covVal, as(xxx,"Rle")) }
    }

   pd <- phenoData(rs)
   nd <- newNuctleotideDistr(covVal, rs@chr, rs@start, rs@end,
  rs@strand, "COV", phenoData = pd)
   nd
}


#### zamiana Rle na macierz

RleList2matrix <- function(list)
{
        stopifnot(is(list,"list"))

        for (i in length(list))
        {
                m <- max(length(list[[i]]))
        }

        macierz <- matrix(0,m,length(list))

        for (e in 1:length(list))
        {
                macierz[,e] <- as.vector(list[[e]])
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


