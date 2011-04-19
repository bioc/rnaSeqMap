########################################################################
#   Transformations of Nucleotide Distribution objects
#   Michal Okoniewski, Anna Lesniewska, 3.09.2010 ostatni update 18.09.2010
########################################################################

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

###############################################################

findRegionsAsIR <- function(nd,mi, minsup=5, exp)
{
	stopifnot ( is( nd, "NucleotideDistr"))
	
	regIR <- NULL
    starts<- NULL
	ends  <- NULL
	
	exps <- length(distribs(nd))
	cc    <- 1	
	
	reg <- .Call("regionmining",as.integer(nd@start:nd@end),as.integer(RleList2matrix(nd@data[exp]))	,as.integer(mi), as.integer(minsup))
		
		if (length(reg) > 1)
		{
			
			while (cc <length(reg))
			{
				starts <- rbind(starts,(reg[cc]))
				ends   <- rbind(ends,(reg[cc+1]))
				cc <- cc+2
			}	
		}
	

	regIR<-IRanges(starts,ends)
	regIR
}

###################################################################
# jeszcze nie poprawiona

regionBasedCoverage <- function(nd, seqq=1:10, maxexp=20, minsup=5)
{
    maxq <- max(seqq)
    mi <- maxexp*seqq[1]/maxq
   
	ndout <- findRegionsAsND(nd, mi)
	#ndout lista z regionami
	
    for (i in 2:length(seqq))
    {
       mi <- maxexp*seqq[i]/maxq
       ndtmp <- findRegionsAsND(nd, mi, minsup)
	   maxdistr <- pmax(RleList2matrix(ndout@data),RleList2matrix(ndtmp@data))
	   
#assign("distribs", maxdistr, envir=ndout@assayData)
# Nie wiem jak obejsc assayData    
	}
  ndout
}

#  returns normalized NucleotideDistr object of coverages
#  r- preprepared vector of global sums of reads :)

##################################################################

normalizeBySum <- function(nd,r=NULL)
{
  stopifnot ( is( nd, "NucleotideDistr"))
  normalVal <- NULL
  normalVal <- RleList2matrix(nd@data)
  normalVall<- NULL	
  exps <- length(distribs(nd))
  query <- "call seq_count_exps;"

  if (is.null(r))
   {

        CountDB <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
        CountDB <- as.vector(CountDB[,1])
   }
  else
  {
  CountDB <- r
  }

        c = max(CountDB[1])
    
	for (i in 1:exps)
                {
                        normalVal[,i] <- round(normalVal[,i]/CountDB[i]*c)
						normalVall <- c(normalVall,as(normalVal[,i],"Rle"))
				}
	
        nd <- newNuctleotideDistr(normalVall,nd@chr, nd@start, nd@end, nd@strand, "NORM", phenoData=nd@phenoData)
        nd
}

getFCFromND <- function(nd, exps)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  stopifnot(  length(exps) == 2 )
	
  dd <- RleList2matrix(nd@data)
  dd <- dd[,exps]
  dfc <- apply(dd, 1, function(x) x[1]-x[2])
  dfcl <- list(as(dfc,"Rle"))	
  ndout <- newNuctleotideDistr(dfcl,nd@chr, nd@start, nd@end, nd@strand, "FC")
  ndout
}


getSIFromND <- function(nd, exps)
{ 
  stopifnot( is( nd, "NucleotideDistr" ) )
  stopifnot( length(exps)==2)
  
  d1 <- RleList2matrix(nd@data[1])
  d2 <- RleList2matrix(nd@data[2])
  s1<-sum(d1)
  s2<-sum(d2)
  c=((s1-s2))

  outp <- .Call("splicingind", as.numeric(d1), as.numeric(d2),c)
  outpl <- list(rle(outp))	
	
  ndout <- newNuctleotideDistr(outpl,nd@chr, nd@start, nd@end, nd@strand, "SI")
  ndout
}

averageND <- function(nd, exps)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  stopifnot( length(exps)>0)
 
  d <- apply(RleList2matrix(nd@data[exps]),1,mean)
  dl <- list(rle(d))
  ndout <- newNuctleotideDistr(dl,nd@chr, nd@start, nd@end, nd@strand, nd@type)
}#by samples


sumND <- function(nd, exps)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  stopifnot( length(exps)>0)
  d <- apply(RleList2matrix(nd@data[exps]),1,sum)
  dl <- list(rle(d))	
  ndout <- newNuctleotideDistr(dl,nd@chr, nd@start, nd@end, nd@strand, nd@type)
  ndout	
}#by samples

combineND <- function(nd1, nd2)
{
  stopifnot( is( nd1, "NucleotideDistr" ) )
  stopifnot( is( nd2, "NucleotideDistr" ) )
  d1 <- apply(RleList2matrix(nd1@data),1,mean)
  d2 <- apply(RleList2matrix(nd2@data),1,mean) 
  com <- list(rle(d1+d2))	
  ndout <- newNuctleotideDistr(com,nd1@chr, nd1@start, nd1@end, nd1@strand, nd1@type)
  ndout	
}


log2ND <- function(nd)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  dl <- list()
	d <- apply(RleList2matrix(nd@data),c(1,2),function(x) log2(x+0.000000001)) 

	for (i in 1:length(distribs(nd)))
	{
	dl <- c(dl,as((d[,i]),"Rle"))
	}
	
  ndout <- newNuctleotideDistr(dl,nd@chr, nd@start, nd@end, nd@strand, nd@type, phenoData=nd@phenoData)
  ndout	
}




