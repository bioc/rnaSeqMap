########################################################################
#   Transformations of Nucleotide Distribution objects
#   Michal Okoniewski, Anna Lesniewska, 3.09.2010 ostatni update 18.09.2010
########################################################################

findRegionsAsND <- function(nd,mi, minsup=5)
{
  stopifnot ( is( nd, "NucleotideDistr"))
        regVal <- NULL

        regVal <- assayData(nd)[["distribs"]]
        regVal <- regVal-regVal

        exps <- dim(assayData(nd)[["distribs"]])[2]

        for (e in 1:exps)
        {
                reg <- .Call("regionmining",as.integer(nd@start:nd@end),as.integer(as.vector(get("distribs", envir=assayData(nd))[,e])),as.integer(mi), as.integer(minsup))
                if (length(reg) > 1)
                {

                        cc <- 1
                        while (cc <length(reg)){
                                regVal[(reg[cc]-nd@start):(reg[cc+1]-nd@start+1),e] <- mi
                                cc <- cc+2
                        }
                }
                else
                {
                }
        }
        nd <- newNuctleotideDistr(regVal,nd@chr, nd@start, nd@end, nd@strand, "REG", phenoData=nd@phenoData)
        nd
}
###############################################################

findRegionsAsIR <- function(nd,mi, minsup=5, exp)
{
	stopifnot ( is( nd, "NucleotideDistr"))
	
	regIR <- NULL
    starts<- NULL
	ends  <- NULL
	exps  <- dim(assayData(nd)[["distribs"]])[2]
    cc    <- 1	
	
	
	
		
		reg <- .Call("regionmining",as.integer(nd@start:nd@end),as.integer(as.vector(get("distribs", envir=assayData(nd))[,exp])),as.integer(mi), as.integer(minsup))
		
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

###############################################################


regionBasedCoverage <- function(nd, seqq=1:10, maxexp=20, minsup=5)
{
    maxq <- max(seqq)
    mi <- maxexp*seqq[1]/maxq
    ndout <- findRegionsAsND(nd, mi)
    for (i in 2:length(seqq))
    {
       mi <- maxexp*seqq[i]/maxq
       ndtmp <- findRegionsAsND(nd, mi, minsup)
       maxdistr <- pmax(distribs(ndout), distribs(ndtmp))
       assign("distribs", maxdistr, envir=ndout@assayData)
    }
  ndout
}


#  returns normalized NucleotideDistr object of coverages
#  r- preprepared vector of global sums of reads :)
normalizeBySum <- function(nd,r=NULL)
{
  stopifnot ( is( nd, "NucleotideDistr"))
  normalVal <- NULL
  normalVal <- assayData(nd)[["distribs"]]
  exps <- dim(assayData(nd)[["distribs"]])[2]
  query <- "call seq_count_exps;"

  if (is.null(r))
   {

        CountDB <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
        CountDB <- as.vector(CountDB[,1])
   }
  else
  {
  CountDB<-r
  }

        c = max(CountDB[1])
    for (i in 1:exps)
                {
                        normalVal[,i] <- round(normalVal[,i]/CountDB[i]*c)
                }


        nd <- newNuctleotideDistr(normalVal,nd@chr, nd@start, nd@end, nd@strand, "NORM", phenoData=nd@phenoData)
        nd
}


getFCFromND <- function(nd, exps)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  stopifnot( length(exps)==2)

  dd <- distribs(nd)
  dd <- dd[,exps]
  dfc <- apply(dd, 1, function(x) x[1]-x[2])
  ndout <- newNuctleotideDistr(dfc,nd@chr, nd@start, nd@end, nd@strand, "FC")
  ndout
}

getSIFromND <- function(nd, exps)
{ 
  stopifnot( is( nd, "NucleotideDistr" ) )
  stopifnot( length(exps)==2)
  
  d1 <- distribs(nd[,exps[1]]) 
  d2 <- distribs(nd[,exps[2]]) 
  s1<-sum(d1)
  s2<-sum(d2)
  c=((s1-s2))

  outp <- .Call("splicingind", as.numeric(d1), as.numeric(d2),c)
  ndout <- newNuctleotideDistr(outp,nd@chr, nd@start, nd@end, nd@strand, "SI")
  ndout
}


averageND <- function(nd, exps)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  stopifnot( length(exps)>0)
  d <- apply(distribs(nd)[,exps],1,mean)
  ndout <- newNuctleotideDistr(d,nd@chr, nd@start, nd@end, nd@strand, nd@type)
}#by samples

sumND <- function(nd, exps)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  stopifnot( length(exps)>0)
  d <- apply(distribs(nd)[,exps],1,sum)
  ndout <- newNuctleotideDistr(d,nd@chr, nd@start, nd@end, nd@strand, nd@type)
}#by samples


combineND <- function(nd1, nd2)
{
  stopifnot( is( nd1, "NucleotideDistr" ) )
  stopifnot( is( nd2, "NucleotideDistr" ) )
  d1 <- apply(distribs(nd1),1,mean)
  d2 <- apply(distribs(nd2),1,mean) 
  ndout <- newNuctleotideDistr(d1+d2,nd1@chr, nd1@start, nd1@end, nd1@strand, nd1@type)
}


log2ND <- function(nd)
{
  stopifnot( is( nd, "NucleotideDistr" ) )
  d <- apply(distribs(nd),c(1,2),function(x) log2(x+0.000000001))
  ndout <- newNuctleotideDistr(d,nd@chr, nd@start, nd@end, nd@strand, nd@type, phenoData=nd@phenoData)
}




