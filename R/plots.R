#plots

#wykres funkcji pokrycia dla genu

plotGeneCoverage <- function(gene_id,ex)
{
	l<-length(ex)
        par(mfcol=c(l,1))

	colz <- topo.colors(length(ex))
	
 if (xmapConnected())
    {
		gd <- gene.details(gene_id)
		
		chr    <- gd$space
		start  <- gd$ranges@start
		end    <- gd$ranges@start+gd$ranges@width
		strand <- gd$strand 
		
     title <- paste(ex,gene_id)

	 for (i in 1:length(ex))
         {
	 		out<-regionCoverage(chr,start,end,strand,ex[i])
			m <- max(out[,2],ns.rm = FALSE) 
			mm <- max(m,na.rm=FALSE)
  		     
			plot(0,0,xlab="Position on the chromosome",
				  ylab="Nr of reads", 
				  xlim=c(start,end), 
				  ylim=c(0,mm) ,
				  bg="grey100", 
				  col.axis="tan4" )
			
			lines(out[,1],out[,2],col=colz[i],type="h",xlim=c(start,end),ylim=c(-1,mm))
			 
		 }
		legend("topleft", legend=c(paste("Gene: ",gene_id),paste("Date: ",format(Sys.time(),"%D")))) 
		legend("topright", legend=c(paste("sample",ex)), fill=(c(col=colz[ex])))

	}
	else
	{
		cat("No connection with database! \n")
    }
	
}



plotRegionCoverage <- function(chr, start, end, strand,ex)
{
   colz <- topo.colors(length(ex))
	l<-length(ex)
    par(mfcol=c(l,1))

	if (xmapConnected())
    {
	 	  
	 for (i in 1:length(ex))
         {
			 out<-regionCoverage(chr,start,end,strand, ex[i])
			 m <- max(out[,2],ns.rm = FALSE) 
			 mm <- max(m,na.rm=FALSE)

  		      plot(0,0,xlab="Position on the chromosome",
  	           ylab="Nr of reads", 
  	           xlim=c(start,end), 
  	           ylim=c(0,m) ,
  	           bg="grey100", 
  	           col.axis="tan4" )

  		     lines(out[,1],out[,2],col=colz[i],type="h")

         	}
      }
      else
      {
       
		  rs <- newSeqReads(chr,start,end,strand)
		  rs <- getBamData(rs,1:length(ex))
		  nd.cov <- getCoverageFromRS(rs,1:length(ex))
		  
		  distrCOVPlot(nd.cov,ex)
       }	
}


plotExonCoverage <- function(exon_id,ex)
{
	colz <- topo.colors(length(ex))
	ed <- exon.details(exon_id)
    l<-length(ex)
    par(mfcol=c(l,1))
	
	
	if (xmapConnected())
    {
	
		chr    <- ed$space
		start  <- ed$ranges@start
		end    <- ed$ranges@start+ed$ranges@width
		strand <- ed$strand 	 
		
	 for (i in 1:length(ex))
         {
			 out<-regionCoverage(chr,start,end,strand, ex[i])
			 m <- max(out[,2],ns.rm = FALSE) 
			 mm <- max(m,na.rm=FALSE)
			 
             plot(0,0,xlab="Position on the chromosome",
  	           ylab="Nr of reads", 
  	           xlim=c(start,end), 
  	           ylim=c(0,m+0.5) ,
  	           bg="grey100", 
  	           col.axis="tan4" )
  	           
  		     lines(out[,1],out[,2],col=colz[i],type="h")

         	}
      }
      else
      {
		cat("No connection with database! \n")
       }
	
	legend("topleft", legend=c(paste("Region - start:",start,"end:",end),paste("Date: ",format(Sys.time(),"%D"))),paste("Exon: ",exon_id)) 
	legend("topright", legend=c(paste("sample",ex)), fill=(c(col=colz[ex])))

}
#histogram 

plotCoverageHistogram <- function(chr,start,end,strand,ex, skip=10)
{
	if (xmapConnected())
    {
			out<-regionCoverage(chr,start,end,strand, ex)
			
			outp<- .Call("ghistogram",as.integer(out[,1]),as.integer(out[,2]),as.integer(skip))
            hNucleotide<- seq(start,end,skip)
			outl<-cbind(hNucleotide,outp)
			m <- max(outl[,2]) 
			
			title <- paste("Histogram for (",ex,chr,start,end,strand,") skip",skip)
			
			plot(0,0,
				 xlab="Position on the chromosome",
				 ylab="Nr of reads",
				 xlim=c(start,end), 
				 ylim=c(0,m) , 
				 main=title,
				 bg="grey100", 
				 col.axis="tan4" )
			
			lines(outl[,1],outl[,2],col="navajowhite4",bg="yellow", type="h",  las=1 ,cex=1.5)
		
	}
	else
	{
		
		rs <- newSeqReads(chr,start,end,strand)
		rs <- getBamData(rs,ex)
		nd.cov <- getCoverageFromRS(rs,ex)
		
		outp<- .Call("ghistogram",as.integer(c(start:end)),as.integer(as.vector(distribs(nd.cov)[[ex]])),as.integer(skip))

		hNucleotide<- seq(start,end,skip)
		outl<-cbind(hNucleotide,outp)
		m <- max(outp) 

		title <- paste("Histogram for (",ex,chr,start,end,strand,") skip",skip)
		
		plot(0,0,
			 xlab="Position on the chromosome",
			 ylab="Nr of reads",
			 xlim=c(start,end), 
			 ylim=c(0,m) , 
			 main=title,
			 bg="grey100", 
			 col.axis="tan4" )
		
		lines(outl[,1],outl[,2],col="navajowhite4",bg="yellow", type="h",  las=1 ,cex=1.5)
	}
}	


plotGeneExonCoverage <- function(gene_id,ex)
{
	l<-length(ex)
    par(mfcol=c(l+1,1))
    par(mar=c(2,2,0.5,0.5))
	
	
	colz <- topo.colors(length(ex))
    
 if (xmapConnected())
    {
              
		gd <- gene.details(gene_id)
		ed <- exon.details(gene.to.exon(gene_id))
		
		chr    <- gd$space
		start  <- gd$ranges@start
		end    <- gd$ranges@start+gd$ranges@width
		strand <- gd$strand
		
		estart <- ed$ranges@start
		ewidth <- ed$ranges@width
		
	   for (i in 1:length(ex))
         {
			 out<-regionCoverage(chr,start,end,strand, ex[i])
			 
  		     m <- max(out[,2],ns.rm = FALSE) 
  		     mm <- max(m,na.rm=FALSE)
			 
  		     plot(0,0,xlab="Position on the chromosome",
  	           ylab="Nr of reads", 
  	           xlim=c(start,end), 
  	           ylim=c(-1.5,mm) ,
  	           #main=title,
  	           bg="grey100", 
  	           col.axis="tan4")
			 lines(out[,1],out[,2],col=colz[i],type="l")
    	 }
   
		
		eksons <- rnaSeqMap:::.munion(gene_id)
	
	plot(0,0,xlab="Position on the chromosome",
		 ylab="Nr of reads", 
		 xlim=c(start,end), 
		 ylim=c(-0.5,1) ,
		 bg="grey100", 
		 col.axis="tan4")
	
  	 for (i in 1:length(eksons))
  	 {
  	 	rect(eksons@start[i],0.5,eksons@start[i]+eksons@width[i],1.5,col="red")
  	 }
	legend("bottomleft", legend=c(paste("Gene: ",gene_id))) 
	legend("bottom", legend=c(paste("Date: ",format(Sys.time(),"%D")))) 
	legend("bottomright", legend=c(paste("sample",ex),"exon"), fill=(c(col=colz[ex],"red"))) 
	}
	else
	{
	cat("No connection with database! \n")
	}
}

#splicing indeks

plotSI <- function(chr,start,end,strand,exp1,exp2)
{ 
	
	title <- paste("Splicing Index for (",exp1,exp2, chr,start,end,strand,")")
	outp <- NULL
	
	if (xmapConnected())
    {
		
		out1<-regionCoverage(chr,start,end,strand,exp1)
		m1 <- max(out1[,2],ns.rm = FALSE) 
		s1<-sum(out1[,2])
		
		out2<-regionCoverage(chr,start,end,strand,exp2)
		m2 <- max(out2[,2],ns.rm = FALSE) 
		s2<-sum(out2[,2])
		m <- max(m1,m2)
		
		c=(s1/s2)
		
		outp <- .Call("splicingind", as.numeric(out1[,2]), as.numeric(out2[,2]),c)		
		
		out <- cbind(out1[,1],outp)
		
		plot(0,0,xlab="Position on the chromosome",
			 ylab="Nr of reads", 
			 xlim=c(start,end), 
			 ylim=c(-1.5,1.5) ,
			 main=title,
			 bg="grey100", 
			 col.axis="tan4" )
		
		out1[,2] <- out1[,2]/m
		out2[,2] <- out2[,2]/m
		
		lines(out1[,1],out1[,2],col="blue",type="l")
		lines(out2[,1],out2[,2],col="red",type="l")
		lines(out[,1],out[,2],col="lawngreen",type="l")
		
		legend("topleft",
			   legend=c(paste("Sample",exp1),paste("Sample",exp2),"Splicing index"),
			   fill=c("blue","red","lawngreen"))
		legend("bottom",legend=c("The coverage value is normalised by max value"))
	}
	else
	{
		
		rs <- newSeqReads(chr,start,end,strand)
		rs <- getBamData(rs,c(exp1,exp2))
		
		nd.cov <- getCoverageFromRS(rs,1:2)
		m1 <- max(as.vector(distribs(nd.cov)[[1]]),ns.rm = FALSE) 
		s1<-sum(as.vector(distribs(nd.cov)[[1]]))
		
		m2 <- max(as.vector(distribs(nd.cov)[[2]]),ns.rm = FALSE) 
		s2<-sum(as.vector(distribs(nd.cov)[[2]]))
		
		m <- max(m1,m2)
		c=(s1/s2)
		
		outp <- .Call("splicingind", as.numeric(as.vector(distribs(nd.cov)[[1]])), as.numeric(as.vector(distribs(nd.cov)[[2]])),c)		
		
		out <- cbind(c(start:end),outp)
		
		plot(0,0,xlab="Position on the chromosome",
			 ylab="Nr of reads", 
			 xlim=c(start,end), 
			 ylim=c(-1.5,1.5) ,
			 main=title,
			 bg="grey100", 
			 col.axis="tan4" )
		
		o1 <- as.numeric(as.vector(distribs(nd.cov)[[1]]))/m
		o2 <- as.numeric(as.vector(distribs(nd.cov)[[2]]))/m
		
		lines(c(start:end),o1,col="blue",type="l")
		lines(c(start:end),o2,col="red",type="l")
		lines(c(start:end),out[,2],col="lawngreen",type="l")
		
		legend("topleft",
			   legend=c(paste("Sample",exp1),paste("Sample",exp2),"Splicing index"),
			   fill=c("blue","red","lawngreen"))
		legend("bottom",legend=c("The coverage value is normalised by max value"))
		
		
	}
}
