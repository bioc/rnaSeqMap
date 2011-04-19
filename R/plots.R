#plots

#wykres funkcji pokrycia dla genu

plotGeneCoverage <- function(gene_id,ex, db = "FALSE")
{
	l<-length(ex)
        par(mfcol=c(l,1))

	colz <- topo.colors(length(ex))
	gd <- gene.details(gene_id)

	chr    <- gd$space
	start  <- gd$ranges@start
	end    <- gd$ranges@start+gd$ranges@width
	strand <- gd$strand 	 
 
 if (db=='FALSE')
    {
     
     title <- paste(ex,gene_id)

	 for (i in 1:length(ex))
         {
	 		out<-regionCoverage(ex[i],chr,start,end,strand)
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
     }
	else
	{
		for (i in 1:length(ex))
         {
  		 query <- paste("call seq_coverage(",ex[i],",",chr,",",start,",",end,",",strand,");", sep="")
  	     out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
    
  		 m <- max(out[,2],ns.rm = FALSE)
   	     plot(0,0,xlim=c(start,end), ylim=c(0,m) )

         lines(out[,1],out[,2],col=colz[i],type="h")
         }
    }
	
	legend("topleft", legend=c(paste("Gene: ",gene_id),paste("Date: ",format(Sys.time(),"%D")))) 
	legend("topright", legend=c(paste("sample",ex)), fill=(c(col=colz[ex])))
}




plotRegionCoverage <- function(chr, start, end, strand,ex, db = "FALSE" )
{
    colz <- topo.colors(length(ex))
	l<-length(ex)
    par(mfcol=c(l,1))
	

	if (db=='FALSE')
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
       for (i in 1:length(ex))
		  {
   			query <- paste("call seq_coverage(",ex[i],",",chr,",",start,",",end,",",strand,");", sep="")
   			out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
   			m <- max(out[,2],na.rm=FALSE)
   			plot(0,0,xlim=c(start,end), ylim=c(0,m) )
   			lines(out[,1],out[,2],col=colz[i],type="h")
  			}
       }
	legend("topleft", legend=c(paste("Region - start:",start,"end:",end),paste("Date: ",format(Sys.time(),"%D")))) 
	legend("topright", legend=c(paste("sample",ex)), fill=(c(col=colz[ex])))
}


plotExonCoverage <- function(exon_id,ex,db = "FALSE")
{
	colz <- topo.colors(length(ex))
	ed <- exon.details(exon_id)
    l<-length(ex)
    par(mfcol=c(l,1))
	
	chr    <- ed$space
	start  <- ed$ranges@start
	end    <- ed$ranges@start+ed$ranges@width
	strand <- ed$strand 	 

	if (db=='FALSE')
    {
	 	  
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
		for (i in 1:length(ex))
  		   {
  			query <- paste("call seq_coverage(",ex[i],",",chr,",",start,",",end,",",strand,");", sep="")
  			out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
   
   			m <- max(out[,2],ns.rm = FALSE)
			mm <- max(m,na.rm=FALSE)
			plot(0,0,xlim=c(start,end), ylim=c(0,m) )
   			lines(out[,1],out[,2],col=colz[i],type="h")
           }
       }
	
	legend("topleft", legend=c(paste("Region - start:",start,"end:",end),paste("Date: ",format(Sys.time(),"%D"))),paste("Exon: ",exon_id)) 
	legend("topright", legend=c(paste("sample",ex)), fill=(c(col=colz[ex])))

}

#histogram 
plotCoverageHistogram <- function(chr,start,end,strand,ex, skip=10,db = "FALSE")
{
	l<-length(ex)
    par(mfcol=c(l,1))
	
	colz <- heat.colors(length(ex))

	if (db=='FALSE')
    {
		for (i in 1:length(ex))
		{
			out<-regionCoverage(chr,start,end,strand, ex[i])
			
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
	}
	else
  {
   query <- paste("call seq_histogram(",ex,",",chr,",",start,",",end,",",strand,",",skip,");", sep="")
    
   out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
   m <- max(out[,2],na.rm = FALSE)
   plot(0,0,xlim=c(start,end), ylim=c(0,m) )
   
   lines(out[,1],out[,2],col="navajowhite4",bg="yellow", type="h",  las=1 ,cex=1.5)
  }
}


plotGeneExonCoverage <- function(gene_id,ex, db = "FALSE")
{
	l<-length(ex)
    par(mfcol=c(l+1,1))
    par(mar=c(2,2,0.5,0.5))
	
	
	colz <- topo.colors(length(ex))
	gd <- gene.details(gene_id)
	ed <- exon.details(gene.to.exon(gene_id))

	chr    <- gd$space
	start  <- gd$ranges@start
	end    <- gd$ranges@start+gd$ranges@width
	strand <- gd$strand

    estart <- ed$ranges@start
 	ewidth <- ed$ranges@width
    
 if (db=='FALSE')
    {
              
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
   
	}
	else
	{ 
	for (i in 1:length(ex))
  		{
  			query <- paste("call seq_coverage(",ex[i],",",chr,",",start,",",end,",",strand,");", sep="")
   			out <- xmapcore:::.xmc.db.call(query, xmapcore:::.xmap.internals$con)
  			m <- max(out[,2],na.rm = FALSE)
  			mm <- max(m,na.rm=FALSE)
  			plot(0,0,xlim=c(start,end), ylim=c(-1,mm) ) 
			lines(out[,1],out[,2],col=colz[i],type="h")
			
  		}
	}
	
  	 eksons <- .munion(gene_id)
	
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

#splicing indeks

plotSI <- function(chr,start,end,strand,exp1,exp2,db = "FALSE" )
{ 
   
	title <- paste("Splicing Index for (",exp1,exp2, chr,start,end,strand,")")
	outp <- NULL
	
	if (db=='FALSE')
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
      	query <- paste("call seq_indeks(",exp1,",",exp2,",",chr,",",start,",",end,",",strand,",",exp1,");", sep="") 
   		outp <- xmapcore:::.xmc.db.call(query,   xmapcore:::.xmap.internals$con)
   
   		query <- paste("call seq_coverage(",exp1,",",chr,",",start,",",end,",",strand,");", sep="")
  		out1 <- xmapcore:::.xmc.db.call(query,   xmapcore:::.xmap.internals$con)
   
   		query <- paste("call seq_coverage(",exp2,",",chr,",",start,",",end,",",strand,");", sep="") 
	    out2 <- xmapcore:::.xmc.db.call(query,   xmapcore:::.xmap.internals$con)
  		
  		m1 <- max(out1[,2],na.rm=FALSE)
  		m2 <- max(out2[,2],na.rm=FALSE)
  		m <- max(m1,m2)
		  
        print(outp,digits=16)

		plot(0,0,xlab="Position on the chromosome",
  	           ylab="Nr of reads", 
 	           xlim=c(start,end), 
 	           ylim=c(0,m) ,
 	           main=title,
 	           bg="grey100", 
 	           col.axis="tan4" )
		  
		lines(out1[,1],out1[,2],col="blue",type="l")
  		lines(out2[,1],out2[,2],col="red",type="l")
        lines(outp[,1],outp[,2],col="lawngreen",type="l")
 	 
		legend("topleft",
				 legend=c(paste("Sample",exp1),paste("Sample",exp2),"Splicing index"),
				 fill=c("blue","red","lawngreen"))
	  }
}
