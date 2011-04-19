#plots for nd with list of Rle
#AL


####################################################################################################################

distrSIPlot <- function(nd, ex1, ex2, mi,minsup=5)
{  
	par(mfcol=c(3,1))
	par(mar=c(2,2,1,1))
	
	nucleotides <- nd@start:nd@end 
	m1 <- max(RleList2matrix(nd@data[ex1]),ns.rm = FALSE)
	m2 <- max(RleList2matrix(nd@data[ex2]),ns.rm = FALSE)
	m  <- max(m1,m2)
	
	plot(0,0,xlab="Position on the chromosome",
		 ylab="Nr of reads", 
	     xlim=c(nd@start,nd@end), 
		 ylim=c(0.0,m) ,
		 bg="grey100", 
		 col.axis="tan4" )
	
	lines(nucleotides, RleList2matrix(nd@data[ex1]), col="blue",type="l")
	lines(nucleotides, RleList2matrix(nd@data[ex2]), col="red",type="l")
	
	s1<-sum(RleList2matrix(nd@data[ex1]))
	s2<-sum(RleList2matrix(nd@data[ex2]))
	c=((s1/s2))
	
	outp <- .Call("splicingind", as.numeric(RleList2matrix(nd@data[ex1])), as.numeric(RleList2matrix(nd@data[ex2])),c)
	out <- cbind(nucleotides,(outp))

	plot(0,0,xlab="Position on the chromosome",
		 ylab="Nr of reads", 
	     xlim=c(nd@start,nd@end), 
		 ylim=c(-1,1) ,
		 bg="grey100", 
		 col.axis="tan4" )
	
	lines(nucleotides,out[,2],col="lawngreen",type="l")
	
	nd.reg <- findRegionsAsND(nd,as.integer(mi),minsup=minsup)

    plot(0,0,xlab="Position on the chromosome",
		 ylab="Nr of reads", 
	     xlim=c(nd@start,nd@end), 
		 ylim=c(0.0,m) ,
		 bg="grey100", 
		 col.axis="tan4" )
	
	lines(nucleotides,as.vector(nd.reg@data[[ex1]]),col="blue",type="h")
	lines(nucleotides,as.vector(nd.reg@data[[ex2]]),col="red",type="h")
	

	legend("topleft",
		   legend=c(paste("Sample",ex1),paste("Sample",ex2),"Splicing index",paste("Region ",ex1,"(",mi,",",minsup,")"),paste("Region ",ex2,"(",mi,",",minsup,")")),
		   fill=c("blue","red","lawngreen","blue","red"))
}

#####################################################################################################################
distrCOVPlot <- function(nd,exps)
{  
 if (xmapConnected())  
	{
	colz <- topo.colors(length(exps))
	genes <- geneInChromosome(nd@chr,nd@start,nd@end,nd@strand)
	
	l<-length(exps)
	par(mfcol=c(l+1,1))
	par(mar=c(2,2,1,1))
	m <- NULL
		
	nucleotides <- nd@start:nd@end 
		
		for (k in 1:length(exps)){
		m <- c(m,max(as.vector(nd@data[[k]])))}
		m<-max(m)
		
	for (ii in 1:length(exps))
	{

		plot(0,0,xlab="Position on the chromosome",
			 ylab="Nr of reads", 
			 xlim=c(nd@start,nd@end), 
			 ylim=c(0,m),
			 col.axis="tan4" 
			 )
		
		lines(nucleotides, RleList2matrix(nd@data[ii]), col=colz[ii],type="l")
		legend("topleft", legend=c(paste("sample",exps[ii])), fill=(c(col=colz[ii])))
	}
	
	if (is.null(genes))
	{
		plot(0,0, xlim=c(nd@start,nd@end), ylim=c(-0.5,1) ,bg="grey100", col.axis="tan4")
		
		legend("center", legend=c(paste("No gene in this region !","\n")))
	}
	else 
	{
		plot(0,0,xlab="Position on the chromosome",	ylab="Nr of reads", xlim=c(nd@start,nd@end), ylim=c(-0.5,dim(genes)[1]) ,bg="grey100", col.axis="tan4")
		
		k <- 0
		for (i in 1:dim(genes)[1]) 
		{
			gd <- gene.details(genes[i,1])
			chr    <- gd$space
			start  <- gd$ranges@start
			end    <- gd$ranges@start+gd$ranges@width
			strand <- gd$strand
			
			transcripts <- gene.to.transcript(genes[i,1])
			
			for (f in 1:length(transcripts)) 
			{
				ed <- exon.details(transcript.to.exon(transcripts[f]))
				estart <- ed$ranges@start
				ewidth <- ed$ranges@width
				
				eksons <- rnaSeqMap:::.tunion(transcripts[f])
				
				for (j in 1:length(eksons))
				{
					rect(eksons@start[j],k,eksons@start[j]+eksons@width[j],k+0.2,col=colz[i])
				}
				k<-k+0.3
			}
		}
	}
    legend("bottom", legend=paste("Coverage for chromosome:",nd@chr," strand:",nd@strand))
    legend("bottomright", legend=c(paste("Gene: ",genes[,2])))
	
	}
else
{
	colz <- topo.colors(length(exps))
	l<-length(exps)
	par(mfcol=c(l,1))
	par(mar=c(2,2,1,1))
	nucleotides <- nd@start:nd@end 
	m <- NULL
	
	for (k in 1:length(exps)){
	m <- c(m,max(as.vector(nd@data[[k]])))}
	m<-max(m)
	
	for (ii in 1:length(exps))
	{		
		plot(0,0,xlab="Position on the chromosome",
			 ylab="Nr of reads", 
			 xlim=c(nd@start,nd@end), 
			 ylim=c(0,m),
			 col.axis="tan4" 
			 )
		
		lines(nucleotides, RleList2matrix(nd@data[ii]), col=colz[ii],type="l")
		legend("topleft", legend=c(paste("sample",exps[ii])), fill=(c(col=colz[ii])))
	}
}}

###################################################################################################


distrCOVPlotg <- function(gene_id,exps)
{  
	if (xmapConnected())  
	{
	colz <- topo.colors(length(exps))
	
	l<-length(exps)
	par(mfcol=c(l+1,1))
	par(mar=c(2,2,1,1))
	
	gd <- gene.details(gene_id)
	ed <- exon.details(gene.to.exon(gene_id))
	
	chr    <- gd$space
	start  <- gd$ranges@start
	end    <- gd$ranges@start+gd$ranges@width
	strand <- gd$strand
	
	estart <- ed$ranges@start
 	ewidth <- ed$ranges@width
	nucleotides <- start:end 
	
	genes <- geneInChromosome(chr,start,end,strand)
	
	for (ii in 1:length(exps))
	{
		out<-regionCoverage(chr,start,end,strand,exps[ii])
		
		m <- max(out[,2],ns.rm = FALSE) 
		mm <- max(m,na.rm=FALSE)
		
		plot(0,0,xlab="Position on the chromosome",
			 ylab="Nr of reads", 
			 xlim=c(start,end), 
			 ylim=c(0,mm) ,
			 bg="grey100", 
			 col.axis="tan4")
		lines(out[,1],out[,2],col=colz[ii],type="l")
	}
	
	if (is.null(genes))
	{
		plot(0,0, xlim=c(start,end), ylim=c(-0.5,1) ,bg="grey100", col.axis="tan4")
		
		legend("center", legend=c(paste("No gene in this region !","\n")))
	}
	else 
	{
		plot(0,0,xlab="Position on the chromosome",	ylab="Nr of reads", xlim=c(start,end), ylim=c(-0.5,dim(genes)[1]) ,bg="grey100", col.axis="tan4")
		
		k <- 0
		for (i in 1:dim(genes)[1]) 
		{
			gd <- gene.details(genes[i,1])
			chr    <- gd$space
			start  <- gd$ranges@start
			end    <- gd$ranges@start+gd$ranges@width
			strand <- gd$strand
			
			transcripts <- gene.to.transcript(genes[i,1])
			
			for (f in 1:length(transcripts)) 
			{
				ed <- exon.details(transcript.to.exon(transcripts[f]))
				estart <- ed$ranges@start
				ewidth <- ed$ranges@width
				
				eksons <- .tunion(transcripts[f])
				
				for (j in 1:length(eksons))
				{
					rect(eksons@start[j],k,eksons@start[j]+eksons@width[j],k+0.2,col=colz[i])
				}
				k<-k+0.3
			}
		}
	}
    legend("bottom", legend=paste("Coverage for chromosome:",chr," strand:",strand))
    legend("bottomright", legend=c(paste("Gene: ",genes[,2])))
	
}
else
	{
	cat("No connection with db")
	}}

################################
