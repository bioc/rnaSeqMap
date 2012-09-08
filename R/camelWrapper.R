# WRAPPER for the main functions.... 
# MO 9.07.2012, thanx to Mark Robinson for the idea of cleaning it up :) 



allCamelMeasuresForRegion <- function(ch, st, en, str, cvd, sample.idx1, sample.idx2, sums=NULL)
{
	    ch <- as.character(ch);
        st <- as.numeric(st);
        en <- as.numeric(en);
        str <- as.numeric(str)
        idx.both <- c(sample.idx1, sample.idx2) 
        
        rs <- newSeqReads(ch,st, en, str);
        rs <- getBamData(rs,idx.both, cvd=cvd)
        nd <- getCoverageFromRS(rs, idx.both)

        if (!is.null(sums))   nd <- normalizeBySum(nd, sums)
        
        nd1o <- averageND(nd, sample.idx1)
        nd2o <- averageND(nd, sample.idx2)
        dd.nonorm <- cbind(as.vector(nd1o@data[[1]]), as.vector(nd2o@data[[1]]))
        covDensC1 <- sum(nd1o@data[[1]])/(en-st)
        covDensC2 <- sum(nd2o@data[[1]])/(en-st)
         
        nd1 <- min_maxNormalize(nd1o)
        nd2 <- min_maxNormalize(nd2o)
        dd.mm <- cbind(as.vector(nd1@data[[1]]), as.vector(nd2@data[[1]]))

        nd1 <- densityNormalize(nd1o)
        nd2 <- densityNormalize(nd2o)
        dd.d <- cbind(as.vector(nd1@data[[1]]), as.vector(nd2@data[[1]]))

        DA <-diff_area(dd.nonorm)
        DA.mm <-diff_area(dd.mm)
        DA.d <- diff_area(dd.d)
        PP <-pp_plot(dd.nonorm)
        PP.mm <-pp_plot(dd.mm)
        PP.d <-pp_plot(dd.d)
        QQ <-qq_plot(dd.nonorm)
        QQ.mm <-qq_plot(dd.mm)
        QQ.d <-qq_plot(dd.d)
        HD1 <-hump_diff1(dd.nonorm)
        HD1.mm <-hump_diff1(dd.mm)
        HD1.d <-hump_diff1(dd.d)
        HD2 <-hump_diff2(dd.nonorm)
        HD2.mm <-hump_diff2(dd.mm)
        HD2.d <-hump_diff2(dd.d)


        out <- c(DA, DA.mm, DA.d, PP, PP.mm, PP.d, QQ,QQ.mm,QQ.d, HD1,HD1.mm,HD1.d,HD2,HD2.mm,HD2.d,covDensC1, covDensC2)
        out[is.nan(out)] <- 0
        names(out) <- c( "DA","DAmm","DA.d","PP","PP.mm","PP.d","QQ","QQ.mm","QQ.d","HD1","HD1.mm","HD1.d","HD2","HD2.mm","HD2.d","covDensC1", "covDensC2")
        out
}


.getAvgND <- function(ch, st, en, str, cvd, sample.idx1, sample.idx2,sums=NULL)
{
	    out <- list()
	    ch <- as.character(ch);
        st <- as.numeric(st);
        en <- as.numeric(en);
        str <- as.numeric(str)
        
        idx.both <- c(sample.idx1, sample.idx2)
        rs <- newSeqReads(ch,st, en, str);
        rs <- getBamData(rs,idx.both, cvd=cvd)
        nd <- getCoverageFromRS(rs, idx.both)

        if (!is.null(sums))   nd <- normalizeBySum(nd, sums)
        
      nd1o <- averageND(nd, sample.idx1)
      nd2o <- averageND(nd, sample.idx2)
      nd1o@data[[2]] <- nd2o@data[[1]]
      nd1o
}

# clean it up, when protocols cleaned!!! MO
.simplePlot <- function (nd, exps, xlab="genome coordinates", ylab="coverage") 
{
   nucleotides <- nd@start:nd@end
   m1 <- max(max(nd@data[[exps[1]]]), max(round(nd@data[[exps[2]]])))
   m2 <- min(min(nd@data[[exps[1]]]), min(round(nd@data[[exps[2]]])))   
   plot(nucleotides, nd@data[[exps[1]]], xlim = c(nucleotides[1],nucleotides[length(nucleotides)]), xlab = xlab, ylab = ylab ,ylim = c(m2, m1),col = "red", type = "l")
   lines(nucleotides, nd@data[[exps[2]]], col = "blue", type = "l")
   if (length(exps)>2) lines(nucleotides, nd@data[[exps[3]]], col = "green", type = "l")
}

simplePlot <- function (nd, exps, xlab="genome coordinates", ylab="coverage") 
{
   nucleotides <- nd@start:nd@end
   m1 <- max(max(nd@data[[exps[1]]]), max(round(nd@data[[exps[2]]])))
   m2 <- min(min(nd@data[[exps[1]]]), min(round(nd@data[[exps[2]]])))   
   plot(nucleotides, nd@data[[exps[1]]], xlim = c(nucleotides[1],nucleotides[length(nucleotides)]), xlab = xlab, ylab = ylab ,ylim = c(m2, m1),col = "red", type = "l")
   lines(nucleotides, nd@data[[exps[2]]], col = "blue", type = "l")
   if (length(exps)>2) lines(nucleotides, nd@data[[exps[3]]], col = "green", type = "l")
}




gRanges2CamelMeasures <- function(gR, cvd,  sample.idx1, sample.idx2, sums=NULL, progress=NULL)
{
	out <- NULL
	for (i in 1:length(gR))
	{
		ch <- as.character(gR@seqnames[i])
		st <- as.numeric(gR@ranges@start[i])
		en <- st +  as.numeric(gR@ranges@width[i]) -1
		strand <- as.character(gR@strand[i])
		strand[strand=="+"] <- 1
		strand[strand=="-"] <- (-1)
		mm <- allCamelMeasuresForRegion(ch, st, en, strand, cvd,  sample.idx1, sample.idx2, sums)
		out <- rbind(out, mm)
		if (!is.null(progress)) 
          if (i%%progress==0) cat(".") 
	}
	rownames(out) <- names(gR)
	out
}
# tmp <- gRanges2CamelMeasures(my.gR, cvd, 1:3, 4:7, progress=1)

# clean it up, when protocols cleaned!!! MO
.fiveCol2GRanges <- function(t)  
# t is a table with id, chr, start, end, strand
{
  ids <- unlist(t[,1])	
  seqnames <- unlist(t[,2])
  ranges <- IRanges(start=unlist(t[,3]), end=unlist(t[,4]))
  strand <- unlist(t[,5])
  strand[strand==1] <- "+"
  strand[strand==-1] <-"-"
  gr <- GRanges(seqnames, ranges, strand)	
  names(gr) <- ids
  gr
}

fiveCol2GRanges <- function(t)  
# t is a table with id, chr, start, end, strand
{
  ids <- unlist(t[,1])	
  seqnames <- unlist(t[,2])
  ranges <- IRanges(start=unlist(t[,3]), end=unlist(t[,4]))
  strand <- unlist(t[,5])
  strand[strand==1] <- "+"
  strand[strand==-1] <-"-"
  gr <- GRanges(seqnames, ranges, strand)	
  names(gr) <- ids
  gr
}

#my.gR <- .fiveCol2GRanges(region.annot[1:5, ])






