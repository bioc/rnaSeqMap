bam2sig <- function(annot, covdesc="covdesc", species="homo_sapiens")
# annot is for now a table ID,chr,  st, end, strand
{

setSpecies(species)
tout <- NULL

numgenes <- dim(annot)[1]
  cvd <- read.table(covdesc)
  numbams <- dim(cvd)[1] 

for (i in 1:numgenes)
  {  
     ch <- rnaSeqMap:::.chr.convert((annot[i,"chr"]));
     st <- as.numeric(annot[i,"start"]);
     en <- as.numeric(annot[i,"end"]);
     str <- as.numeric(annot[i,"strand"]);
     rs <- newSeqReads(ch,st, en, str);
     rs <- getBamData(rs,1:numbams)
     counts <- lapply(rs@data, .countz) 
     tout <- rbind(tout, counts)
  } 
rownames(tout) <- annot[,"name"]
tout
}

.countz <- function(x)
{
  out <- dim(x)[1]
  if (out==1 & x[1,2]==0) out=0  # case of (0,0) read
  out
}


