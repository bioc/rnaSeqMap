#######################################################################
#  Analysis pipelines for RNAseq and not only
#  Anna Lesniewska, Michal Okoniewski,   ostatni update 18.09.2010
#######################################################################
# Building an object for edgeR
buildDESeq <- function(genes,exps,conds=NULL)
{
        if (is.null(conds)) stop ("Vector of conditions not specified")
	tabin <- NULL
	namesr <- genes
	for (i in 1:length(genes))
	{
		out <- .rsCount(genes[i],exps)
		tabin <- rbind(tabin,c(out))
	}
	rownames(tabin) <- namesr
	conds <- factor(conds)

	cds <- newCountDataSet(tabin,conds)
	cds
}

# Building an object for edgeR
buildDGEList <- function(genes,exps,conds=NULL)
{
        if (is.null(conds)) stop ("Vector of conditions not specified")

        tabin <- NULL
        namesr <- genes

        for (i in 1:length(genes))
        {
                out <- .rsCount(genes[i],exps)
                tabin <- rbind(tabin,c(out))
        }
        rownames(tabin) <- namesr
        conds <- factor(conds)

        f <- calcNormFactors(tabin)
        f <- f/exp(mean(log(f)))

        g <- gsub("R[1-2]L[1-8]", "", colnames(tabin))

        d <- DGEList(tabin,conds,lib.size=colSums(tabin)*f)
        d
}



