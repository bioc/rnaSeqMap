#########################################################
# measures of difference between RNAseq coverage regions
# MO, AL, AS  lipiec 2011
#########################################################

ks_test <- function(dd)
{
	out <- ks.test(dd[,1],dd[,2])$p.value
	out
}


diff_area <- function(dd, cconst = NULL)
{
	if (is.null(cconst)) cconst <- max(abs(dd[,1:2]))
	out <- sum(abs(dd[,1]-dd[,2]))
	n <- dim(dd)[1]
	out <- out/(cconst*n)
	out
}

diff_derivative_area <- function(dd, cconst = NULL)
{
	dd <- .pochodna(dd)
	if (is.null(cconst)) cconst <- max(abs(dd[,1:2]))
	out <- sum(abs(dd[,1]-dd[,2]))
	n <- dim(dd)[1]
	out <- out/(cconst*n)
	out
}

qq_plot <- function(dd)
{
	qqresult<-qqplot(dd[,1], dd[,2], ,plot=F)
	out<-.rmse(qqresult)
	out
}

qq_derivative_plot <- function(dd)
{
	dd <- .pochodna(dd)
	qqDerResult<-qqplot(dd[,1], dd[,2], ,plot=F)
	out <- .rmse(qqDerResult)
	out
}	

pp_plot <- function(dd) {
	dd[,1]<-cumsum(dd[,1])
	dd[,2]<-cumsum(dd[,2])
	out <- sqrt(mean((dd[,1]-dd[,2])^2)/2)
	out
} 

pp_derivative_plot <- function(dd) {
	dd <- .pochodna(dd)
	dd[,1]<-cumsum(dd[,1])
	dd[,2]<-cumsum(dd[,2])
	out <- sqrt(mean((dd[,1]-dd[,2])^2)/2)
	out
} 

hump_diff1 <- function(dd)
{ 
    garby1<-.localMax(dd[,1])
    garby2<-.localMax(dd[,2])
    n1<-length(garby1)
    for (i in 1:n1)
    {
        if (any(garby2==garby1[i])) garby2<-garby2[-which(garby2==garby1[i])]
    }
    out<-0
    for (i in 1:n1)
    {
        out<-out+abs(dd[garby1[i],1]-dd[garby1[i],2])
    }
    n2<-length(garby2)
    if (n2!=0)
    {
        for (i in 1:n2)
        {
            out<-out+abs(dd[garby2[i],1]-dd[garby2[i],2])
        }
    }
    out<-out/(n1+n2)
    if (length(out)==0) out <- 0
	out
}

hump_diff2 <- function(dd)
{ 
    garby1<-.localMax(dd[,1])  
    garby2<-.localMax(dd[,2])
    n<-2*min(length(garby1),length(garby2))
    n1<-length(garby1)
    for (i in 1:n1)
    {
        if (any(garby2==garby1[i])) garby2<-garby2[-which(garby2==garby1[i])]
    }
    n2<-length(garby2)
    out<-0
    for (i in 1:n1)
    {
        out<-out+abs(dd[garby1[i],1]-dd[garby1[i],2])
    }
	
    if (n2!=0)
    {
        for (i in 1:n2)
        {
            out<-out+abs(dd[garby2[i],1]-dd[garby2[i],2])
        }
    }
	if (n==0) n=n1+n2 #in case of no extrema, switches to hump_diff1
    out<-out/n
    if (length(out)==0) out <- 0
	out
}

#########################################################
#utilities
#########################################################
.derywatywka <- function(x){
# calculates pseudo - derivative
	a<-length(x)
	y<-x[2:a]-x[1:(a-1)]
}


.pochodna <- function(dd1)
{
	a <- dim(dd1)[1]
	#cat(dim(dd1))
	dd <- matrix(,a[1]-1,2)
	for (i in 1:2){
		dd[,i]<-.derywatywka(dd1[,i])
	}
	dd[1,1] <- 0
	dd[1,2] <- 0
	pdd	<- dd
    pdd
}

.rmse <- function(obs){
# liczenie bledu srednio-kwadratowego  
	sqrt(mean((obs$x-obs$y)^2)/2)
}    


.localMax <- function(v)
{
#finds local maxima and returns them as index
	vl <- length(v)
	idx <-NULL
	if (v[1]>v[2]) idx <- 1
	for (i in 2:(vl-1))
	{
		if (v[i]>v[i-1] & v[i]>=v[i+1]) idx <- c(idx, i)
	}
	idx
}
