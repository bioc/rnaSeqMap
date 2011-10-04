generatorAddSquare <- function(nd, deg, length.prop=0.5)
# generator of additive degeneration - adds a square
{
    nddegen <- nd
	d1 <- nd@data[[1]]
	mm <- max(d1)
	ll <- length(d1)
	sr <- trunc(ll*(1-length.prop))
	degenerat <- Rle(c(0,mm*deg), c(sr,ll-sr))
	
	nddegen@data[[2]] <- d1 + degenerat 
	#wpisuje na druga pozycje zdegradowane rle
	
	nddegen
}

generatorAdd <- function(nd, deg, length.prop=0.5)
# generator of additive degeneration - adds another coverage profile
{
    nddegen <- nd
	d1 <- nd@data[[1]]
	mm <- max(d1)
	ll <- length(d1)
	sr <- trunc(ll*(1-length.prop))
	degenerat <- Rle(c(0,mm*deg), c(sr,ll-sr))
	
	nddegen@data[[2]] <- d1 + degenerat 
	#wpisuje na druga pozycje zdegradowane rle
	
	nddegen
}


generatorMultiply <- function(nd, deg, length.prop=0.5)
# multiplies a portion of the profile
{
    nddegen <- nd
	d1 <- nd@data[[1]]
	ll <- length(d1)
	sr <- trunc(ll*(1-length.prop))
	v1 <- vector(length=sr)
	v2 <- vector(length=ll-sr)
	v1[1:sr] <- as.vector(d1)[1:sr]
	v2[1:(ll-sr)] <- as.vector(d1)[(sr+1):ll]
	v2 <- v2*(1+deg)
	degenerat <- Rle(c(v1,v2))
	
	nddegen@data[[2]] <- degenerat 
	#wpisuje na druga pozycje zdegradowane rle
	
	nddegen
}

generatorTrunc <- function(nd,deg)
# generator of profile truncation
{
	nddegen <- nd
	d1 <- nd@data[[1]]
	ll <- length(d1)
	sr <- trunc(ll*deg)
	#cat(ll, " ",sr,"\n")
	if (sr>0) d1[1:sr] <- 0

	nddegen@data[[2]] <- d1
	
	nddegen 
}


generatorSynth <- function(nd, deg, length.prop=0.5)
{
	nddegen <- nd
    d1 <- nd@data[[1]]
	mm <- max(d1)
	if (mm==0) mm <- 10
	ll <- length(d1)
	d1 <- round(mm*abs(sin((1:ll)*2*pi/ll)))
	nddegen@data[[1]] <- Rle(d1)

	# now the part from genMultiply
	nddegen <- nd
	sr <- trunc(ll*(1-length.prop))
	v1 <- vector(length=sr)
	v2 <- vector(length=ll-sr)
	v1[1:sr] <- as.vector(d1)[1:sr]
	v2[1:(ll-sr)] <- as.vector(d1)[(sr+1):ll]
	v2 <- v2*(1+deg)
	degenerat <- Rle(c(v1,v2))
	
	nddegen@data[[1]] <- Rle(d1)
	nddegen@data[[2]] <- degenerat 
	nddegen
}


generatorPeak <- function(nd, deg, sr=10, mult=10)
# generator that makes a peak
{
    nddegen <- nd
	d1 <- nd@data[[1]]
	mm <- max(d1)*mult
	ll <- length(d1)
	if (ll >50+sr) degenerat <- Rle(c(0,mm*deg,0), c(sr,50, ll-sr-50))
	else degenerat <- Rle(mm*deg,ll)
	
	nddegen@data[[2]] <- d1 + degenerat 
	nddegen
}

