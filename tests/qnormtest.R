library(preprocessCore)

err.tol <- 10^-8

x <- matrix(c(100,15,200,250,110,16.5,220,275,120,18,240,300),ncol=3)
x
normalize.quantiles(x)

x.norm.truth <- matrix(rep(c(110.0,16.5,220,275.0),3),ncol=3)

if (all(abs(x.norm.truth - normalize.quantiles(x)) < err.tol) != TRUE){
	stop("Disagreement in normalize.quantiles(x)")
}

normalize.quantiles.determine.target(x)

x.norm.target.truth <- c(16.5,110.0,220.0,275.0)

if (all(abs(x.norm.target.truth - normalize.quantiles.determine.target(x)) < err.tol) != TRUE){
	stop("Disagreement in normalize.quantiles.determine.target(x)")
}


y <- x
y[2,2] <- NA
y
normalize.quantiles(y)

y.norm.target.truth <- c(47.6666666666667,134.4444444444444,226.1111111111111,275.0000000000000)

y.norm.truth <- matrix(c(134.4444444444444,  47.6666666666667, 134.4444444444444,
                         47.6666666666667,                NA,  47.6666666666667,
                        226.1111111111111, 180.2777777777778, 226.1111111111111,
                        275.0000000000000, 275.0000000000000, 275.0000000000000),byrow=TRUE,ncol=3)


if (all(abs(y.norm.truth - normalize.quantiles(y)) < err.tol,na.rm=TRUE) != TRUE){
	stop("Disagreement in normalize.quantiles(y)")
}



if (all(abs(y.norm.target.truth - normalize.quantiles.determine.target(y)) < err.tol) != TRUE){
	stop("Disagreement in normalize.quantiles.determine.target(y)")
}



if (all(abs(normalize.quantiles.use.target(y,y.norm.target.truth) - y.norm.truth) < err.tol,na.rm=TRUE) != TRUE){
		stop("Disagreement in normalize.quantiles.use.target(y)")
}


x <- matrix(c(100,15,200,250,110,16.5,220,275,120,18,240,300),ncol=3)
rownames(x) <- letters[1:4]
colnames(x) <- LETTERS[1:3]
y <- normalize.quantiles(x, keep.names = TRUE)
if(!all(colnames(x)==colnames(y))){
    stop("Disagreement between initial and final column names despite keep.names=TRUE")
}
if(!all(rownames(x)==rownames(y))){
    stop("Disagreement between initial and final row names despite keep.names=TRUE")
}
