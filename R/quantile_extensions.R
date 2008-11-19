

normalize.quantiles.determine.target <- function(x,target.length=NULL,subset=NULL){

  if (!is.matrix(x)){
    stop("This function expects supplied argument to be matrix")
  }
  if (!is.numeric(x)){
    stop("Supplied argument should be a numeric matrix")
  }
  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (is.integer(x)){
    x <- matrix(as.double(x), rows, cols)
  }
  
  if (is.null(target.length)){
    target.length <- rows
  }
  
  if (target.length <= 0){
    stop("Need positive length for target.length")
  }

  if (is.null(subset)){
    return(.Call("R_qnorm_determine_target",x,target.length,PACKAGE="preprocessCore"))
  } else {
    if (length(subset) != rows){
       stop("subset should have same length as nrows(x)")
    }
    subset <- as.integer(subset)
    return(.Call("R_qnorm_determine_target_via_subset",x, subset,target.length,PACKAGE="preprocessCore"))			
  }

}



normalize.quantiles.use.target <- function(x,target,copy=TRUE,subset=NULL){

  if (!is.matrix(x)){
    stop("This function expects supplied argument to be matrix")
  }
  if (!is.numeric(x)){
    stop("Supplied argument should be a numeric matrix")
  }
  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (is.integer(x)){
    x <- matrix(as.double(x), rows, cols)
  }
  
  if (!is.vector(target)){
     stop("This function expects target to be vector")
  }
  if (!is.numeric(target)){
    stop("Supplied target argument should be a numeric vector")
  }

  if (is.integer(target)){
     target <- as.double(target)
  }
  if (is.null(subset)){	
     return(.Call("R_qnorm_using_target",x,target,copy,PACKAGE="preprocessCore"))
  } else {
    if (length(subset) != rows){
       stop("subset should have same length as nrows(x)")
    }
    subset <- as.integer(subset)
    return(.Call("R_qnorm_using_target_via_subset",x, subset, target, copy, PACKAGE="preprocessCore"))			
  }


}



normalize.quantiles.in.blocks <- function(x,blocks,copy=TRUE){

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (rows != length(blocks)){
    stop("blocks is not vector of correct length")
  }

  if (is.factor(blocks)){
    blocks <- as.integer(blocks)
  }

  if (!is.numeric(blocks)){
    stop("non-numeric vector used for blocks")
  }


  return(.Call("R_qnorm_within_blocks",x,blocks,copy,PACKAGE="preprocessCore"))



}




