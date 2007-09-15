


rcModelPLM <- function(y){
  if (!is.matrix(y))
    stop("argument should be matrix")
  PsiCode <- 0
  PsiK <- 1.345
  .Call("R_rlm_rma_default_model",y,PsiCode,PsiK,PACKAGE="preprocessCore")
}



rcModelWPLM <- function(y, w){
  if (!is.matrix(y))
    stop("argument should be matrix")
  if (is.vector(w)){
    if (length(w) != prod(dim(y))){
      stop("weights are not correct length")
    }
  } else if (is.matrix(w)){
    if (!all(dim(w) == dim(y))){
      stop("weights should be same dimension as input matrix")
    }

  }
  if (any(w < 0)){
    stop("weights should be no negative")
  }

  
    
  PsiCode <- 0
  PsiK <- 1.345
  .Call("R_wrlm_rma_default_model",y,PsiCode,PsiK,as.double(w),PACKAGE="preprocessCore")

}



rcModelMedianPolish <- function(y){
  if (!is.matrix(y))
    stop("argument should be matrix")
  PsiCode <- 0
  PsiK <- 1.345
  .Call("R_medianpolish_rma_default_model",y,PACKAGE="preprocessCore")
}
