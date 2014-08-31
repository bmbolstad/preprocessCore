


rcModelPLM <- function(y,row.effects=NULL, input.scale=NULL){
  if (!is.matrix(y))
    stop("argument should be matrix")
  PsiCode <- 0
  PsiK <- 1.345
  if (is.null(row.effects)){
    .Call("R_rlm_rma_default_model",y,PsiCode,PsiK,input.scale,PACKAGE="preprocessCore")
  } else {
    if (length(row.effects) != nrow(y)){
       stop("row.effects parameter should be same length as number of rows")
    }  
    if (abs(sum(row.effects)) > length(row.effects)*.Machine$double.eps){
       stop("row.effects should sum to zero")
    }
    .Call("R_rlm_rma_given_probe_effects",y,as.double(row.effects),PsiCode,PsiK,input.scale,PACKAGE="preprocessCore") 
  }	
}



rcModelWPLM <- function(y, w, row.effects=NULL, input.scale=NULL){
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
  if (is.null(row.effects)){
     .Call("R_wrlm_rma_default_model",y,PsiCode,PsiK,as.double(w),input.scale,PACKAGE="preprocessCore")
  } else {
    if (length(row.effects) != nrow(y)){
       stop("row.effects parameter should be same length as number of rows")
    }  
    if (abs(sum(row.effects)) > length(row.effects)*.Machine$double.eps){
       stop("row.effects should sum to zero")
    }
    .Call("R_wrlm_rma_given_probe_effects",y,as.double(row.effects),PsiCode,PsiK,as.double(w),input.scale,PACKAGE="preprocessCore") 
  }	

}



rcModelMedianPolish <- function(y){
  if (!is.matrix(y))
    stop("argument should be matrix")
  PsiCode <- 0
  PsiK <- 1.345
  .Call("R_medianpolish_rma_default_model",y,PACKAGE="preprocessCore")
}




subrcModelMedianPolish <- function(y,group.labels){

  if (!is.matrix(y))
    stop("argument should be matrix")

  if (!is.double(y) & is.numeric(y))
    y <- matrix(as.double(y),dim(y)[1],dim(y)[2])
  else if (!is.numeric(y))
    stop("argument should be numeric matrix")

  rowIndexList <- convert.group.labels(group.labels)
  
  x <- .Call("R_sub_rcModelSummarize_medianpolish", y, rowIndexList,PACKAGE="preprocessCore")

  names(x) <- names(rowIndexList)
  x
}





subrcModelPLM <- function(y,group.labels,row.effects=NULL, input.scale=NULL){

  if (!is.matrix(y))
    stop("argument should be matrix")  

  if (!is.double(y) & is.numeric(y))
    y <- matrix(as.double(y),dim(y)[1],dim(y)[2])
  else if (!is.numeric(y))
    stop("argument should be numeric matrix")
   
  rowIndexList <- convert.group.labels(group.labels)
 
  PsiCode <- 0
  PsiK <- 1.345

  if (is.null(row.effects)){
    x <- .Call("R_sub_rcModelSummarize_plm", y, rowIndexList, PsiCode, PsiK, input.scale,PACKAGE="preprocessCore")
    names(x) <- names(rowIndexList)
    x

  } else {
    stop("row.effects not yet implemented for subrcModelPLM")
    if (length(row.effects) != nrow(y)){
       stop("row.effects parameter should be same length as number of rows")
    }  
    if (abs(sum(row.effects)) > 10*.Machine$double.eps){
       stop("row.effects should sum to zero")
    }
    .Call("R_rlm_rma_given_probe_effects",y,as.double(row.effects),PsiCode,PsiK,input.scale,PACKAGE="preprocessCore") 
  }	
}

