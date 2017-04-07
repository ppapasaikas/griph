library(Rcpp)
sourceCpp("/work/gbioinfo/papapana/FMI_groups/SingleCell_Analysis/JSDmat.cpp")

CalcJSDivergence <- function(x){

  if( ! is.numeric(x) | ! is.matrix(x)  ){
    stop("x must be a numeric matrix. ")
  }
  
    x <- t(x)
    
    result <- JSDmat(x) # this function only calculates the upper triangle
    
    result <- result + t(result)
    
    colnames(result) <- rownames(x)
    rownames(result) <- rownames(x)
  result
}





sourceCpp("/work/gbioinfo/papapana/FMI_groups/SingleCell_Analysis/HellingerMat.cpp")
CalcHellingerDist <- function(x){

  if( ! is.numeric(x) | ! is.matrix(x)  ){
    stop("x must be a numeric matrix. ")
  }
    
    x <- t(x)
 
    result <- HellingerMat(x) # this function only calculates the upper triangle
    result <- result + t(result)
    
    colnames(result) <- rownames(x)
    rownames(result) <- rownames(x)
  result
}





CalcKSDist <- function(M){
GeneCounts=rowSums(M)
M<-sweep(M,2,colSums(M),FUN="/")
pairs=combn(c(1:ncol(M)),2,simplify=F)
result=diag(nrow=ncol(M))
result[lower.tri(result)]=sapply(pairs,function(x) max(abs(cumsum(M[order(GeneCounts,decreasing = T),x[1]])-cumsum(M[order(GeneCounts,decreasing = T),x[2]])   ))      )
diag(result)=0
result=result+t(result)
}



