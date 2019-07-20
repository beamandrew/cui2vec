#' Constructs PCA embedding based raw co-occurrence counts
#'
#' @param coccur The co-occurrence matrix
#' @param dim_size The embedding dimension
#' @param iters The number of iterations
#' @param remove_empty Flag to remove vectors that are identically 0
#'
#' @return PCA embedding
#' @export
construct_pca_embedding <- function(coccur,dim_size=100,iters=25,remove_empty=TRUE) {
  fit <- irlba::irlba(coccur,nv=dim_size,maxit=iters,verbose=TRUE)
  W <- fit$u
  vecs <- W
  rownames(vecs) <- rownames(coccur)
  if(remove_empty) {
    ## Remove empty word vectors ##
    vecs <- vecs[which(rowSums(vecs) != 0),]
  }
  return(vecs)
}
