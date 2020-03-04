#' Construct pointwise mutual information matrix
#'
#' @param cooccur The co-occurrence matrix
#' @param singletons The dataframe of CUIs and their counts
#' @param N The number of bins
#' @param smooth The smoothing factor
#'
#' @return A dataframe of pointwise mutual information
#' @export
construct_pmi <- function(cooccur,singletons,N,smooth=0.75) {
  singletons$Count <- singletons$Count^smooth/N^smooth
  concept_list <- row.names(cooccur)
  nz <- Matrix::which(cooccur != 0, arr.ind = TRUE)
  # masking the lower half of the matrix because cooccur will always be symmetric
  # don't want to double count CUI1-CUI2 as CUI2-CUI1
  nz <- nz[which(nz[,1] <= nz[,2]),]
  pmi_df <- data.frame(Concept_1 = concept_list[nz[,1]], Concept_2 = concept_list[nz[,2]], JointProb = cooccur[nz]/N, stringsAsFactors = F)
  pmi_df <- pmi_df %>%
    dplyr::inner_join(singletons,by=c("Concept_1" = "CUI")) %>%
    dplyr::rename(Concept_1_Prob=.data$Count) %>%
    dplyr::inner_join(singletons,by=c("Concept_2" = "CUI")) %>%
    dplyr::rename(Concept_2_Prob=.data$Count) %>%
    dplyr::mutate(PMI = log(.data$JointProb/(.data$Concept_1_Prob  * .data$Concept_2_Prob))) %>%
    dplyr::select(.data$Concept_1, .data$Concept_2, .data$PMI)
  return(pmi_df)
}

#' Construct the shifted positive pointwise mutual information matrix
#'
#' @param pmi Pointwise mutual information matrix from \code{\link{construct_pmi}}
#' @param k Shift in SSPMI formula
#'
#' @return The shifted positive pointwise mutual information matrix
#' @export
construct_sppmi <- function(pmi,k=1) {
  sppmi_df <- pmi %>%
    dplyr::mutate(SPPMI = pmax(.data$PMI - log(k),0)) %>%
    dplyr::select(.data$Concept_1,.data$Concept_2,.data$SPPMI)

  all_words <- unique(c(sppmi_df$Concept_1,sppmi_df$Concept_2))
  word_2_index <- 1:length(all_words)
  names(word_2_index) <- all_words

  i <- as.numeric(word_2_index[as.character(sppmi_df$Concept_1)])
  j <- as.numeric(word_2_index[as.character(sppmi_df$Concept_2)])
  x <- as.numeric(sppmi_df$SPPMI)

  ## Remove 0s ##
  non_zero <- which(x != 0)
  i <- i[non_zero]
  j <- j[non_zero]
  x <- x[non_zero]

  ism <- c(i,j)
  jsm <- c(j,i)
  xsm <- c(x,x)
  sppmi <- Matrix::sparseMatrix(i=ism,j=jsm,x=xsm)
  rownames(sppmi) <- all_words
  colnames(sppmi) <- all_words
  return(sppmi)
}

#' Construct word2vec embedding
#'
#' @param sppmi SPPMI matrix from \code{\link{construct_sppmi}}
#' @param dim_size The embedding dimension
#' @param iters The number of iterations
#' @param remove_empty Flag to remove vectors that are identically 0
#' @param use_sum Flag to add \eqn{V_d\sqrt{\Sigma_d}} to embedding calculation
#'
#' @return word2vec
#' @export
construct_word2vec_embedding <- function(sppmi,dim_size=100,iters=25,remove_empty=TRUE,use_sum=TRUE) {
  fit <- irlba::irlba(sppmi,nv=dim_size,maxit=iters,verbose=TRUE)
  W <- fit$u %*% diag(sqrt(fit$d))
  vecs <- W
  if(use_sum) {
    C <- fit$v %*% diag(sqrt(fit$d))
    vecs <- vecs + C
  }
  rownames(vecs) <- rownames(sppmi)
  if(remove_empty) {
    ## Remove empty word vectors ##
    vecs <- vecs[which(rowSums(vecs) != 0),]
  }
  return(vecs)
}

