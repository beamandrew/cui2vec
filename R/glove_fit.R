#' Constructs GloVe embeddings from co-occurrence matrix
#'
#' @param coccur The co-occurrence matrix
#' @param dim_size The embedding dimension
#' @param iters The number of iterations
#'
#' @return GloVe embedding
#' @export
construct_glove_embedding <- function(coccur,dim_size=100,iters=25) {
  glove = text2vec::GlobalVectors$new(word_vectors_size = dim_size, vocabulary = row.names(coccur), x_max = 10, learning_rate = 0.01)
  word_vectors <- glove$fit_transform(coccur, n_iter = iters, convergence_tol = 1e-4)
  return(word_vectors)
}
