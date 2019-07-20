#' Bootstrap samples
#'
#' @param query_db Query database
#' @param results_db Results database
#' @param bootstraps Number of bootstraps to perform
#' @param eps Value to add to denominator to prevent division by zero
#'
#' @return Bootstrapped scores
#' @export
bootstrap_samples <- function(query_db,results_db,bootstraps=10000,eps=1e-6) {
  query_rows <- sample(x=1:nrow(query_db),size=bootstraps,replace=TRUE)
  results_rows <- sample(x=1:nrow(results_db),size=bootstraps,replace=TRUE)
  X <- query_db[query_rows,]
  Y <- results_db[results_rows,]

  # prevent division by zero by assigning at least eps of value
  t1  <- pmax(sqrt(apply(X, 1, crossprod)),eps)
  t2 <- pmax(sqrt(apply(Y, 1, crossprod)),eps)

  X <- X/t1
  Y <- Y/t2

  bootstrap_scores <- rowSums(X * Y)
  return(bootstrap_scores)
}


#' Bind Semantic Types
#'
#' @param embedding_df The embedding dataframe
#'
#' @return Embedding dataframe with CUI, SemanticType, and String columns prepended
#' @export
bind_semantic_types <- function(embedding_df){

  semantic_embedding_df <- data.frame(embedding_df)


  if(all(c("CUI", "SemanticType", "String") %in% colnames(semantic_embedding_df))){
    return(dplyr::select(.data$CUI,.data$SemanticType,.data$String,dplyr::everything()))
  }

  if(!("CUI" %in% colnames(semantic_embedding_df))){
    semantic_embedding_df <- semantic_embedding_df %>%
      dplyr::mutate(CUI=rownames(semantic_embedding_df)) %>%
      dplyr::select(.data$CUI, dplyr::everything())
  }

  data(semantic_type)

  semantic_embedding_df <- semantic_embedding_df %>%
    dplyr::inner_join(semantic_type,by='CUI') %>%
    dplyr::select(.data$CUI,.data$SemanticType,.data$String,dplyr::everything())
  return(semantic_embedding_df)
}


#' Load an embedding from a csv file.
#'
#' @param filename File name to be loaded. CUIs should be in first column.
#' @param join_semantic_type Flag to join semantic type
#'
#' @return Dataframe containing the embeddings (and semnatic type)
#' @export
load_embeddings <- function(filename,join_semantic_type=TRUE) {
  embedding_df <- readr::read_csv(filename,col_names=TRUE)
  embedding_df[,-1] <- apply(embedding_df[,-1],2,as.numeric)
  colnames(embedding_df)[1] <- "CUI"
  if(join_semantic_type){
    embedding_df <- bind_semantic_types(embedding_df)
  }
  return(embedding_df)
}

#' Checking embedding semantic type columns
#'
#' @param embedding_df The embedding dataframe
#'
#' @return True/False that the first 3 columns are CUI, SemanticType, String
#' @export
check_embedding_semantic_columns <- function(embedding_df){
  if(is.null(colnames(embedding_df))){
    return(FALSE)
  }
  # Check that an embedding dataframe starts with CUI, SemanticType, String columns
  return(all(c("CUI", "SemanticType", "String")==colnames(embedding_df[1:3])))
}




#' Similarity between embedding vectors and query vector
#'
#' @param word_vectors The embedding vectors
#' @param query The query vector
#' @param sort_result Flag whether to sort by cosine similarity
#'
#' @return List of cosine similarities to query vector for each word in the embedding
#' @export
query_similarity <- function(word_vectors,query,sort_result=TRUE) {
  word_vectors <- as.matrix(word_vectors)
  #word_vectors_norm <- sqrt(rowSums(word_vectors^2))
  query_vec <- word_vectors[query,,drop=FALSE]
  cos_sim <- text2vec::sim2(query_vec,
                            word_vectors,
                            method = 'cosine',
                            norm = 'l2'
  )

  cos_sim <- cos_sim[1,]
  names(cos_sim) <- row.names(word_vectors)
  if(sort_result) {
    cos_sim <- sort(cos_sim,decreasing = TRUE)
  }
  ## Remove query vec ##
  #cos_dist <- cos_dist[which(names(cos_dist) != query)]
  return(cos_sim)
}

#' Cosine similarity between two vectors.
#'
#' @param vec1 First vector
#' @param vec2 Second vector
#'
#' @return The cosine similarity between two vectors
#' @export
cosine_similarity <- function(vec1,vec2){
  sim <- sum(vec1*vec2)/(sqrt(sum(vec1^2))*sqrt(sum(vec2^2)))
  #If vec1 or vec2 is 0, then R normally returns NA because a 0 division error, this catches this
  if(is.na(sim)){return(0)}
  return(sim)
}

