#' Run paired comparisons
#'
#' @param embedding1_df The first embedding's dataframe with bound semantic type
#' @param embedding2_df The second embedding's dataframe with bound semantic type
#' @param benchmark Which benchmarks to perform. Valid inputs are c("all", "comorbidities", "causative", "ndfrt", "semantic_type", "similarity")
#'
#' @return A list with two entries of result dataframes, one for each embedding.
#' @export
compare_embeddings <- function(embedding1_df, embedding2_df,
                               benchmark=c("all", "comorbidities", "causitive", "ndfrt", "semantic_type", "similarity")){
  concepts_in_common <- intersect(embedding1_df$CUI, embedding2_df$CUI)


  if(benchmark=="all"){
    results_1 <- run_all_benchmarks(embedding1_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
    results_2 <- run_all_benchmarks(embedding2_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
  }
  else if(benchmark=="comorbidities"){
    results_1 <- benchmark_comorbidities(embedding1_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
    results_2 <- benchmark_comorbidities(embedding2_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
  }
  else if(benchmark=="causative"){
    results_1 <- benchmark_causative(embedding1_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
    results_2 <- benchmark_causative(embedding2_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
  }
  else if(benchmark=="ndfrt"){
    results_1 <- benchmark_ndf_rt(embedding1_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
    results_2 <- benchmark_ndf_rt(embedding2_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
  }
  else if(benchmark=="semantic_type"){
    results_1 <- benchmark_semantic_type(embedding1_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
    results_2 <- benchmark_semantic_type(embedding2_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
  }
  else if(benchmark=="similarity"){
    results_1 <- benchmark_similarity(embedding1_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
    results_2 <- benchmark_similarity(embedding2_df %>% dplyr::filter(.data$CUI %in% concepts_in_common))
  }
  else{
    stop("Invalid benchmark")
  }

  return(list(embedding_1_results=results_1, embedding_2_results=results_2))
}

