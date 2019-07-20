#' Run all benchmarks
#'
#' @param embedding_df The embedding dataframe with bound semantic types
#'
#' @return Results of all benchmarks
#' @export
run_all_benchmarks <- function(embedding_df) {
  if(!check_embedding_semantic_columns(embedding_df)){
    print("Binding semantic type information...")
    embedding_df <- bind_semantic_types(embedding_df)
  }
  comorbidities <- benchmark_comorbidities(embedding_df)
  comorbidities <- comorbidities %>% dplyr::filter(.data$Best == TRUE) %>% dplyr::select(.data$Test, .data$File, .data$Num_Positive, .data$Total, .data$Score)
  causative <- benchmark_causative(embedding_df)
  ndfrt <- benchmark_ndf_rt(embedding_df)
  semantic <- benchmark_semantic_type(embedding_df)
  human_sim <- benchmark_similarity(embedding_df)
  results <- rbind(comorbidities,causative,ndfrt,semantic,human_sim)
  return(results)
}



