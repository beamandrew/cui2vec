#' Benchmarking comorbidities task
#'
#' The benchmarking strategy leverages previously published ‘known’ relationships between medical concepts.
#' We compare how similar the embeddings for a pair of concepts are by computing the
#' cosine similarity of their corresponding vectors,
#' and use this similarity to assess whether or not the two concepts are related.
#' \code{\link{benchmark_comorbidities}} focuses on an embedding's ability to identify comorbidities.
#' A comorbidity is a disease or condition that frequently accompanies a  primary  diagnosis.
#'
#' @param embedding_df The embedding dataframe with bound semantic types
#' @param sig_level The significance level threshold
#' @param bootstraps The number of bootstraps to perform
#'
#' @return Dataframe of performance on this benchmark per disease
#' @export
benchmark_comorbidities <- function(embedding_df,sig_level=0.05, bootstraps=10000) {
  if(!check_embedding_semantic_columns(embedding_df)){
    print("Binding semantic type information...")
    embedding_df <- bind_semantic_types(embedding_df)
  }
  #Create the data frame that we will append all scores to
  all_files <- list.files(system.file("extdata", "benchmarks", "comorbidities", package = "cui2vec"))
  df <- NULL
  #Looping over all the comorbidities
  for(file in all_files){
    #Load the file
    # TODO: change path
    comorbidity <- utils::read.delim(system.file("extdata", "benchmarks", "comorbidities", file, package = "cui2vec"),
                                     stringsAsFactors = F)

    #Get only the concepts and associations that are in the embedding
    concepts <- intersect(comorbidity$CUI[which(comorbidity$Type=='Concept')],embedding_df$CUI)
    if(length(concepts)==0){
      stop("No concepts in the intersection of comorbidity benchmark and embedding")
    }
    strings <- comorbidity %>%  dplyr::filter(.data$CUI %in% concepts) %>% dplyr::select(.data$String)
    associations <- intersect(comorbidity$CUI[which(comorbidity$Type=='Association')],embedding_df$CUI)
    association_semantic_type <- embedding_df %>%
                                    dplyr::filter(.data$CUI %in% associations) %>%
                                    dplyr::select(.data$SemanticType) %>%
                                    unique()
    association_semantic_type <- as.character(association_semantic_type$SemanticType)
    local_df <- NULL
    for(i in 1:length(concepts)) {
      current_concept <- concepts[i]
      concept_semantic_type <- embedding_df %>%
                                dplyr::filter(.data$CUI %in% current_concept) %>%
                                dplyr::select(.data$SemanticType) %>%
                                unique() %>%
                                as.character()

      query_db <-  embedding_df %>% dplyr::filter(.data$SemanticType %in% concept_semantic_type) %>% data.frame()
      result_db <- embedding_df %>% dplyr::filter(.data$SemanticType %in% association_semantic_type | .data$CUI == current_concept) %>% data.frame()

      query_word_vectors <- as.matrix(query_db[,c(-1,-2,-3)])
      row.names(query_word_vectors) <- query_db$CUI

      word_vectors <- as.matrix(result_db[,c(-1,-2,-3)])
      row.names(word_vectors) <- result_db$CUI

      sim_scores <- query_similarity(word_vectors,current_concept)
      sim_scores <- sim_scores[-1] #remove self

      null_scores <- bootstrap_samples(query_word_vectors,word_vectors, bootstraps=bootstraps)
      observed_score <- sim_scores[associations]
      sig_threshold <- stats::quantile(null_scores,p=1-sig_level)

      num_positive <- length(which(observed_score > sig_threshold))
      power <- num_positive/length(associations)

      local_df <- rbind(local_df,data.frame(Test="Comorbidity",File=file,Concept=strings$String[i],Num_Positive=num_positive,Total=length(associations),Score=power,Best=FALSE))
    }
    best_score <- max(local_df$Score)
    local_df[which(local_df$Score == best_score),"Best"] <- TRUE
    df <- rbind(df,local_df)
  }
  return(df)
}

#' Benchmarking causative task
#'
#' The benchmarking strategy leverages previously published ‘known’ relationships between medical concepts.
#' We compare how similar the embeddings for a pair of concepts are by computing the
#' cosine similarity of their corresponding vectors,
#' and use this similarity to assess whether or not the two concepts are related.
#' \code{\link{benchmark_causative}} assesses an embedding's ability to recover causes from
#' the UMLS' table (MRREL) of entities known to be the cause of a certain result.
#'
#' @param embedding_df The embedding data frame with bound semantic type
#' @param sig_level The significance level threshold
#' @param bootstraps The number of boostraps to perform
#' @param verbose Flag for verbosity
#'
#' @return Dataframe of scores on this task per causative relationship
#' @export
benchmark_causative <- function(embedding_df,sig_level = 0.05, bootstraps=10000, verbose = TRUE){
  if(!check_embedding_semantic_columns(embedding_df)){
    print("Binding semantic type information...")
    embedding_df <- bind_semantic_types(embedding_df)
  }
  df <- NULL
  #Looping over all the causative files
  all_files <- list.files(system.file("extdata", "benchmarks", "causative", package = "cui2vec"))
  for(file in all_files){
    print(paste0("Benchmarking using ",file,"..."))
    # TODO: change path
    causative <- utils::read.delim(system.file("extdata", "benchmarks", "causative", file, package = "cui2vec"),
                                   stringsAsFactors = F)
    #If there are no reference CUIs we simply take all the CUIs in the embedding
    #Need to get only the cause:result pairs where both cause AND result are in the embedding
    cuis_in_both <- NULL
    for(j in 1:nrow(causative)) {
      if(causative$CUI_Cause[j] %in% embedding_df$CUI & causative$CUI_Result[j] %in% embedding_df$CUI) {
        cuis_in_both <- c(cuis_in_both,j)
      }
    }
    cuis <- causative[cuis_in_both,]
    if(nrow(cuis) > 0) {
      pb <- utils::txtProgressBar(min=1,max=nrow(cuis),style=3)
      local_df <- NULL
      for(i in 1:nrow(cuis)) {
        current_concept = cuis$CUI_Cause[i]
        current_result = cuis$CUI_Result[i]

        cause_semantic_type <- embedding_df %>% dplyr::filter(.data$CUI %in% current_concept) %>% dplyr::select(.data$SemanticType) %>% unique()
        cause_semantic_type <- as.character(cause_semantic_type$SemanticType)

        result_semantic_type <- embedding_df %>% dplyr::filter(.data$CUI %in% current_result) %>% dplyr::select(.data$SemanticType) %>% unique()
        result_semantic_type <- as.character(result_semantic_type$SemanticType)

        cause_db <-  embedding_df %>% dplyr::filter(.data$SemanticType %in% cause_semantic_type) %>% data.frame()
        result_db <- embedding_df %>% dplyr::filter(.data$SemanticType %in% result_semantic_type | .data$CUI %in% cuis$CUI_Cause) %>% data.frame()

        cause_word_vectors <- as.matrix(cause_db[,c(-1,-2,-3)])
        row.names(cause_word_vectors) <- cause_db$CUI

        word_vectors <- as.matrix(result_db[,c(-1,-2,-3)])
        row.names(word_vectors) <- result_db$CUI

        observed_score <- cosine_similarity(cause_word_vectors[current_concept,],word_vectors[current_result,])

        null_scores <- bootstrap_samples(cause_word_vectors,word_vectors,bootstraps = bootstraps)
        sig_threshold <- stats::quantile(null_scores,p=1-sig_level)

        is_positive <- as.numeric(observed_score > sig_threshold)

        local_df <- rbind(local_df,data.frame(Test="Comorbidity",File=file,Called_Sig=is_positive))
        utils::setTxtProgressBar(pb,i)
      }
      df <- rbind(df,data.frame(Test="Causative",File=file,Num_Positive=sum(local_df$Called_Sig),Total=nrow(local_df),Score=sum(local_df$Called_Sig)/nrow(local_df)))
    }
  }
  return(df)
}

#' Benchmarking NDF RT task
#'
#' The benchmarking strategy leverages previously published ‘known’ relationships between medical concepts.
#' We compare how similar the embeddings for a pair of concepts are by computing the
#' cosine similarity of their corresponding vectors,
#' and use this similarity to assess whether or not the two concepts are related.
#' \code{\link{benchmark_ndf_rt}} assesses an embedding's ability to power to detect "may treat" and "may prevent"
#' relationships using bootstrap scores of random drug-disease pairs.
#'
#' @param embedding_df The embedding data frame with bound semantic type
#' @param sig_level The significance level threshold
#' @param bootstraps The number of bootstraps to perform
#'
#' @return Dataframe of performance on this task
#' @export
benchmark_ndf_rt <- function(embedding_df, sig_level = 0.05, bootstraps=10000){
  if(!check_embedding_semantic_columns(embedding_df)){
    print("Binding semantic type information...")
    embedding_df <- bind_semantic_types(embedding_df)
  }
  #Generate the data frame we will return
  df <- NULL
  #Loop over all the NDF RT files
  all_files <- list.files(system.file("extdata", "benchmarks", "ndf_rt", package = "cui2vec"))
  for(file in all_files){
    # TODO: change path
    ndfrt <- utils::read.delim(system.file("extdata", "benchmarks", "ndf_rt", file, package = "cui2vec"),
                               stringsAsFactors = F)
    #Need to get a list of all treatment cuis with at least one cui in the condition column. Initializing this list as null
    cuis <- NULL
    for (i in 1:length(ndfrt$Treatment)){
      #If the treatment cui is in the embedding and has at least one condition cui, then we include this treatment cui in the list of valid cuis
      if((ndfrt$Treatment[i] %in% embedding_df$CUI) & length(intersect(strsplit(as.character(ndfrt$Condition[i]),','),embedding_df$CUI))>0){
        #Appending the valid cui to the list
        cuis <- c(cuis,i)
      }
    }
    if(length(cuis) > 0) {
      relevent_ndfrt <- ndfrt %>% dplyr::filter(.data$Treatment %in% ndfrt$Treatment[cuis])
      pb <- utils::txtProgressBar(min=1,max=nrow(relevent_ndfrt),style=3)
      local_df <- NULL
      for(i in 1:nrow(relevent_ndfrt)){
        current_treatment = relevent_ndfrt$Treatment[i]
        current_condition = relevent_ndfrt$Condition[i]
        current_condition <- stringr::str_split(stringr::str_trim(current_condition),',')[[1]]
        current_condition <- current_condition[which(current_condition != "")]

        treatment_semantic_type <- embedding_df %>% dplyr::filter(.data$CUI %in% current_treatment) %>% dplyr::select(.data$SemanticType) %>% unique()
        treatment_semantic_type <- as.character(treatment_semantic_type$SemanticType)

        condition_semantic_type <- embedding_df %>% dplyr::filter(.data$CUI %in% current_condition) %>% dplyr::select(.data$SemanticType) %>% unique()
        condition_semantic_type <- as.character(condition_semantic_type$SemanticType)

        treatment_db <-  embedding_df %>% dplyr::filter(.data$SemanticType %in% treatment_semantic_type) %>% data.frame()
        condition_db <- embedding_df %>% dplyr::filter(.data$SemanticType %in% condition_semantic_type) %>% data.frame()

        treatment_word_vectors <- as.matrix(treatment_db[,c(-1,-2,-3)])
        row.names(treatment_word_vectors) <- treatment_db$CUI

        condition_word_vectors <- as.matrix(condition_db[,c(-1,-2,-3)])
        row.names(condition_word_vectors) <- condition_db$CUI

        observed_score <- NULL
        if(length(current_condition) == 1) {
          observed_score <- cosine_similarity(treatment_word_vectors[current_treatment,],condition_word_vectors[current_condition,])
        } else {
          for(j in 1:length(current_condition)) {
            if(current_condition[j] %in% row.names(condition_word_vectors)) {
              observed_score <- c(observed_score,cosine_similarity(treatment_word_vectors[current_treatment,],condition_word_vectors[current_condition[j],]))
            }
          }
        }

        null_scores <- bootstrap_samples(treatment_word_vectors,condition_word_vectors,bootstraps=bootstraps)
        sig_threshold <- stats::quantile(null_scores,p=1-sig_level)

        is_positive <- as.numeric(observed_score > sig_threshold)
        local_df <- rbind(local_df,data.frame(Test="Comorbidity",File=file,Called_Sig=is_positive))
        utils::setTxtProgressBar(pb,i)
      }
      df <- rbind(df,data.frame(Test="NDFRT",File=file,Num_Positive=sum(local_df$Called_Sig),Total=nrow(local_df),Score=sum(local_df$Called_Sig)/nrow(local_df)))
    }
  }
  return(df)
}

#' Benchmarking semantic type task
#'
#' The benchmarking strategy leverages previously published ‘known’ relationships between medical concepts.
#' We compare how similar the embeddings for a pair of concepts are by computing the
#' cosine similarity of their corresponding vectors,
#' and use this similarity to assess whether or not the two concepts are related.
#' \code{\link{benchmark_semantic_type}} assesses an ability to identify semantic types. Semantic types are
#' meta-information about which category a concept belongs to, and these categories are arranged in a hierarchy.
#'
#' @param embedding_df The embedding dataframe with bound semantic types
#' @param sig_level The significance level threshold
#' @param bootstraps The number of bootstraps to perform
#'
#' @return Dataframe of performance on this task per semantic type
#' @export
benchmark_semantic_type <- function(embedding_df, sig_level = 0.05, bootstraps=10000){
  if(!check_embedding_semantic_columns(embedding_df)){
    print("Binding semantic type information...")
    embedding_df <- bind_semantic_types(embedding_df)
  }
  df <- NULL
  pb <- utils::txtProgressBar(min=1,max=length(unique(embedding_df$SemanticType)),style=3)
  counter = 1
  for(semantic_type in unique(embedding_df$SemanticType)) {
    utils::setTxtProgressBar(pb,counter)
    semantic_db <- embedding_df %>% dplyr::filter(.data$SemanticType == semantic_type)
    not_semantic_db <-  embedding_df %>% dplyr::filter(.data$SemanticType != semantic_type)

    semantic_word_vectors <- as.matrix(semantic_db[,c(-1,-2,-3)])
    row.names(semantic_word_vectors) <- semantic_db$CUI

    not_semantic_word_vectors <- as.matrix(not_semantic_db[,c(-1,-2,-3)])
    row.names(not_semantic_word_vectors) <- not_semantic_db$CUI

    sim_scores <- text2vec::sim2(semantic_word_vectors,
                       method = 'cosine',
                       norm = 'l2')
    observed_score <- sim_scores[upper.tri(sim_scores)]

    null_scores <- bootstrap_samples(semantic_word_vectors,not_semantic_word_vectors,bootstraps = bootstraps)
    sig_threshold <- stats::quantile(null_scores,p=1-sig_level)

    num_positive <- length(which(observed_score > sig_threshold))
    power <- num_positive/length(observed_score)
    df <- rbind(df, data.frame(Test="Semantic Type",File=semantic_type,Num_Positive = num_positive, Total = length(observed_score), Score = power))
    counter <- counter + 1
  }
  return(df)
}



#' Benchmarking similarity task
#'
#' The benchmarking strategy leverages previously published ‘known’ relationships between medical concepts.
#' We compare how similar the embeddings for a pair of concepts are by computing the
#' cosine similarity of their corresponding vectors,
#' and use this similarity to assess whether or not the two concepts are related.
#' \code{link{benchmark_similarity}} reports the spearman correlation between the human assessment scores
#' and cosine similarity from the embeddings.
#'
#' @param embedding_df The embedding dataframe
#'
#' @return Dataframe of relatedness
#' @export
benchmark_similarity <- function(embedding_df){
  if(!check_embedding_semantic_columns(embedding_df)){
    print("Binding semantic type information...")
    embedding_df <- bind_semantic_types(embedding_df)
  }
  results = data.frame(Similarity = 0, Relatedness = 0)
  #Take every CUI in the embedding if no reference specified
  #Load the file with concept pairs and mean resident scores
  # TODO: change path

  similar <- utils::read.csv(system.file("extdata", "benchmarks", "similarity_and_relatedness", "UMNSRS_similarity.csv", package = "cui2vec"), header=TRUE,
                             stringsAsFactors = F)
  #Generate a data frame to store the cosine similarity and mean resident similarity for concept pairs
  df <- NULL
  for(i in 1:length(similar$CUI1)){
    #dplyr::selecting only concept pairs where both CUIs are in the embedding
    if((similar$CUI1[i] %in% embedding_df$CUI) & (similar$CUI2[i] %in% embedding_df$CUI)){
      cui_1 <- embedding_df %>% dplyr::filter(.data$CUI == similar$CUI1[i]) %>% data.frame()
      cui_1 <- cui_1[1,c(-1,-2,-3)]
      cui_2 <- embedding_df %>% dplyr::filter(.data$CUI == similar$CUI2[i]) %>% data.frame()
      cui_2 <- cui_2[1,c(-1,-2,-3)]
      df <- rbind(df,data.frame(Similarity=similar$Mean[i],Cosine=cosine_similarity(cui_1,cui_2)))
    }
  }
  r <- stats::cor(df$Similarity,df$Cosine,method="spearman")
  results = data.frame(Test="Human Assessment",File="Similarity",Num_Positive=NA,Total = nrow(df),Score=r)
  relatedness <- utils::read.csv(system.file("extdata", "benchmarks", "similarity_and_relatedness", "UMNSRS_relatedness.csv", package = "cui2vec"), header=TRUE,
                                 stringsAsFactors = F)
  #Generate a data frame to store the cosine similarity and mean resident similarity for concept pairs
  df <- NULL
  for(i in 1:length(relatedness$CUI1)){
    #dplyr::selecting only concept pairs where both CUIs are in the embedding
    if((relatedness$CUI1[i] %in% embedding_df$CUI) & (relatedness$CUI2[i] %in% embedding_df$CUI)){
      #Write to the data frame
      cui_1 <- embedding_df %>% dplyr::filter(.data$CUI == relatedness$CUI1[i]) %>% data.frame()
      cui_1 <- cui_1[1,c(-1,-2,-3)]
      cui_2 <- embedding_df %>% dplyr::filter(.data$CUI == relatedness$CUI2[i]) %>% data.frame()
      cui_2 <- cui_2[1,c(-1,-2,-3)]
      df <- rbind(df,data.frame(Relatedness=relatedness$Mean[i],Cosine=cosine_similarity(cui_1,cui_2)))
    }
  }

  if(is.null(df)){
    stop("No concept pairs where both CUIs are in the embedding")
  }

  r <- stats::cor(df$Relatedness,df$Cosine,method="spearman")
  results = rbind(results,data.frame(Test="Human Assessment",File="Relatedness",Num_Positive=NA,Total = nrow(df),Score=r))
  return(results)
}

