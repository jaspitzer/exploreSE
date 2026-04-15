.check_precomputed_de <- function(SE){
  if(!is.null(SE) && methods::is(SE, "DeeDeeExperiment") && !is.null(DeeDeeExperiment::getDEANames(SE))){
    return(TRUE)
  }else{
    !is.null(S4Vectors::metadata(SE)$de_results) &&
      is.list(S4Vectors::metadata(SE)$de_results) &&
      length(S4Vectors::metadata(SE)$de_results) > 0
  }
}

.check_precomputed_fe <- function(SE){
  if(!is.null(SE) && methods::is(SE, "DeeDeeExperiment") && !is.null(DeeDeeExperiment::getFEANames(SE))){
    return(TRUE)
  }else{
    !is.null(S4Vectors::metadata(SE)$fe_results) &&
      is.list(S4Vectors::metadata(SE)$fe_results) &&
      length(S4Vectors::metadata(SE)$fe_results) > 0
  }
}

.de_results_names <- function(SE){
  if(methods::is(SE, "DeeDeeExperiment")){
    DeeDeeExperiment::getDEANames(SE)
  }else{
    names(S4Vectors::metadata(SE)$de_results)
  }
}

.fe_results_names <- function(SE, NAME){
  if(methods::is(SE, "DeeDeeExperiment")){
    return(names(DeeDeeExperiment::getFEAList(SE, NAME)))
  }
  return(names(S4Vectors::metadata(SE)$fe_results[[NAME]]))
}


.de_result <- function(SE, NAME){
  if(methods::is(SE, "DeeDeeExperiment")){
    res <- DeeDeeExperiment::getDEA(SE, NAME, format = "original")
  }else{
    res <- S4Vectors::metadata(SE)$de_results[[NAME]]
  }
}

.fe_result <- function(SE, NAME){
  fe_list <- DeeDeeExperiment::getFEAList(SE, NAME, format = "original") %>% purrr::map(BiocGenerics::as.data.frame)
  return(fe_list)
}
