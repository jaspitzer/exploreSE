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


.plot_fe <- function(FE, NAME = "up_go", padj_CO = 0.05, N_terms = 5, COLS = c("Up" = "darkred", "Down" = "blue")){


  padj_col <- ifelse("gs_p.adjust" %in% names(FE), "gs_p.adjust", "p.adjust")

  if(sum(stringr::str_detect(names(FE), "NES")) > 0){

    #write code for gsea style enrichment
    gsea_df <- FE %>%
      dplyr::filter(!!rlang::sym(padj_col) < padj_CO) %>%
      dplyr::group_by(sign(NES)) %>%
      dplyr::slice_max(abs(NES), n = round(N_terms / 2)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Description = stringr::str_remove(Description, "HALLMARK_") %>%
                      stringr::str_remove("REACTOME") %>%
                      stringr::str_replace_all("_", " ") %>%
                      forcats::fct_reorder(NES),
                    dir = ifelse(NES > 0, "Up", "Down"))

    p <- ggplot2::ggplot(gsea_df, ggplot2::aes(NES, Description))+
      ggplot2::geom_col(ggplot2::aes(fill = dir))+
      ggplot2::scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 60))+
      ggplot2::scale_fill_manual(values = COLS)+
      ggplot2::theme_light(base_size = 14)+
      ggplot2::labs(x = "NES (Normalized Enrichment Score)",
                    title = if (stringr::str_detect(name, "HALLMARK")) {
                      paste("Top", N_terms, "enriched hallmark gene sets")
                    }else {
                      paste("Top", N_terms, "enriched Reactome gene sets")
                    })+
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "none")

    p_p <- plotly::ggplotly(p, tooltip = c("text")) %>%
      plotly::layout(hovermode = "closest") %>%
      plotly::style(textposition = "right")
  }else{
    #write code for GO enrichments

    filter_col <- dplyr::case_when(
      "FoldEnrichment" %in% names(FE) ~ "FoldEnrichment",
      "GeneRatio" %in% names(FE) & "BgRatio" %in% names(FE)  ~ "FoldEnrichment",
      "gs_p.adjust" %in% names(FE) ~ "gs_p.adjust",
      T ~ "p.adjust"
    )

    if(filter_col == "FoldEnrichment" & !("FoldEnrichment" %in% names(FE))){
      FE <- FE %>%
        dplyr::rowwise() %>%
        dplyr::mutate(GR = stringr::str_split(GeneRatio, "/") %>% unlist() %>% as.numeric() %>% (\(x) x[[1]] / x[[2]]),
                      BR = stringr::str_split(BgRatio, "/") %>% unlist() %>% as.numeric() %>% (\(x) x[[1]] / x[[2]])) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(FoldEnrichment = GR/BR) %>%
        dplyr::select(-GR, -BR)
    }

    if(filter_col == "FoldEnrichment"){
      df <- FE %>%
        dplyr::filter(!!rlang::sym(padj_col) < padj_CO) %>%
        dplyr::slice_max(FoldEnrichment, n = N_terms) %>%
        dplyr::mutate(Description = forcats::fct_reorder(Description, FoldEnrichment),
                      dir = ifelse(stringr::str_detect(NAME, "[Uu][Pp]"),
                                   "Up", "Down"))
      p <- ggplot2::ggplot(df, ggplot2::aes(FoldEnrichment,
                                            Description, text = stringr::str_wrap(stringr::str_replace_all(geneID, "\\/", ", "),
                                                                                  width = 60)))+
        ggplot2::geom_col(ggplot2::aes(fill = dir))+
        ggplot2::scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 60))+
        ggplot2::scale_fill_manual(values = COLS)+
        ggplot2::theme_light(base_size = 14)+
        ggplot2::labs(x = "Fold Enrichment over Background",
                      title = if (stringr::str_detect(NAME, "up")) {
                        paste("Top", N_terms, "upregulated GO Terms")
                      }else {
                        paste("Top", N_terms, "downregulated GO Terms")
                      })+
        ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "none")

      p_p <- plotly::ggplotly(p, tooltip = c("text")) %>%
        plotly::layout(hovermode = "closest") %>%
        plotly::style(textposition = "right")
    }else{
      df <- FE %>%
        dplyr::filter(!!rlang::sym(padj_col) < padj_CO) %>%
        dplyr::slice_min(!!rlang::sym(filter_col), n = N_terms) %>%
        dplyr::mutate(padj = -log10(!!rlang::sym(filter_col)),
                      Description = forcats::fct_reorder(Description, padj),
                      dir = ifelse(stringr::str_detect(NAME, "[Uu][Pp]"),
                                   "Up", "Down"))
      p <- ggplot2::ggplot(df, ggplot2::aes(padj,
                                            Description, text = stringr::str_wrap(stringr::str_replace_all(geneID, "\\/", ", "),
                                                                                  width = 60)))+
        ggplot2::geom_col(ggplot2::aes(fill = dir))+
        ggplot2::scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 60))+
        ggplot2::scale_fill_manual(values = COLS)+
        ggplot2::theme_light(base_size = 14)+
        ggplot2::labs(x = "-log10 adjusted p-Value",
                      title = if (stringr::str_detect(NAME, "[Uu][Pp]")) {
                        paste("Top", N_terms, "upregulated GO Terms")
                      }else {
                        paste("Top", N_terms, "downregulated GO Terms")
                      })+
        ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "none")

      p_p <- plotly::ggplotly(p, tooltip = c("text")) %>%
        plotly::layout(hovermode = "closest") %>%
        plotly::style(textposition = "right")
    }
    #
    # df <- FE %>%
    #   dplyr::filter(!!rlang::sym(padj_col) < padj_CO) %>%
    #   dplyr::slice_max(FoldEnrichment, n = N_terms) %>%
    #   dplyr::mutate(Description = forcats::fct_reorder(Description, FoldEnrichment),
    #                 dir = stringr::str_detect(NAME, "up"))
    # p <- ggplot2::ggplot(df, ggplot2::aes(FoldEnrichment,
    #                                           Description, text = stringr::str_wrap(stringr::str_replace_all(geneID, "\\/", ", "),
    #                                                                                 width = 60)))+
    #   ggplot2::geom_col(ggplot2::aes(fill = dir))+
    #   ggplot2::scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 60))+
    #   ggplot2::scale_fill_manual(values = c("Up" = "darkred", "Down" = "blue"))+
    #   ggplot2::theme_light(base_size = 14)+
    #   ggplot2::labs(x = "Fold Enrichment over Background",
    #                 title = if (stringr::str_detect(NAME, "up")) {
    #                   paste("Top", N_terms, "upregulated GO Terms")
    #                 }else {
    #                   paste("Top", N_terms, "downregulated GO Terms")
    #                 })+
    #   ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "none")
    #
    # p_p <- plotly::ggplotly(p, tooltip = c("text")) %>%
    #   plotly::layout(hovermode = "closest") %>%
    #   plotly::style(textposition = "right")
  }
  return(p_p)
}
