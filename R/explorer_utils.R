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
    res <- DeeDeeExperiment::getDEA(SE, NAME, type = "data.frame") %>%
      dplyr::rename_with(.fn = \(x) stringr::str_remove(x, paste0(NAME, "_")))
  }else{
    res <- S4Vectors::metadata(SE)$de_results[[NAME]]
  }
}

.fe_result <- function(SE, NAME){
  fe_list <- DeeDeeExperiment::getFEAList(SE, NAME, format = "original") %>% purrr::map(BiocGenerics::as.data.frame)
  return(fe_list)
}


.plot_des <- function(RES, NAME, padj_CO = 0.05, fc_CO = 1, COLS = c("darkred", "darkblue")){
  p <- RES %>%
    dplyr::filter(padj < padj_CO,
                  abs(fc_CO) >= 1) %>%
    dplyr::mutate(direction = ifelse(log2FoldChange > fc_CO, "Up", "Down")) %>%
    dplyr::count(direction) %>%
    ggplot2::ggplot(ggplot2::aes(direction, y = n, fill = direction))+
    ggplot2::geom_col()+
    ggplot2::geom_label(ggplot2::aes(label = n), fill = "white", show.legend = F)+
    ggplot2::scale_fill_manual(values = COLS)+
    ggplot2::labs(title = paste("DE genes in", NAME),)+
    ggplot2::theme_light(base_size = 14)+
    ggplot2::theme(legend.position = "none")

  # p_p <- plotly::ggplotly(p)

  return(p)
}

.plot_fe <- function(FE, NAME = "up_go", padj_CO = 0.05, N_terms = 5, COLS = c("Up" = "darkred", "Down" = "blue")){


  # this can probably go when the DeeDee devel branch hits main, but for now needs to stay
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

    # this can probably go, when the DeeDee devel branch hits main, but for now needs to stay
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
  }
  return(p_p)
}


.plot_volcano <- function(RES, NAME, padj_CO = 0.05, fc_CO = 1, highlights = NULL,
                          COLS, LABEL_TOP = T, TOPN= 5){
  # Prepare data
  volcano_df <- RES %>%
    dplyr::filter(!is.na(padj) & !is.na(log2FoldChange) & !is.infinite(log2FoldChange)) %>%
    dplyr::mutate(
      highlighted = gene_id %in% highlights,
      sig = dplyr::case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "Up",
        padj < 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE ~ "NS"
      )) %>%
    dplyr::mutate(
      display_category = dplyr::case_when(
        highlighted ~ "Highlighted",
        TRUE ~ sig
      )
    )

  # Identify top genes to label (excluding highlighted genes since they'll be labeled anyway)
  genes_to_label <- data.frame()
  if (LABEL_TOP && TOPN > 0) {
    top_genes <- volcano_df %>%
      dplyr::filter(sig != "NS" & !highlighted) %>%
      dplyr::arrange(padj) %>%
      utils::head(TOPN)
    genes_to_label <- rbind(genes_to_label, top_genes)
  }

  # Always label highlighted genes
  if (length(highlights) > 0) {
    highlighted_genes <- volcano_df %>%
      dplyr::filter(highlighted)
    genes_to_label <- rbind(genes_to_label, highlighted_genes)
  }

  # Reorder so highlighted genes are plotted on top
  volcano_df <- volcano_df %>%
    dplyr::arrange(highlighted)

  # Create plot
  p <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = log2FoldChange, y = -log10(padj),
                                                color = display_category, text = paste("log2 FC:", round(log2FoldChange, 2), "\nGene:", gene_id, "\nadjusted p-value:", format(padj, digits=4)))) +
    ggplot2::geom_point(ggplot2::aes(size = highlighted, alpha = ifelse(highlighted, 1, 0.6))) +
    ggplot2::scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.5), guide = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_color_manual(values = COLS,
                                name = "Category",
                                breaks = c("Up", "Down", "Highlighted", "NS"),
                                labels = c("Up-regulated", "Down-regulated", "Highlighted", "Not significant")) +
    ggplot2::geom_vline(xintercept = c(-fc_CO, fc_CO),
                        linetype = "dashed", color = "grey30", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(padj_CO),
                        linetype = "dashed", color = "grey30", linewidth = 0.5) +
    ggplot2::labs(
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      title = if (!is.null(NAME)) {
        paste("Volcano Plot:", NAME)
      } else {
        "Volcano Plot"
      }
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  # Add labels
  if (nrow(top_genes) > 0) {
    p <- p + ggplot2::geom_text(
      data = top_genes,
      ggplot2::aes(label = gene_id),
      size = 5,
      show.legend = FALSE
    )
  }

  p_p <- plotly::ggplotly(p, tooltip = c("text")) %>%
    plotly::layout(hovermode = "closest") %>%
    plotly::style(textposition = "right")

  return(p_p)
}
