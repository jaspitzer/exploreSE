

#' Title
#'
#' @param NAME name of the de comparison
#' @param obj the object
#' @param species the species hs/mm
#' @param gene_type the type of gene identifier
#'
#' @returns a DeeDeeExperiment or SummarizedExperiment
#' @export
#'
#' @examples
#' TRUE
get.gos <- function(NAME, obj = dds, species = "hs", gene_type = "SYMBOL"){
  if(species == "hs"){
    DB <- org.Hs.eg.db::org.Hs.eg.db
  }else {
    DB <- org.Mm.eg.db::org.Mm.eg.db
  }

  res <- .de_result(obj, NAME)

  up_genes <- res %>% BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(padj < 0.05, log2FoldChange > 1) %>% dplyr::pull(gene)
  dn_genes <- res %>% BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(padj < 0.05, log2FoldChange < -1) %>% dplyr::pull(gene)

  universe <- res %>% BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(padj)) %>% dplyr::pull(gene)

  up_go <- clusterProfiler::enrichGO(up_genes, DB, keyType = gene_type, ont = "BP", universe = universe) %>%
    BiocGenerics::as.data.frame()
  dn_go <- clusterProfiler::enrichGO(dn_genes, DB, keyType = gene_type, ont = "BP", universe = universe) %>%
    BiocGenerics::as.data.frame()

  if(!methods::is(obj, "DeeDeeExperiment")){
    obj <- DeeDeeExperiment::DeeDeeExperiment(obj)
  }

  if(length(up_genes) > 0){
    up_go <- clusterProfiler::enrichGO(up_genes, DB, keyType = gene_type, ont = "BP", universe = universe)
    if(!is.null(up_go)){
      obj <- DeeDeeExperiment::addFEA(obj, up_go, NAME)
      obj <- DeeDeeExperiment::renameFEA(obj, "up_go", paste0(NAME, "_up_go"))
    }
  }

  if(length(dn_genes) > 0){
    dn_go <- clusterProfiler::enrichGO(dn_genes, DB, keyType = gene_type, ont = "BP", universe = universe)
    if(!is.null(dn_go)){
      obj <- DeeDeeExperiment::addFEA(obj, dn_go, NAME)
      obj <- DeeDeeExperiment::renameFEA(obj, "dn_go", paste0(NAME, "_down_go"))
    }
  }

  # if(!is.list(metadata(obj)$fe_results)){
  #   S4Vectors::metadata(obj)$fe_results <- list()
  # }
  # if(!is.list(metadata(obj)$fe_results[[NAME]])){
  #   S4Vectors::metadata(obj)$fe_results[[NAME]] <- list()
  # }
  #
  # S4Vectors::metadata(obj)$fe_results[[NAME]][["up_go"]] <- up_go
  # S4Vectors::metadata(obj)$fe_results[[NAME]][["dn_go"]] <- dn_go
  return(obj)
}
## EXAMPLE
# dds_2 <- get.gos("KO effect in untreated")



#' Title
#'
#' @param NAME name of the DE comparison
#' @param obj the SE/DeeDeeExperiemnt
#' @param type hallmark or reactome
#' @param conditions which conditions are compared
#' @param species what species
#' @param condition_var whats the variable name for this comparison
#'
#' @returns a DeeDeeExperiment or SummarizedExperiment
#' @export
#'
#' @examples
#' TRUE
get.gsea <- function(NAME, obj = dds, type = "HALLMARK", conditions, species = "hs", condition_var = condition){
  if(species == "hs"){
    SPECIES <- "HS"
  }else {
    SPECIES <- "MM"
  }

  if(type == "HALLMARK"){
    gene_sets <- msigdbr::msigdbr(collection = "H", db_species = SPECIES)
  }else if(type == "REACTOME"){
    gene_sets <- msigdbr::msigdbr(collection = "C2", subcollection = "CP:REACTOME", db_species = SPECIES)
  }

  gene_sets <- gene_sets %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::distinct()

  rankings <- BiocGenerics::counts(obj, normalized = T) %>%
    BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene) %>%
    dplyr::left_join(BiocGenerics::as.data.frame(SummarizedExperiment::colData(obj)), dplyr::join_by(name == sample)) %>%
    dplyr::mutate(condition = !!rlang::sym(condition_var)) %>%
    dplyr::filter(condition %in% conditions)  %>%
    dplyr::select(gene, condition, value) %>%
    dplyr::group_by(condition, gene) %>%
    dplyr::summarise(sd = stats::sd(value, na.rm = T), mean_expr = mean(value, na.rm = T)) %>%
    dplyr::mutate(sd = ifelse(sd >= 0.2 * abs(mean_expr), sd,  0.2 * abs(mean_expr)),
           mean_expr = ifelse(mean_expr == 0, 1, mean_expr))%>%
    dplyr::mutate(condition = ifelse(condition == conditions[2], "group2", "group1")) %>%
    dplyr::select(condition, mean_expr:sd, gene) %>%
    tidyr::pivot_wider(names_from = condition,  values_from = mean_expr:sd) %>%
    dplyr::filter(dplyr::if_all(sd_group1:sd_group2,\(x) x > 0),
                  dplyr::if_any(mean_expr_group1:mean_expr_group2,\(x) x > 0)) %>%
    dplyr::mutate(ranking = (mean_expr_group1 - mean_expr_group2) / (sd_group1 + sd_group2),
           gene = forcats::fct_reorder(gene, ranking, .desc = T)) %>%
    dplyr::filter(!is.na(ranking))%>%
    dplyr::arrange(dplyr::desc(ranking)) %>%
    dplyr::pull(ranking, name = gene)

  gsea <- clusterProfiler::GSEA(rankings, TERM2GENE = gene_sets)

  if(!methods::is(obj, "DeeDeeExperiment")){
    obj <- DeeDeeExperiment::DeeDeeExperiment(obj)
  }

  obj <- DeeDeeExperiment::addFEA(obj, gsea, NAME)
  obj <- DeeDeeExperiment::renameFEA(obj, "gsea", paste0(NAME, "_gsea_", type))


  return(obj)
}

## Example
#
# dds_3 <- get.gsea(NAME = "KO effect in untreated", obj = dds_2, type = "HALLMARK", conditions = c("ko_untreated", "wt_untreated"))
# dds_3 <- get.gsea(NAME = "KO effect in untreated", obj = dds_3, type = "REACTOME", conditions = c("ko_untreated", "wt_untreated"))
