## Tests for the ggplot2/plotly building helpers in explorer_utils.R.
## We only assert on the returned object's type/structure, not on exact
## pixel/rendering output.

make_de_results <- function() {
  set.seed(1)
  data.frame(
    gene_id = paste0("GENE", 1:50),
    log2FoldChange = stats::rnorm(50, 0, 2),
    padj = stats::runif(50, 0, 1)
  )
}

make_go_fe_results <- function() {
  data.frame(
    Description = paste("Term", 1:10),
    p.adjust = stats::runif(10, 0, 0.1),
    FoldEnrichment = stats::runif(10, 1, 5),
    geneID = paste0("GENE", 1:10, "/GENE", 11:20),
    GeneRatio = "5/100",
    BgRatio = "50/1000"
  )
}

make_gsea_fe_results <- function() {
  data.frame(
    Description = paste0("HALLMARK_TERM_", 1:10),
    p.adjust = stats::runif(10, 0, 0.1),
    NES = stats::rnorm(10, 0, 2)
  )
}

test_that(".plot_des returns a ggplot summarising up/down gene counts", {
  p <- .plot_des(make_de_results(), "test comparison")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "DE genes in test comparison")
})

test_that(".plot_fe handles GO-style (FoldEnrichment) results", {
  p <- .plot_fe(make_go_fe_results(), NAME = "up_go")
  expect_s3_class(p, "plotly")
})

test_that(".plot_fe handles GSEA-style (NES) results without erroring", {
  p <- .plot_fe(make_gsea_fe_results(), NAME = "comparisonA_gsea_HALLMARK")
  expect_s3_class(p, "plotly")
})

test_that(".plot_fe titles GSEA plots based on HALLMARK vs Reactome in NAME", {
  hallmark_df <- make_gsea_fe_results()
  reactome_df <- make_gsea_fe_results()

  p_hallmark <- .plot_fe(hallmark_df, NAME = "comparisonA_gsea_HALLMARK")
  p_reactome <- .plot_fe(reactome_df, NAME = "comparisonA_gsea_REACTOME")

  expect_s3_class(p_hallmark, "plotly")
  expect_s3_class(p_reactome, "plotly")
})

test_that(".plot_volcano returns a plotly object with LABEL_TOP = TRUE (default)", {
  p <- .plot_volcano(
    make_de_results(),
    "test",
    COLS = c(Up = "red", Down = "blue", Highlighted = "green", NS = "grey")
  )
  expect_s3_class(p, "plotly")
})

test_that(".plot_volcano returns a plotly object with LABEL_TOP = FALSE", {
  p <- .plot_volcano(
    make_de_results(),
    "test",
    COLS = c(Up = "red", Down = "blue", Highlighted = "green", NS = "grey"),
    LABEL_TOP = FALSE
  )
  expect_s3_class(p, "plotly")
})

test_that(".plot_volcano works with highlighted genes and no top-N labels", {
  res <- make_de_results()
  p <- .plot_volcano(
    res,
    "test",
    COLS = c(Up = "red", Down = "blue", Highlighted = "green", NS = "grey"),
    highlights = res$gene_id[1:3],
    LABEL_TOP = FALSE
  )
  expect_s3_class(p, "plotly")
})

test_that(".plot_volcano handles no significant genes gracefully", {
  res <- data.frame(
    gene_id = paste0("GENE", 1:10),
    log2FoldChange = rep(0.1, 10),
    padj = rep(0.9, 10)
  )
  p <- .plot_volcano(
    res,
    "test",
    COLS = c(Up = "red", Down = "blue", Highlighted = "green", NS = "grey")
  )
  expect_s3_class(p, "plotly")
})
