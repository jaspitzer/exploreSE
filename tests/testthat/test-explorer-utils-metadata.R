## Tests for the "plain SummarizedExperiment + metadata()$de_results / $fe_results"
## code path (as opposed to the DeeDeeExperiment path, covered separately).

make_se_with_metadata <- function(de_results = NULL, fe_results = NULL) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2"))))
  )
  if (!is.null(de_results)) {
    S4Vectors::metadata(se)$de_results <- de_results
  }
  if (!is.null(fe_results)) {
    S4Vectors::metadata(se)$fe_results <- fe_results
  }
  se
}

test_that(".check_precomputed_de is FALSE for NULL, empty and missing metadata", {
  expect_false(.check_precomputed_de(NULL))
  expect_false(.check_precomputed_de(make_se_with_metadata()))
  expect_false(.check_precomputed_de(make_se_with_metadata(de_results = list())))
})

test_that(".check_precomputed_de is TRUE when de_results is a non-empty list", {
  de <- list(comparisonA = data.frame(padj = 0.01, log2FoldChange = 2))
  expect_true(.check_precomputed_de(make_se_with_metadata(de_results = de)))
})

test_that(".check_precomputed_fe is FALSE for NULL, empty and missing metadata", {
  expect_false(.check_precomputed_fe(NULL))
  expect_false(.check_precomputed_fe(make_se_with_metadata()))
  expect_false(.check_precomputed_fe(make_se_with_metadata(fe_results = list())))
})

test_that(".check_precomputed_fe is TRUE when fe_results is a non-empty list", {
  fe <- list(comparisonA = list(up_go = data.frame(Description = "term")))
  expect_true(.check_precomputed_fe(make_se_with_metadata(fe_results = fe)))
})

test_that(".de_results_names reads names from metadata$de_results", {
  de <- list(comparisonA = data.frame(), comparisonB = data.frame())
  se <- make_se_with_metadata(de_results = de)
  expect_equal(.de_results_names(se), c("comparisonA", "comparisonB"))
})

test_that(".de_results_names is empty/NULL when no de_results are present", {
  expect_null(.de_results_names(make_se_with_metadata()))
})

test_that(".fe_results_names reads names from metadata$fe_results[[NAME]]", {
  fe <- list(comparisonA = list(up_go = data.frame(), dn_go = data.frame()))
  se <- make_se_with_metadata(fe_results = fe)
  expect_equal(.fe_results_names(se, "comparisonA"), c("up_go", "dn_go"))
})

test_that(".fe_results_names is NULL for an unknown comparison", {
  fe <- list(comparisonA = list(up_go = data.frame()))
  se <- make_se_with_metadata(fe_results = fe)
  expect_null(.fe_results_names(se, "does_not_exist"))
})

test_that(".de_result returns the stored data frame unchanged on the metadata path", {
  res_df <- data.frame(padj = c(0.01, 0.2), log2FoldChange = c(2, -0.5))
  se <- make_se_with_metadata(de_results = list(comparisonA = res_df))

  expect_equal(.de_result(se, "comparisonA"), res_df)
})
