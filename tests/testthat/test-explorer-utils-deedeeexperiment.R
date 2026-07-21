## Tests for the DeeDeeExperiment code path of the explorer_utils.R helpers.
## clusterProfiler/DeeDeeExperiment internals aren't exercised here - we only need
## a real DeeDeeExperiment object so that methods::is(SE, "DeeDeeExperiment") is
## TRUE, and mock the getDEA*/getFEA* accessors to keep the tests fast and focused
## on the logic in explorer_utils.R itself.

make_empty_dde <- function() {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2"))))
  )
  DeeDeeExperiment::DeeDeeExperiment(se)
}

test_that(".check_precomputed_de is FALSE for a DeeDeeExperiment with no DEA results", {
  expect_false(.check_precomputed_de(make_empty_dde()))
})

test_that(".check_precomputed_de is TRUE when getDEANames() reports results", {
  testthat::local_mocked_bindings(
    getDEANames = function(...) c("comparisonA"),
    .package = "DeeDeeExperiment"
  )
  expect_true(.check_precomputed_de(make_empty_dde()))
})

test_that(".check_precomputed_fe is FALSE for a DeeDeeExperiment with no FEA results", {
  expect_false(.check_precomputed_fe(make_empty_dde()))
})

test_that(".check_precomputed_fe is TRUE when getFEANames() reports results", {
  testthat::local_mocked_bindings(
    getFEANames = function(...) c("comparisonA_up_go"),
    .package = "DeeDeeExperiment"
  )
  expect_true(.check_precomputed_fe(make_empty_dde()))
})

test_that(".de_results_names delegates to getDEANames() for a DeeDeeExperiment", {
  testthat::local_mocked_bindings(
    getDEANames = function(...) c("comparisonA", "comparisonB"),
    .package = "DeeDeeExperiment"
  )
  expect_equal(.de_results_names(make_empty_dde()), c("comparisonA", "comparisonB"))
})

test_that(".fe_results_names delegates to names(getFEAList()) for a DeeDeeExperiment", {
  testthat::local_mocked_bindings(
    getFEAList = function(...) list(up_go = data.frame(), dn_go = data.frame()),
    .package = "DeeDeeExperiment"
  )
  expect_equal(.fe_results_names(make_empty_dde(), "comparisonA"), c("up_go", "dn_go"))
})

test_that(".de_result strips the 'NAME_' column prefix returned by getDEA()", {
  fake_res <- data.frame(
    comparisonA_padj = c(0.01, 0.2),
    comparisonA_log2FoldChange = c(2, -0.5)
  )
  testthat::local_mocked_bindings(
    getDEA = function(...) fake_res,
    .package = "DeeDeeExperiment"
  )

  result <- .de_result(make_empty_dde(), "comparisonA")

  expect_named(result, c("padj", "log2FoldChange"))
  expect_equal(result$padj, c(0.01, 0.2))
  expect_equal(result$log2FoldChange, c(2, -0.5))
})

test_that(".de_result only strips the matching NAME prefix, not unrelated columns", {
  fake_res <- data.frame(
    comparisonA_padj = 0.01,
    unrelated_column = "kept as-is"
  )
  testthat::local_mocked_bindings(
    getDEA = function(...) fake_res,
    .package = "DeeDeeExperiment"
  )

  result <- .de_result(make_empty_dde(), "comparisonA")

  expect_named(result, c("padj", "unrelated_column"))
})

test_that(".fe_result converts every entry of the FEA list to a data frame", {
  fake_fea <- list(
    up_go = data.frame(Description = "term A", p.adjust = 0.01),
    dn_go = data.frame(Description = "term B", p.adjust = 0.02)
  )
  testthat::local_mocked_bindings(
    getFEAList = function(...) fake_fea,
    .package = "DeeDeeExperiment"
  )

  result <- .fe_result(make_empty_dde(), "comparisonA")

  expect_named(result, c("up_go", "dn_go"))
  expect_true(all(vapply(result, is.data.frame, logical(1))))
  expect_equal(result$up_go$Description, "term A")
})
