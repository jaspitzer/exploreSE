test_that(".create_demo_data returns a well-formed SummarizedExperiment", {
  se <- .create_demo_data()

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(1000, 12))
  expect_true("counts" %in% SummarizedExperiment::assayNames(se))
})

test_that(".create_demo_data counts are non-negative integers with the expected dimnames", {
  se <- .create_demo_data()
  counts <- SummarizedExperiment::assay(se, "counts")

  expect_true(all(counts >= 0))
  expect_equal(counts, round(counts))
  expect_equal(rownames(counts), paste0("GENE", seq_len(1000)))
  expect_equal(colnames(counts), paste0("Sample", seq_len(12)))
})

test_that(".create_demo_data sets up condition/batch/replicate colData", {
  se <- .create_demo_data()
  cd <- SummarizedExperiment::colData(se)

  expect_equal(colnames(cd), c("condition", "batch", "replicate"))
  expect_equal(as.character(cd$condition), rep(c("Control", "Treatment"), each = 6))
  expect_equal(as.character(cd$batch), rep(c("A", "B"), times = 6))
  expect_equal(cd$replicate, rep(seq_len(6), times = 2))
})

test_that(".create_demo_data sets up gene_id/gene_name rowData", {
  se <- .create_demo_data()
  rd <- SummarizedExperiment::rowData(se)

  expect_equal(rd$gene_id, rownames(se))
  expect_equal(rd$gene_name, paste0("Gene_", seq_len(1000)))
})
