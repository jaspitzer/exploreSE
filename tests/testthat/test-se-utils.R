## Tests for the exported orchestration functions get.gos()/get.gsea().
##
## Both functions delegate the actual enrichment computation to
## clusterProfiler (enrichGO()/GSEA()), external annotation databases and
## DeeDeeExperiment (addFEA()/renameFEA()). Rather than depending on real
## biological results (slow, non-deterministic, and would really be testing
## clusterProfiler rather than this package), we mock those calls and assert
## on the package's own logic: which genes get classified as up/down/universe,
## how the GSEA ranking statistic is computed, and how results get attached
## to the returned object.

make_de_se <- function(res_df) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      1,
      nrow = nrow(res_df),
      ncol = 1,
      dimnames = list(rownames(res_df), "s1")
    ))
  )
  S4Vectors::metadata(se)$de_results <- list(comparisonA = res_df)
  se
}

test_that("get.gos classifies up/down/universe genes and forwards them to enrichGO()", {
  res_df <- data.frame(
    padj = c(0.001, 0.001, 0.5, 0.001, NA),
    log2FoldChange = c(3, -3, 0.1, 3, 5),
    row.names = c("UP1", "DOWN1", "NS1", "UP2", "NA_GENE")
  )
  se <- make_de_se(res_df)

  enrichGO_calls <- list()
  testthat::local_mocked_bindings(
    enrichGO = function(gene, OrgDb, keyType, ont, universe) {
      enrichGO_calls[[length(enrichGO_calls) + 1]] <<- list(
        gene = gene,
        keyType = keyType,
        ont = ont,
        universe = universe
      )
      if (length(gene) == 0) {
        return(NULL)
      }
      data.frame(ID = paste0("GO:", seq_along(gene)), Description = "term")
    },
    .package = "clusterProfiler"
  )

  addFEA_calls <- list()
  renameFEA_calls <- list()
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) {
      addFEA_calls[[length(addFEA_calls) + 1]] <<- list(name = name)
      x
    },
    renameFEA = function(x, old_name, new_name) {
      renameFEA_calls[[length(renameFEA_calls) + 1]] <<- list(
        old_name = old_name,
        new_name = new_name
      )
      x
    },
    .package = "DeeDeeExperiment"
  )

  result <- get.gos("comparisonA", obj = se, species = "hs", gene_type = "SYMBOL")

  expect_true(methods::is(result, "DeeDeeExperiment"))

  # up/down classification is based on padj < 0.05 & |log2FC| > 1, NAs excluded
  up_calls <- Filter(function(x) setequal(x$gene, c("UP1", "UP2")), enrichGO_calls)
  dn_calls <- Filter(function(x) setequal(x$gene, "DOWN1"), enrichGO_calls)
  expect_true(length(up_calls) > 0)
  expect_true(length(dn_calls) > 0)

  # universe excludes genes with NA padj
  expect_true(all(vapply(
    enrichGO_calls,
    function(x) setequal(x$universe, c("UP1", "DOWN1", "NS1", "UP2")),
    logical(1)
  )))

  expect_true(all(vapply(enrichGO_calls, function(x) x$keyType == "SYMBOL", logical(1))))
  expect_true(all(vapply(enrichGO_calls, function(x) x$ont == "BP", logical(1))))

  expect_equal(length(addFEA_calls), 2)
  expect_true(all(vapply(addFEA_calls, function(x) x$name == "comparisonA", logical(1))))

  expect_equal(
    vapply(renameFEA_calls, function(x) x$old_name, character(1)),
    c("up_go", "dn_go")
  )
  expect_equal(
    vapply(renameFEA_calls, function(x) x$new_name, character(1)),
    c("comparisonA_up_go", "comparisonA_down_go")
  )
})

test_that("get.gos skips addFEA/renameFEA when a direction has no significant genes", {
  res_df <- data.frame(
    padj = c(0.001, 0.9),
    log2FoldChange = c(3, 0.1),
    row.names = c("UP1", "NS1")
  )
  se <- make_de_se(res_df)

  testthat::local_mocked_bindings(
    enrichGO = function(gene, OrgDb, keyType, ont, universe) {
      if (length(gene) == 0) {
        return(NULL)
      }
      data.frame(ID = paste0("GO:", seq_along(gene)))
    },
    .package = "clusterProfiler"
  )

  addFEA_calls <- list()
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) {
      addFEA_calls[[length(addFEA_calls) + 1]] <<- list(name = name)
      x
    },
    renameFEA = function(x, old_name, new_name) x,
    .package = "DeeDeeExperiment"
  )

  get.gos("comparisonA", obj = se, species = "hs", gene_type = "SYMBOL")

  # only the "up" direction has significant genes, so addFEA should be called once
  expect_equal(length(addFEA_calls), 1)
})

make_dds_fixture <- function() {
  counts <- matrix(
    c(
      # GENE1: clearly higher in groupA
      100, 110, 90, 10, 12, 8,
      # GENE2: clearly higher in groupB
      10, 12, 8, 100, 110, 90,
      # GENE3: similar in both groups
      50, 52, 48, 50, 51, 49
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- paste0("GENE", 1:3)
  colnames(counts) <- paste0("S", 1:6)

  coldata <- data.frame(
    sample = colnames(counts),
    condition = factor(rep(c("groupA", "groupB"), each = 3)),
    row.names = colnames(counts)
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~condition)
  DESeq2::estimateSizeFactors(dds)
}

test_that("get.gsea builds a descending gene ranking and forwards it to GSEA()", {
  dds <- make_dds_fixture()

  msigdbr_calls <- list()
  testthat::local_mocked_bindings(
    msigdbr = function(collection, db_species, ...) {
      msigdbr_calls[[length(msigdbr_calls) + 1]] <<- list(
        collection = collection,
        db_species = db_species
      )
      data.frame(
        gs_name = c("SET_A", "SET_A", "SET_B"),
        gene_symbol = c("GENE1", "GENE2", "GENE3")
      )
    },
    .package = "msigdbr"
  )

  gsea_calls <- list()
  testthat::local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      gsea_calls[[length(gsea_calls) + 1]] <<- list(
        geneList = geneList,
        TERM2GENE = TERM2GENE
      )
      structure(list(), class = "mock_gsea_result")
    },
    .package = "clusterProfiler"
  )

  addFEA_calls <- list()
  renameFEA_calls <- list()
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) {
      addFEA_calls[[length(addFEA_calls) + 1]] <<- list(name = name)
      x
    },
    renameFEA = function(x, old_name, new_name) {
      renameFEA_calls[[length(renameFEA_calls) + 1]] <<- list(
        old_name = old_name,
        new_name = new_name
      )
      x
    },
    .package = "DeeDeeExperiment"
  )

  result <- get.gsea(
    NAME = "comparisonA",
    obj = dds,
    type = "HALLMARK",
    conditions = c("groupA", "groupB"),
    species = "hs",
    condition_var = "condition"
  )

  expect_true(methods::is(result, "DeeDeeExperiment"))

  expect_equal(length(msigdbr_calls), 1)
  expect_equal(msigdbr_calls[[1]]$collection, "H")
  expect_equal(msigdbr_calls[[1]]$db_species, "HS")

  expect_equal(length(gsea_calls), 1)
  rankings <- gsea_calls[[1]]$geneList

  expect_true(is.numeric(rankings))
  expect_true(!is.null(names(rankings)))
  # groupA = "group1", groupB = "group2" -> ranking = group1 - group2, scaled
  expect_true(rankings[["GENE1"]] > rankings[["GENE2"]])
  # results are returned sorted from most to least "up"
  expect_true(all(diff(rankings) <= 0))

  expect_equal(length(addFEA_calls), 1)
  expect_equal(addFEA_calls[[1]]$name, "comparisonA")
  expect_equal(renameFEA_calls[[1]]$old_name, "gsea")
  expect_equal(renameFEA_calls[[1]]$new_name, "comparisonA_gsea_HALLMARK")
})

test_that("get.gsea requests the Reactome collection for type = 'REACTOME'", {
  dds <- make_dds_fixture()

  msigdbr_calls <- list()
  testthat::local_mocked_bindings(
    msigdbr = function(collection, subcollection = NULL, db_species, ...) {
      msigdbr_calls[[length(msigdbr_calls) + 1]] <<- list(
        collection = collection,
        subcollection = subcollection,
        db_species = db_species
      )
      data.frame(gs_name = "SET_A", gene_symbol = "GENE1")
    },
    .package = "msigdbr"
  )
  testthat::local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) structure(list(), class = "mock_gsea_result"),
    .package = "clusterProfiler"
  )
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) x,
    renameFEA = function(x, old_name, new_name) x,
    .package = "DeeDeeExperiment"
  )

  get.gsea(
    NAME = "comparisonA",
    obj = dds,
    type = "REACTOME",
    conditions = c("groupA", "groupB"),
    species = "mm",
    condition_var = "condition"
  )

  expect_equal(msigdbr_calls[[1]]$collection, "C2")
  expect_equal(msigdbr_calls[[1]]$subcollection, "CP:REACTOME")
  expect_equal(msigdbr_calls[[1]]$db_species, "MM")
})

test_that("get.gsea's default condition_var works without being supplied explicitly", {
  dds <- make_dds_fixture()

  testthat::local_mocked_bindings(
    msigdbr = function(...) data.frame(gs_name = "SET_A", gene_symbol = "GENE1"),
    .package = "msigdbr"
  )
  gsea_calls <- list()
  testthat::local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      gsea_calls[[length(gsea_calls) + 1]] <<- list(geneList = geneList)
      structure(list(), class = "mock_gsea_result")
    },
    .package = "clusterProfiler"
  )
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) x,
    renameFEA = function(x, old_name, new_name) x,
    .package = "DeeDeeExperiment"
  )

  # condition_var intentionally omitted - relies on the "condition" default
  expect_no_error(
    get.gsea(
      NAME = "comparisonA",
      obj = dds,
      type = "HALLMARK",
      conditions = c("groupA", "groupB"),
      species = "hs"
    )
  )

  expect_equal(length(gsea_calls), 1)
  expect_true(is.numeric(gsea_calls[[1]]$geneList))
})

test_that("get.gsea accepts condition_var as either a bare symbol or a string", {
  dds <- make_dds_fixture()

  testthat::local_mocked_bindings(
    msigdbr = function(...) data.frame(gs_name = "SET_A", gene_symbol = "GENE1"),
    .package = "msigdbr"
  )
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) x,
    renameFEA = function(x, old_name, new_name) x,
    .package = "DeeDeeExperiment"
  )
  gsea_calls <- list()
  testthat::local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      gsea_calls[[length(gsea_calls) + 1]] <<- list(geneList = geneList)
      structure(list(), class = "mock_gsea_result")
    },
    .package = "clusterProfiler"
  )

  expect_no_error(
    get.gsea(
      NAME = "comparisonA",
      obj = dds,
      type = "HALLMARK",
      conditions = c("groupA", "groupB"),
      species = "hs",
      condition_var = condition # bare symbol, not a string
    )
  )
  expect_no_error(
    get.gsea(
      NAME = "comparisonA",
      obj = dds,
      type = "HALLMARK",
      conditions = c("groupA", "groupB"),
      species = "hs",
      condition_var = "condition"
    )
  )

  expect_equal(length(gsea_calls), 2)
  expect_equal(gsea_calls[[1]]$geneList, gsea_calls[[2]]$geneList)
})
