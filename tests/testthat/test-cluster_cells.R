context("test-cluster_cells")

set.seed(123)
x <- cluster_cells(scmisc::sce, method = "kmeans", ncluster = 5)

test_that("cluster_cells works", {
  expect_false(is.null(x[["cluster"]]))
  expect_equal(nlevels(x[["cluster"]]), 5)
  expect_equal(as.character(x[["cluster"]][1]), "4")
})
