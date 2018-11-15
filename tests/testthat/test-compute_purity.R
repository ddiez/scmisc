context("test-compute_purity")

set.seed(123)
x <- cluster_cells(scmisc::sce, method = "kmeans", ncluster = 5, column.name = "kmeans")
x <- cluster_cells(x, method = "hclust", ncluster = 5, column.name = "hclust")
d <- compute_purity(x, "kmeans", "hclust")

test_that("compute_purity works", {
  expect_is(d, "data.frame")
  expect_equal(dim(d), c(13, 4))
  expect_equal(colnames(d), c("kmeans", "hclust", "count", "purity"))
  expect_equal(d[["purity"]][1], 0.0175, tolerance = 0.0001)
})
