context("test-compute_jaccard")

set.seed(123)
x <- cluster_cells(scmisc::sce, method = "kmeans", ncluster = 5, column.name = "kmeans")
x <- cluster_cells(x, method = "hclust", ncluster = 5, column.name = "hclust")
d <- compute_jaccard(x, "kmeans", "hclust")

test_that("compute_jaccard works", {
  expect_is(d, "data.frame")
  expect_equal(dim(d), c(11, 4))
  expect_equal(colnames(d), c("kmeans", "hclust", "count", "jaccard"))
  expect_equal(d[["jaccard"]][1], 0.0211, tolerance = 0.0001)
})
