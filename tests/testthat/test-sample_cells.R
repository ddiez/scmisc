context("test-sample_cells")


x <- sample_cells(scmisc::sce, n = 20)

test_that("sample_cells works", {
  expect_equal(ncol(x), 20)
})

x <- cluster_cells(scmisc::sce, method = "kmeans", ncluster = 3)
x <- sample_cells(x, n = 20, group = "cluster")

test_that("sample_cells works", {
  expect_equal(ncol(x), 60)
})
