context("reduce_dim")

x <- as.matrix(mtcars)
y <- reduce_dim(x)

test_that("reduce_dim works", {
  expect_equal(dim(y), c(11, 2))
  expect_equal(colnames(y), c("dim1", "dim2"))
  expect_equal(rownames(y), colnames(x))
})

y <- reduce_dim(x, dims = 3)

test_that("reduce_dim works", {
  expect_equal(dim(y), c(11, 3))
  expect_equal(colnames(y), c("dim1", "dim2", "dim3"))
  expect_equal(rownames(y), colnames(x))
})

y <- reduce_dim(x, "tsne", perplexity = 3)

test_that("reduce_dim works", {
  expect_equal(dim(y), c(11, 2))
  expect_equal(colnames(y), c("dim1", "dim2"))
  expect_equal(rownames(y), colnames(x))
})
