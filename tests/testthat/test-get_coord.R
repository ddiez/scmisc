context("test-get_coord")

x <- scmisc::sce
x <- reduce_dim(x, dims = 10)
y <- get_coord(x, "PCA", annotate = FALSE)

test_that("get_coord works", {
  expect_equal(dim(y), c(493, 2))
})
