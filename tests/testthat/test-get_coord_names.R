context("test-get_coord_names")

x <- scmisc::sce
y1 <- get_coord_names(x)
x <- reduce_dim(x)
y2 <- get_coord_names(x)

test_that("get_coord_names works", {
  expect_equal(length(y1), 0)
  expect_equal(y1, NULL)
  expect_equal(length(y2), 1)
  expect_equal(y2, "pca")
})
