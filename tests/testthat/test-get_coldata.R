context("test-get_coldata")

x <- get_coldata(scmisc::sce)

test_that("get_coldata works", {
  expect_equal(dim(x), c(493, 3))
})
