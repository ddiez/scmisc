context("test-get_rowdata")

x <- get_rowdata(scmisc::sce)

test_that("get_rowdata works", {
  expect_equal(dim(x), c(434, 4))
})
