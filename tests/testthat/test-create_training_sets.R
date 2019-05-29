context("test-create_training_sets")

set.seed(123)
x <- create_training_sets(scmisc::sce)

test_that("create_training_set works", {
  expect_equal(length(x), 2)
  expect_named(x)
  expect_equal(names(x), c("train", "test"))
  expect_equal(dim(x[["train"]]), c(461, 444))
  expect_equal(dim(x[["test"]]), c(461, 49))
  expect_equal(colnames(x[["train"]])[1], c("S00001_274363"))
  expect_equal(colnames(x[["test"]])[1], c("S00001_260765"))
})
