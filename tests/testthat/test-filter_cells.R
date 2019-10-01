context("filter_cells")

x <- scmisc::sce
x <- filter_cells(x)

test_that("filter_cells works", {
  expect_equal(dim(x), c(434, 493))
})

x <- filter_cells(x, min.genes = 100)

test_that("filter_cells works", {
  expect_equal(dim(x), c(434, 430))
})

x <- filter_cells(x, min.count = 10)

test_that("filter_cells works", {
  expect_equal(dim(x), c(434, 427))
})
