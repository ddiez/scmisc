context("filter_cells")

x <- scmisc::sce
x <- filter_cells(x)

test_that("filter_cells works", {
  expect_equal(dim(x), c(461, 493))
})

x <- filter_cells(x, min.gene = 100)

test_that("filter_cells works", {
  expect_equal(dim(x), c(461, 430))
})

x <- filter_cells(x, min.count = 10)

test_that("filter_cells works", {
  expect_equal(dim(x), c(461, 427))
})
