context("filter_genes")

x <- scmisc::sce
x <- filter_genes(x)

test_that("filter_genes works", {
  expect_equal(dim(x), c(387, 493))
})

x <- filter_genes(x, min.cell = 100)

test_that("filter_genes works", {
  expect_equal(dim(x), c(210, 493))
})

x <- filter_genes(x, min.count = 10)

test_that("filter_genes works", {
  expect_equal(dim(x), c(94, 493))
})
