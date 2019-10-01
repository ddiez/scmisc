context("compute_density")

x <- compute_density(mtcars[, 1], mtcars[, 2])

test_that("compute_density works", {
  expect_is(x, "numeric")
  expect_equal(length(x), nrow(mtcars))
})
