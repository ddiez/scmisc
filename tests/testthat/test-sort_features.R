context("test-sort_features")

features <- sort_features(scmisc::sce)

test_that("sort_features works", {
  expect_equal(length(features), 434)
  expect_equal(features[1:3], c("Tlr9", "Birc3", "Rel"))
})
