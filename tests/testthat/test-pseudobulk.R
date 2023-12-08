context("test-pseudobulk")

m <- matrix(c(1:24), nrow=3, dimnames=list(paste0("gene-", 1:3), paste0("cell-", 1:8)))
samplename <- rep(c("D1", "D2"), 4)
celltype <- rep(c("A", "B"), each=4)
meta <- data.frame(samplename=samplename, celltype=celltype)

x <- CreateSeuratObject(as.sparse(m), meta=meta)

samples <- data.frame(sex=c("F", "M"), row.names=c("D1", "D2"))

y <- pseudobulk(x, split.by="celltype", group.by="samplename", assay="RNA", samples=samples)

test_that("sort_features works", {
  expect_equal(length(y), 2)
  expect_equal(dim(y[[1]]), c(3, 2))
  expect_equal(dim(y[[2]]), c(3, 2))
  expect_equal(dim(y[[1]]$samples), c(2, 5))
  expect_equal(dim(y[[2]]$samples), c(2, 5))
  expect_null(y[[1]]$genes)
  expect_null(y[[2]]$genes)
})

genes <- data.frame(row.names=rownames(m))
y <- pseudobulk(x, split.by="celltype", group.by="samplename", assay="RNA", genes=genes)

test_that("sort_features works", {
  expect_equal(length(y), 2)
  expect_equal(dim(y[[1]]), c(3, 2))
  expect_equal(dim(y[[2]]), c(3, 2))
  expect_equal(dim(y[[1]]$samples), c(2, 4))
  expect_equal(dim(y[[2]]$samples), c(2, 4))
  expect_equal(dim(y[[1]]$genes), c(3, 0))
  expect_equal(dim(y[[2]]$genes), c(3, 0))
})
