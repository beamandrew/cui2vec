context("Cosine similarity")

library(cui2vec)
test_that("cosine_similarity is correct", {
  expect_equal(cosine_similarity(c(1,1), c(1,1)), 1.0)
  expect_equal(cosine_similarity(c(1,1), c(-1,-1)), -1.0)
  expect_equal(cosine_similarity(c(1,0), c(0,1)), 0.0)
})