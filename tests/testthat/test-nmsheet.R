test_that("slope", {
  expect_equal(slope(1:5, 1:5 * 1.5 + 3), 1.5)
})

test_that("intercept", {
  expect_equal(intercept(1:5, 1:5 * 1.5 + 3), 3)
})

test_that("scoreConveter", {
  expect_equal(scoreConverter(10, 10), 21)
})
