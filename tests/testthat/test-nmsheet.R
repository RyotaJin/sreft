test_that("slope", {
  expect_equal(slope(1:5, 1:5 * 1.5 + 3), 1.5)
})

test_that("intercept", {
  expect_equal(intercept(1:5, 1:5 * 1.5 + 3), 3)
})
