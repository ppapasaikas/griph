test_that("HellingerMat works properly", {
    x <- matrix(0:8, nrow = 3)
    expect_is(HellingerMat(x), "matrix")
    expect_equal(sum(HellingerMat(x)), 1.617188192, tolerance = 1.0e-8)
    expect_equal(sum(is.finite(HellingerMat(x - 3))), 9)
})
