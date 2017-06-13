test_that("PHellingerMat works properly", {
    x <- matrix(0:8, nrow = 3)
    res1 <- PHellingerMat(x, x)
    res2 <- PHellingerMat(x - 3, x - 3)
    expect_is(res1, "matrix")
    expect_equal(sum(res1), 1.617188192, tolerance = 1.0e-8)
    expect_equal(sum(is.finite(res2)), 9)
})

test_that("PHellingerMatOMP works properly", {
    x <- matrix(0:8, nrow = 3)
    res1 <- PHellingerMat(x, x)
    res2 <- PHellingerMatOMP(x, x)
    expect_identical(res1, res2)
})

test_that("PCanberraMat works properly", {
    x <- matrix(0:8, nrow = 3)
    res1 <- PCanberraMat(x, x)
    res2 <- PCanberraMat(x - 3, x - 3)
    expect_is(res1, "matrix")
    expect_equal(sum(res1), 10.4308025308, tolerance = 1.0e-8)
    expect_equal(sum(is.finite(res2)), 7)
})

test_that("PCanberraMatOMP works properly", {
    x <- matrix(0:8, nrow = 3)
    res1 <- PCanberraMat(x, x)
    res2 <- PCanberraMatOMP(x, x)
    expect_identical(res1, res2)
})
