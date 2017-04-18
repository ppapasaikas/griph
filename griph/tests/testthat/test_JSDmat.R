test_that("JSDmat works properly", {
    x <- matrix(0:8, nrow = 3)
    expect_is(JSDmat(x), "matrix")
    expect_equal(sum(JSDmat(x)), 0.4412547, tolerance = 1.0e-8)
    expect_equal(sum(is.finite(JSDmat(x - 3))), 9)
})
