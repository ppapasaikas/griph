test_that("classError and mapLabelsGreedy works properly", {
    a <- rep(1:3, c(10,20,30))
    b <- a
    set.seed(0)
    b[sample(60, 20, TRUE)] <- sample(1:3, 20, TRUE)
    b[sample(60, 10, TRUE)] <- 4
    map <- mapLabelsGreedy(a, b)
    expect_is(map, "character")
    expect_equal(classError(a, b, map), 22/60, tolerance = 1.0e-8)
    expect_equal(classError(a, b, exhaustive = TRUE), classError(a, b))
})
