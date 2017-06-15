test_that("get.knn works properly", {
    set.seed(0)
    M <- as(matrix(rnorm(100), nrow=10), "sparseMatrix")
    res1 <- get.knn(M)
    res2 <- get.knn(M, 2)
    expect_is(res1, "matrix")
    expect_equal(sum(res1), 159)
    expect_equal(res1[-1,], res2)
})

test_that("keep.mknn works properly", {
    set.seed(0)
    M <- as(matrix(rnorm(100), nrow=10), "sparseMatrix")
    res1 <- keep.mknn(M, 3)
    res2 <- keep.mknn(M, 4)
    expect_is(res1, "Matrix")
    expect_equal(sum(res1 > 0), 12)
    expect_true(all(which(as.logical(res1 > 0)) %in% which(as.logical(res2 > 0))))
})

test_that("sparsify works properly", {
    set.seed(0)
    M <- as(matrix(abs(rnorm(100)), nrow=10), "sparseMatrix")
    res <- sparsify(M, 0.5)
    expect_is(res, "Matrix")
    expect_equal(sum(res > 0), 50)
})

test_that("sparse.cor works properly", {
    set.seed(7)
    M <- as(matrix(abs(rnorm(200)), nrow=10), "sparseMatrix")
    res1 <- sparse.cor(M)
    res2 <- cor(as.matrix(M))
    expect_is(res1, "matrix")
    expect_equal(dim(res1), dim(res2))
    expect_true(all(abs(res1-res2) < 1e-6))
})

test_that("psparse.cor works properly", {
    set.seed(3)
    M1 <- as(matrix(abs(rnorm(50)), nrow=10), "sparseMatrix")
    M2 <- as(matrix(abs(rnorm(50)), nrow=10), "sparseMatrix")
    res1 <- psparse.cor(M1, M2)
    res2 <- cor(as.matrix(M1), as.matrix(M2))
    expect_is(res1, "matrix")
    expect_equal(dim(res1), dim(res2))
    expect_true(all(abs(res1-res2) < 1e-6))
})

test_that("gg_color_hue works properly", {
    expect_is(gg_color_hue(10), "character")
    expect_length(gg_color_hue(7), 7)
    expect_match(gg_color_hue(10), "^#[0-9A-F]{6}$", all = TRUE)
})

test_that("adjust.color works properly", {
    cols <- c("red","green","blue","orange")
    cols2 <- adjust.color(cols, 1.5)
    expect_is(cols2, "character")
    expect_length(cols2, 4)
    expect_identical(cols2, c("#FFAAAA","#AAFFAA","#AAAAFF","#FFE1AA"))
})

test_that("RandString works properly", {
    expect_is(RandString(n=1, len = 10), "character")
    expect_length(RandString(7), 7)
    expect_match(RandString(10, 3), "^[0-9A-Za-z]{3}$", all = TRUE)
})

