test_that("griph_cluster works properly", {
    set.seed(1)
    M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
    label <- attr(M, "label")
    i <- sample(ncol(M), 50)
    # serial
    set.seed(2)
    res <- griph_cluster(M[,i], ClassAssignment = label[i], ref.iter = 1,
                         use.par = FALSE, filter = TRUE, plot_ = FALSE)
    expect_is(res, "list")
    expect_is(res$GRAO, "igraph")
    expect_length(res$MEMB, length(i))
    expect_length(res$MEMB.true, length(i))
    # expect_lt(res$miscl, 0.3)
    expect_equal(dim(res$ConfMatrix), c(length(unique(res$MEMB)), nlevels(label)))
    # parallel
    set.seed(2)
    resP <- griph_cluster(M[,i], ClassAssignment = label[i], ref.iter = 1,
                          use.par = TRUE, filter = TRUE, plot_ = FALSE)
    expect_equal(res$MEMB, resP$MEMB)
})
