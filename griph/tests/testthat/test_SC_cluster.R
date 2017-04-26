test_that("SC_cluster works properly", {
    M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
    lab <- attr(M, "label")
    res1 <- SC_cluster(M, ClassAssignment = lab, plotG = FALSE)
    expect_is(res1, "list")
    expect_equal(griph:::classError(res1$MEMB, lab), 50/288)
    # currently, check does freeze if running with use.par=TRUE
    #res2 <- SC_cluster(M, ClassAssignment = lab, plotG = FALSE, use.par = TRUE, ncores = 4)
    #expect_identical(res1$MEMB, res2$MEMB)
    #expect_equal(res1$DISTM, res2$DISTM)
    #expect_true(igraph::isomorphic(res1$GRAO, res2$GRAO))
})
