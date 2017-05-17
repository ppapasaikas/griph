test_that("griph_cluster works properly", {
    M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
    lab <- attr(M, "label")
    res1 <- griph_cluster(M, ClassAssignment = lab, plotG = FALSE, ref.iter=2)
    expect_is(res1, "list")
    #Inonsistencies in code below due to random sampling
    #expect_equal(griph:::classError(res1$MEMB, lab), 61/288)
    # currently, check does freeze if running with use.par=TRUE
    #res2 <- griph_cluster(M, ClassAssignment = lab, plotG = FALSE, ref.iter=2, use.par = TRUE, ncores = 4)
    #expect_identical(res1$MEMB, res2$MEMB)
    #expect_equal(res1$DISTM, res2$DISTM)
    #expect_true(igraph::isomorphic(res1$GRAO, res2$GRAO))
})
