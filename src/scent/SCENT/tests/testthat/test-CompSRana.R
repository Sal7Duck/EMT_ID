test_that("CompSRana produces valid output", {
  skip_on_cran()
  
  data(net13Jun12.m, package = "SCENT")
  
  set.seed(42)
  exp <- matrix(rpois(5500 * 20, 5), nrow = 5500, ncol = 20)
  rownames(exp) <- head(rownames(net13Jun12.m), 5500)
  colnames(exp) <- paste0("Cell_", 1:20)
  
  integ <- DoIntegPPI(exp, net13Jun12.m)
  sr <- CompSRana(integ, local = FALSE)
  
  expect_type(sr, "list")
  expect_named(sr, c("SR", "inv", "locS", "nlocS"))
  expect_length(sr$SR, 20)
  expect_true(all(sr$SR >= 0 & sr$SR <= 1))
  expect_true(all(!is.na(sr$SR)))
})


test_that("CompSRana with local=TRUE returns nlocS", {
  skip_on_cran()
  
  data(net13Jun12.m, package = "SCENT")
  
  set.seed(42)
  exp <- matrix(rpois(5500 * 10, 5), nrow = 5500, ncol = 10)
  rownames(exp) <- head(rownames(net13Jun12.m), 5500)
  
  integ <- DoIntegPPI(exp, net13Jun12.m)
  sr <- CompSRana(integ, local = TRUE)
  
  expect_false(is.null(sr$nlocS))
  expect_true(all(is.finite(sr$nlocS)))
  expect_true(all(sr$nlocS >= 0))
})


test_that("CompCCAT produces valid output", {
  skip_on_cran()
  
  data(net13Jun12.m, package = "SCENT")
  
  set.seed(42)
  exp <- matrix(rpois(5500 * 50, 5), nrow = 5500, ncol = 50)
  rownames(exp) <- head(rownames(net13Jun12.m), 5500)
  
  ccat <- CompCCAT(exp, net13Jun12.m)
  
  expect_length(ccat, 50)
  expect_true(all(ccat >= -1 & ccat <= 1))
  expect_true(all(!is.na(ccat)))
})


test_that("DoIntegPPI produces valid output", {
  skip_on_cran()
  
  data(net13Jun12.m, package = "SCENT")
  
  set.seed(42)
  exp <- matrix(rpois(5500 * 10, 5), nrow = 5500, ncol = 10)
  rownames(exp) <- head(rownames(net13Jun12.m), 5500)
  
  integ <- DoIntegPPI(exp, net13Jun12.m)
  
  expect_type(integ, "list")
  expect_named(integ, c("expMC", "adjMC"))
  expect_equal(nrow(integ$expMC), nrow(integ$adjMC))
  expect_equal(nrow(integ$adjMC), ncol(integ$adjMC))
})


test_that("SR and CCAT are positively correlated", {
  skip_on_cran()
  
  data(net13Jun12.m, package = "SCENT")
  
  set.seed(123)
  exp <- matrix(rpois(5500 * 100, 5), nrow = 5500, ncol = 100)
  rownames(exp) <- head(rownames(net13Jun12.m), 5500)
  
  integ <- DoIntegPPI(exp, net13Jun12.m)
  sr <- CompSRana(integ)
  ccat <- CompCCAT(exp, net13Jun12.m)
  
  correlation <- cor(sr$SR, ccat)
  expect_true(correlation > 0.5)  # Should be positively correlated
})


test_that("Input validation works correctly", {
  skip_on_cran()
  
  data(net13Jun12.m, package = "SCENT")
  
  # Test missing rownames
  exp_no_names <- matrix(1:100, nrow = 10)
  expect_error(CompCCAT(exp_no_names, net13Jun12.m), "rownames")
  expect_error(DoIntegPPI(exp_no_names, net13Jun12.m), "rownames")
  
  # Test invalid integ.l input
  expect_error(CompSRana(list(wrong = 1)), "integ.l")
})


test_that("Sparse matrix support works", {
  skip_on_cran()
  
  data(net13Jun12.m, package = "SCENT")
  
  set.seed(42)
  exp <- matrix(rpois(5500 * 20, 5), nrow = 5500, ncol = 20)
  rownames(exp) <- head(rownames(net13Jun12.m), 5500)
  
  # Dense
  ccat_dense <- CompCCAT(exp, net13Jun12.m)
  
  # Sparse
  exp_sparse <- Matrix::Matrix(exp, sparse = TRUE)
  ccat_sparse <- CompCCAT(exp_sparse, net13Jun12.m)
  
  expect_equal(ccat_dense, ccat_sparse, tolerance = 1e-10)
})
