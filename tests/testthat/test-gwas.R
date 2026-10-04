test_that("GWAS columns are recognised (GAPIT style)", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("SNP,Chr,Pos,P.value,MAF", "S1_100,1,100,0.5,0.2", "S1_200,1,200,1e-8,0.3",
               "S2_50,2,50,0.01,0.1"), f)
  gw <- read_gwas(f, verbose = FALSE)
  expect_s3_class(gw, "qp_gwas")
  expect_equal(nrow(gw), 3L)
  expect_equal(unname(attr(gw, "columns")), c("SNP", "Chr", "Pos", "P.value"))
  expect_equal(gw$logp[gw$snp == "S1_200"], 8)
  th <- gwas_threshold(gw)
  expect_equal(th$p, 0.05 / 3)
})
