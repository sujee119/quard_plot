test_that("genotype r2 equals the squared correlation of dosages", {
  set.seed(1)
  D <- matrix(sample(0:2, 40 * 6, replace = TRUE), nrow = 6)
  g <- as_genotypes(D, chr = 1, pos = seq(100, 600, by = 100), samples = paste0("s", 1:40))
  ld <- calc_ld(g)
  expect_equal(unname(unclass(ld)[1:6, 1:6]), unname(stats::cor(t(D))^2), tolerance = 1e-8,
               ignore_attr = TRUE)
  dp <- calc_ld(g, stat = "dprime")
  expect_true(all(dp >= 0 & dp <= 1 + 1e-9, na.rm = TRUE))
})

test_that("LD blocks are found and the region is the block of the target", {
  ex <- qp_example()
  g <- genotype_qc(read_vcf(ex$vcf, chr = 9, verbose = FALSE), verbose = FALSE)
  bl <- find_ld_blocks(g, verbose = FALSE)
  expect_gt(nrow(bl), 0L)
  expect_true(all(bl$end >= bl$start))
  reg <- define_region(9, ex$pos, blocks = bl, geno = g, verbose = FALSE)
  expect_true(reg$core[1] <= ex$pos && reg$core[2] >= ex$pos)
})
