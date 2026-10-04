# Output regression: the simulated example data (fixed seed) must keep giving
# the same region, blocks and tables. A change here means that the results of
# the package changed and should be explained in NEWS.

test_that("the example data give the documented region and tables", {
  ex <- qp_example()
  on.exit(set_species(NULL, verbose = FALSE))
  set_species("rice", verbose = FALSE)
  expect_equal(ex$pos, 12559349L)
  g <- genotype_qc(read_vcf(ex$vcf, chr = 9, verbose = FALSE), verbose = FALSE)
  bl <- find_ld_blocks(g, verbose = FALSE)
  expect_equal(nrow(bl), 7L)
  expect_equal(bl$nsnps, c(43L, 16L, 33L, 61L, 14L, 21L, 45L))
  out <- file.path(tempdir(), "qp_regression")
  res <- quard_plot(ex$gwas, ex$gff, ex$vcf, chr = ex$chr, pos = ex$pos, panels = "none",
                    tables = c("snps", "ld"), output = file.path(out, "reg"), verbose = FALSE)
  expect_equal(c(res$region$core[1], res$region$core[2]), c(12514676, 12602793))
  expect_equal(res$lead, "S9_12559349")
  expect_equal(nrow(res$snp_table), 67L)
  expect_equal(sum(res$snp_table$R2_WITH_LEAD >= 0.8, na.rm = TRUE), 30L)
  expect_equal(as.vector(table(factor(res$snp_table$LOCATION,
                                      c("exon", "intron", "upstream", "downstream", "intergenic")))),
               c(11L, 7L, 17L, 4L, 28L))
  readme <- readLines(file.path(out, "reg_README.txt"))
  expect_true(any(grepl("^Made with: +quard_plot .* \\(R [0-9.]+\\)", readme)))
  expect_true(any(grepl("^Settings: +region = ld_blocks; min_maf = 0.05", readme)))
})
