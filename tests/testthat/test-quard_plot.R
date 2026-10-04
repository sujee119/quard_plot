test_that("quard_plot writes the figure and the tables", {
  ex <- qp_example()
  on.exit(set_species(NULL, verbose = FALSE))
  set_species("rice", verbose = FALSE)
  out <- file.path(tempdir(), "qp_test_out")
  res <- quard_plot(ex$gwas, ex$gff, ex$vcf, chr = ex$chr, pos = ex$pos,
                    output = file.path(out, "locus.pdf"), verbose = FALSE)
  expect_s3_class(res, "quard_result")
  expect_true(all(file.exists(res$files)))
  expect_true(all(c("pdf", "genes", "snps", "ld_matrix", "readme") %in% names(res$files)))
  expect_true(all(c("SNP", "P", "R2_WITH_LEAD", "GENE", "LOCATION") %in% names(res$snp_table)))
  expect_true(ex$pos >= res$region$start && ex$pos <= res$region$end)
})

test_that("the LD table does not need the LD heatmap, and tables can be written alone", {
  ex <- qp_example()
  out <- file.path(tempdir(), "qp_test_out2")
  r1 <- quard_plot(ex$gwas, ex$gff, ex$vcf, chr = ex$chr, pos = ex$pos,
                   panels = c("manhattan", "regional", "genes"), output = file.path(out, "a.pdf"),
                   verbose = FALSE)
  expect_true(file.exists(file.path(out, "a_LD_r2.csv")))
  expect_null(r1$panels$ld)
  r2 <- quard_plot(ex$gwas, ex$gff, ex$vcf, chr = ex$chr, pos = ex$pos, panels = "none",
                   tables = c("snps", "ld_pairs"), output = file.path(out, "b"), verbose = FALSE)
  expect_false(file.exists(file.path(out, "b.pdf")))
  expect_true(file.exists(file.path(out, "b_LD_r2_pairs.csv")))
  skip_if_not_installed("writexl")
  r3 <- quard_plot(ex$gwas, ex$gff, ex$vcf, chr = ex$chr, pos = ex$pos, table_format = "xlsx",
                   output = file.path(out, "c.pdf"), verbose = FALSE)
  expect_true(file.exists(file.path(out, "c_tables.xlsx")))
})

test_that("clear errors for wrong options", {
  ex <- qp_example()
  expect_error(quard_plot(ex$gwas, output = "x.png", chr = 9, pos = ex$pos, verbose = FALSE), "PDF")
  expect_error(quard_plot(ex$gwas, panels = "heatmap", chr = 9, pos = ex$pos, verbose = FALSE), "Unknown panel")
})

test_that("check_inputs finds no problem in well-formed files", {
  ex <- qp_example()
  on.exit(set_species(NULL, verbose = FALSE))
  set_species("rice", verbose = FALSE)
  r <- check_inputs(list(gwas = ex$gwas, gff = ex$gff, vcf = ex$vcf, hapmap = ex$complete_hapmap),
                    verbose = FALSE)
  expect_false(any(r$status == "PROBLEM"))
})
