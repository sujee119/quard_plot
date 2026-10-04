test_that("SNPs are placed in exons, introns, promoters and between genes", {
  f <- qp_hand_gff()
  pos <- c(10500, 11500, 13000, 19500, 8000, 6999, 20500, 30200, 29500, 41000, 39800, 61000)
  a <- annotate_variants(pos, f, chr = "chr1")
  expect_equal(a$location, c("exon", "exon", "intron", "exon", "upstream", "intergenic", "downstream",
                             "exon", "downstream", "upstream", "exon", "genic"))
  expect_equal(a$feature[c(1, 2, 4, 8, 11)], c("5'UTR", "CDS", "3'UTR", "3'UTR", "5'UTR"))
  expect_equal(a$distance_bp[c(5, 6, 7, 9, 10)], c(2000, 3001, 500, 500, 1000))
  b <- annotate_variants(8000, f, chr = "chr1", upstream_bp = 1000)
  expect_equal(b$location, "intergenic")
})

test_that("a GFF with gene lines only can be read", {
  f <- tempfile(fileext = ".gff3")
  writeLines(c("##gff-version 3", "Chr1\tx\tgene\t100\t900\t.\t+\t.\tID=g1"), f)
  g <- read_gff(f, verbose = FALSE)
  expect_equal(nrow(g$genes), 1L)
  expect_equal(annotate_variants(500, g, chr = 1)$location, "genic")
})
