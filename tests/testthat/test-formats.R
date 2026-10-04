# Reading the input formats listed in Table 1 of the paper.

test_that("GTF files give the same genes as the equivalent GFF3", {
  g <- function(type, s, e, str, attr) paste("chr1", "test", type, s, e, ".", str, ".", attr, sep = "\t")
  gtf <- c(g("gene", 10000, 20000, "+", 'gene_id "geneA"; gene_name "alpha";'),
           g("transcript", 10000, 20000, "+", 'gene_id "geneA"; transcript_id "geneA.1";'),
           g("exon", 10000, 12000, "+", 'gene_id "geneA"; transcript_id "geneA.1";'),
           g("CDS", 11000, 12000, "+", 'gene_id "geneA"; transcript_id "geneA.1";'),
           g("exon", 18000, 20000, "+", 'gene_id "geneA"; transcript_id "geneA.1";'),
           g("CDS", 18000, 19000, "+", 'gene_id "geneA"; transcript_id "geneA.1";'),
           g("gene", 30000, 40000, "-", 'gene_id "geneB";'),
           g("transcript", 30000, 40000, "-", 'gene_id "geneB"; transcript_id "geneB.1";'),
           g("exon", 30000, 31000, "-", 'gene_id "geneB"; transcript_id "geneB.1";'),
           g("exon", 39000, 40000, "-", 'gene_id "geneB"; transcript_id "geneB.1";'))
  f <- tempfile(fileext = ".gtf")
  writeLines(gtf, f)
  gn <- read_gff(f, verbose = FALSE)
  expect_s3_class(gn, "qp_genes")
  expect_setequal(gn$genes$gene_id, c("geneA", "geneB"))
  expect_equal(nrow(gn$transcripts), 2L)
  a <- annotate_variants(data.frame(chr = "chr1", pos = c(11500, 15000, 39500)), gn)
  expect_equal(a$gene_id, c("geneA", "geneA", "geneB"))
  expect_equal(a$location, c("exon", "intron", "exon"))
  expect_equal(a$feature[1], "CDS")
})

test_that("block files from PLINK 1.9 and simple CHR/START/END tables are read", {
  det <- tempfile(fileext = ".blocks.det")
  writeLines(c(" CHR          BP1          BP2           KB  NSNPS SNPS",
               "   9     12526365     12531036        4.672      4 S9_12526365|S9_12526435|S9_12528167|S9_12531036",
               "   9     12540000     12550000       10.001      2 S9_12540000|S9_12550000"), det)
  b <- read_blocks(det)
  expect_equal(nrow(b), 2L)
  expect_equal(b$start, c(12526365, 12540000))
  expect_equal(b$end, c(12531036, 12550000))
  expect_equal(b$nsnps, c(4L, 2L))
  tab <- tempfile(fileext = ".txt")
  writeLines(c("CHR START END", "Chr9 100 900", "chr09 1000 2000"), tab)
  b2 <- read_blocks(tab)
  expect_equal(b2$chr, c("9", "9"))
  expect_equal(b2$end, c(900, 2000))
})

test_that("gzip-compressed PLINK results keep only the ADD rows", {
  f <- tempfile(fileext = ".assoc.linear.gz")
  con <- gzfile(f, "w")
  writeLines(c(" CHR        SNP         BP   A1       TEST    NMISS       BETA         STAT            P",
               "   1       rs1       1000    A        ADD      100      0.10        1.0        0.30",
               "   1       rs1       1000    A       COV1      100      2.00        5.0       1e-06",
               "   1       rs2       2000    G        ADD      100      0.50        4.0       1e-04",
               "   1       rs2       2000    G       COV1      100      2.00        5.0       1e-06"), con)
  close(con)
  gw <- read_gwas(f, verbose = FALSE)
  expect_equal(gw$snp, c("rs1", "rs2"))
  expect_equal(gw$p, c(0.30, 1e-04))
})

test_that("GEMMA and rMVP column names are recognised", {
  gemma <- data.frame(chr = c(1, 1), rs = c("a", "b"), ps = c(10, 20), n_miss = 0, allele1 = "A",
                      allele0 = "G", af = 0.3, beta = 0.1, se = 0.05, p_wald = c(0.01, 0.5))
  gw <- read_gwas(gemma, verbose = FALSE)
  expect_equal(gw$snp, c("a", "b"))
  expect_equal(gw$p, c(0.01, 0.5))
  rmvp <- data.frame(SNP = c("s1", "s2"), CHROM = c("Chr2", "Chr2"), POS = c(5, 50), REF = "A", ALT = "T",
                     Effect = 0.2, SE = 0.1, Height.MLM = c(1e-6, 0.2), check.names = FALSE)
  gw2 <- read_gwas(rmvp, verbose = FALSE)
  expect_equal(gw2$p, c(1e-6, 0.2))
  expect_equal(gw2$chr, c("2", "2"))
})

test_that("LD tables written by PLINK are read into a matrix", {
  f <- tempfile(fileext = ".ld")
  writeLines(c(" CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2",
               "     1          100           v1      1          200           v2         0.81",
               "     1          100           v1      1          300           v3         0.04",
               "     1          200           v2      1          300           v3          0.1"), f)
  m <- read_plink_ld(f)
  expect_equal(dim(m), c(3L, 3L))
  expect_equal(unname(m["v1", "v2"]), 0.81)
  expect_equal(unname(m["v3", "v2"]), 0.1)
  expect_equal(unname(diag(m)), c(1, 1, 1))
})

test_that("a tabix-indexed VCF gives the same genotypes as reading the whole file", {
  skip_if(!nzchar(Sys.which("tabix")) || !nzchar(Sys.which("bgzip")), "tabix/bgzip not installed")
  ex <- qp_example()
  plain <- tempfile(fileext = ".vcf")
  writeLines(readLines(gzfile(ex$vcf)), plain)
  system2("bgzip", c("-f", plain))
  gz <- paste0(plain, ".gz")
  system2("tabix", c("-f", "-p", "vcf", gz))
  a <- read_vcf(gz, chr = 9, start = 12500000, end = 12600000, method = "tabix", verbose = FALSE)
  b <- read_vcf(gz, chr = 9, start = 12500000, end = 12600000, method = "stream", verbose = FALSE)
  expect_equal(a$info$pos, b$info$pos)
  expect_equal(a$dosage, b$dosage)
})
