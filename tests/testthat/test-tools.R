test_that("tables are written as CSV (and Excel when writexl is installed)", {
  d <- data.frame(gene = c("A", "B"), length = c(1200, 3400))
  f <- tempfile(fileext = ".csv")
  write_table(d, f)
  expect_equal(utils::read.csv(f), d)
  skip_if_not_installed("writexl")
  x <- tempfile(fileext = ".xlsx")
  write_table(list(one = d, two = d), x)
  expect_true(file.exists(x))
})

test_that("HapMap to VCF keeps the genotypes", {
  ex <- qp_example()
  v <- tempfile(fileext = ".vcf")
  hapmap_to_vcf(ex$hapmap, v, verbose = FALSE)
  a <- read_hapmap(ex$hapmap, verbose = FALSE)
  b <- read_vcf(v, verbose = FALSE)
  common <- intersect(a$info$pos, b$info$pos)
  expect_gt(length(common), 100L)
  da <- a$dosage[match(common, a$info$pos), ]
  db <- b$dosage[match(common, b$info$pos), ]
  same <- (da == db) | (da == 2 - db)        # REF/ALT may be swapped without a reference genome
  expect_true(all(same | is.na(da) | is.na(db)))
})

test_that("gene lists and gene structure plots work", {
  ex <- qp_example()
  gg <- get_genes(ex$gff, chr = 9, start = 12500000, end = 12600000, verbose = FALSE)
  expect_gt(nrow(gg), 0L)
  expect_true(all(gg$end >= 12500000 & gg$start <= 12600000))
  p <- plot_gene_structure(ex$gff, gg$gene_id[1])
  expect_s3_class(p, "ggplot")
})

test_that("tree PDFs show the sample names", {
  ex <- qp_example()
  out <- file.path(tempdir(), "qp_tree")
  tr <- hapmap_tree(ex$genome_hapmap, file.path(out, "tree.pdf"), k = 3, verbose = FALSE)
  expect_true(all(file.exists(tr$files)))
  has_names <- function(p) {
    any(vapply(p$layers, function(l) inherits(l$geom, "GeomText") && is.data.frame(l$data) &&
                 "lab" %in% names(l$data) && all(tr$tree$tip.label %in% l$data$lab), NA))
  }
  expect_true(has_names(tr$plot))
  expect_length(attr(tr$plot, "size_mm"), 2L)
  expect_true(has_names(plot_tree(tr$tree, layout = "rectangular")))
  expect_true(has_names(plot_tree(tr$tree, layout = "unrooted")))
  expect_false(has_names(plot_tree(tr$tree, layout = "rectangular", labels = FALSE)))
  p3 <- plot_tree(tr$tree, size_mm = c(150, 170))
  expect_lte(attr(p3, "label_size"), 8)
})

test_that("karyotype plot counts SNPs and indels", {
  ex <- qp_example()
  k <- plot_karyotype(ex$genome_hapmap, chrom_sizes = "rice", bin_size = 2e6, verbose = FALSE)
  d <- attr(k, "density")
  expect_true(all(c("snp", "indel") %in% names(d)))
  expect_gt(sum(d$snp), 0)
  expect_gt(sum(d$indel), 0)
})
