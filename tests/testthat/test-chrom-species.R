test_that("chromosome names are matched across naming styles", {
  expect_equal(normalize_chr(c("Chr1", "chr01", "1", "CHR_X", "ChrM", "ChrC", "SL4.0ch01", "ch12")),
               c("1", "1", "1", "X", "MT", "Pt", "1", "12"))
  expect_equal(sort_chr(c("10", "2", "X", "1", "MT")), c("1", "2", "10", "X", "MT"))
})

test_that("built-in chromosome lengths and the species setting work", {
  on.exit(set_species(NULL, verbose = FALSE))
  cs <- chrom_sizes("rice")
  expect_equal(nrow(cs), 12L)
  expect_equal(cs$length[cs$chr == "9"], 23012720)
  expect_true(all(c("rice", "arabidopsis", "tomato", "human", "mouse") %in% list_species()$species))
  expect_message(set_species("tomato"), "Tomato")
  expect_equal(get_species()$name, "tomato")
  expect_equal(nrow(chrom_sizes()), 12L)
  set_species("my_plant", chrom_sizes = c(Chr1 = 1e6, Chr2 = 2e6), verbose = FALSE)
  expect_equal(chrom_sizes()$length, c(1e6, 2e6))
  set_species(NULL, verbose = FALSE)
  expect_null(get_species())
})

test_that("chromosome lengths are read from a FASTA file", {
  fa <- tempfile(fileext = ".fa")
  writeLines(c(">chr1 test", strrep("A", 60), strrep("C", 40), ">chr2", strrep("G", 25)), fa)
  cs <- chrom_sizes(fa)
  expect_equal(cs$length, c(100, 25))
})
