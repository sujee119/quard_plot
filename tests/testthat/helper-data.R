# Small simulated data set, made once per test run.
qp_example <- local({
  ex <- NULL
  function() {
    if (is.null(ex)) ex <<- suppressMessages(simulate_example_data(file.path(tempdir(), "qp_test_data")))
    ex
  }
})

# A small hand-made GFF3 with known gene structures.
qp_hand_gff <- function() {
  g <- function(type, s, e, str, attr) paste("chr1", "test", type, s, e, ".", str, ".", attr, sep = "\t")
  L <- c("##gff-version 3",
         g("gene", 10000, 20000, "+", "ID=geneA;Note=alpha"), g("mRNA", 10000, 20000, "+", "ID=geneA.1;Parent=geneA"),
         g("exon", 10000, 12000, "+", "Parent=geneA.1"), g("five_prime_UTR", 10000, 10999, "+", "Parent=geneA.1"),
         g("CDS", 11000, 12000, "+", "Parent=geneA.1"), g("exon", 14000, 15000, "+", "Parent=geneA.1"),
         g("CDS", 14000, 15000, "+", "Parent=geneA.1"), g("exon", 18000, 20000, "+", "Parent=geneA.1"),
         g("CDS", 18000, 19000, "+", "Parent=geneA.1"), g("three_prime_UTR", 19001, 20000, "+", "Parent=geneA.1"),
         g("gene", 30000, 40000, "-", "ID=geneB"), g("mRNA", 30000, 40000, "-", "ID=geneB.1;Parent=geneB"),
         g("exon", 30000, 31000, "-", "Parent=geneB.1"), g("three_prime_UTR", 30000, 30499, "-", "Parent=geneB.1"),
         g("CDS", 30500, 31000, "-", "Parent=geneB.1"), g("exon", 35000, 36000, "-", "Parent=geneB.1"),
         g("CDS", 35000, 36000, "-", "Parent=geneB.1"), g("exon", 39000, 40000, "-", "Parent=geneB.1"),
         g("CDS", 39000, 39499, "-", "Parent=geneB.1"), g("five_prime_UTR", 39500, 40000, "-", "Parent=geneB.1"),
         g("gene", 60000, 62000, "+", "ID=geneC"))
  f <- tempfile(fileext = ".gff3")
  writeLines(L, f)
  f
}
