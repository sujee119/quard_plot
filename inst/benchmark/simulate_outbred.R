# Simulated outbred diploid population for the PLINK comparison
#
# Writes outbred.vcf.gz (bgzip + tabix): 300 individuals, 2,000 variants on
# 1 Mb of chromosome 9 with haplotype-block structure. Each individual
# carries two haplotypes drawn from 8 founder haplotypes per block, so many
# genotypes are heterozygous and phase is unknown; 2% of calls are missing.
#
# Usage: Rscript simulate_outbred.R [n_individuals] [n_variants]
# Then:  Rscript compare_plink.R outbred.vcf.gz Chr9 START END

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 300L
m <- if (length(args) >= 2) as.integer(args[2]) else 2000L
set.seed(7)
pos <- sort(sample(12000000:13000000, m))
bsize <- pmax(5L, round(stats::rexp(m / 20, 1 / 30)))
blk <- rep(seq_along(bsize), bsize)[seq_len(m)]
G <- matrix(0L, m, n)
for (b in unique(blk)) {
  k <- which(blk == b)
  f <- matrix(stats::runif(8 * length(k)) < rep(stats::runif(length(k), 0.1, 0.9), each = 8), 8)
  h1 <- sample.int(8L, n, replace = TRUE)
  h2 <- sample.int(8L, n, replace = TRUE)
  G[k, ] <- t(f[h1, , drop = FALSE] + f[h2, , drop = FALSE])
}
G[matrix(stats::runif(length(G)) < 0.02, m)] <- NA
lab <- matrix(c("0/0", "0/1", "1/1")[G + 1L], m)
lab[is.na(G)] <- "./."
bases <- c("A", "C", "G", "T")
ref <- sample(bases, m, replace = TRUE)
alt <- vapply(ref, function(z) sample(setdiff(bases, z), 1), "")
fixed <- paste("Chr9", pos, sprintf("S9_%d", pos), ref, alt, ".", "PASS", ".", "GT", sep = "\t")
body <- do.call(paste, c(list(fixed), as.data.frame(lab, stringsAsFactors = FALSE), sep = "\t"))
writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=Chr9,length=23012720>",
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                     sprintf("Ind%03d", seq_len(n))), collapse = "\t"), body), "outbred.vcf")
system2("bgzip", c("-f", "outbred.vcf"))
system2("tabix", c("-f", "-p", "vcf", "outbred.vcf.gz"))
cat("outbred.vcf.gz:", n, "individuals,", m, "variants, Chr9:", min(pos), "-", max(pos), "\n")
