# Compare the LD values and LD blocks of quardplot with PLINK 1.9
#
# Uses the same samples and variants: the region is cut from the VCF, PLINK
# filters it with --maf 0.05 --geno 0.2 and quardplot with genotype_qc()
# (min_maf = 0.05, max_missing = 0.2). Reports the number of variants kept,
# the maximum and mean absolute difference of r2 (genotype correlation,
# plink --r2), of haplotype r2 and D' (plink --r2 dprime), and how many LD
# blocks are identical (plink --blocks).
#
# Usage:
#   Rscript compare_plink.R file.vcf.gz CHR START END [path/to/plink]
# The two comparisons in the paper (Table 2):
#   Rscript compare_plink.R bench_data/bench_n500.vcf.gz Chr9 11456615 11556615   # inbred lines
#   Rscript compare_plink.R outbred.vcf.gz Chr9 12400000 12500000                 # outbred individuals
# (bench_data comes from benchmark_quardplot.R, outbred.vcf.gz from simulate_outbred.R.)
# Writes plink_comparison.csv and session_info.txt (R, package and PLINK versions).
# Needs: quardplot, PLINK 1.9 and tabix (htslib) on the PATH.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript compare_plink.R file.vcf.gz CHR START END [plink]")
vcf <- args[1]; chr <- args[2]; start <- as.numeric(args[3]); end <- as.numeric(args[4])
plink <- if (length(args) >= 5) args[5] else "plink"
suppressPackageStartupMessages(library(quardplot))
dir <- tempfile("plinkcmp")
dir.create(dir)
reg_vcf <- file.path(dir, "region.vcf")
system2("tabix", c("-h", vcf, sprintf("%s:%.0f-%.0f", chr, start, end)), stdout = reg_vcf)
run <- function(...) {
  st <- system2(plink, c(..., "--allow-extra-chr", "--silent"), stdout = FALSE, stderr = FALSE)
  if (!identical(st, 0L)) stop("PLINK failed: ", paste(..., collapse = " "))
}
b <- file.path(dir, "reg")
run("--vcf", reg_vcf, "--double-id", "--maf", "0.05", "--geno", "0.2", "--make-bed", "--out", b)
run("--bfile", b, "--r2", "square", "--out", paste0(b, "_sq"))
run("--bfile", b, "--r2", "dprime", "--ld-window", "100000", "--ld-window-kb", "100000",
    "--ld-window-r2", "0", "--out", paste0(b, "_tab"))
run("--bfile", b, "--blocks", "no-pheno-req", "--blocks-max-kb", "200", "--out", paste0(b, "_blk"))

g <- genotype_qc(read_vcf(reg_vcf, verbose = FALSE), min_maf = 0.05, max_missing = 0.2, verbose = FALSE)
bim <- utils::read.table(paste0(b, ".bim"))
r2_g <- compare_ld(calc_ld(g), read_plink_ld(paste0(b, "_sq.ld"), bim = paste0(b, ".bim")))
r2_h <- compare_ld(calc_ld(g, method = "haplotype"), read_plink_ld(paste0(b, "_tab.ld")))
dp <- compare_ld(calc_ld(g, stat = "dprime"), read_plink_ld(paste0(b, "_tab.ld"), stat = "DP"))
ours <- find_ld_blocks(g, max_kb = 200, verbose = FALSE)
theirs <- read_blocks(paste0(b, "_blk.blocks.det"))
key <- function(x) paste(x$start, x$end)
out <- data.frame(
  measure = c("variants kept (quardplot / PLINK)", "pairs compared",
              "r2 (genotype) max |difference|", "r2 (genotype) mean |difference|",
              "r2 (haplotype) max |difference|", "r2 (haplotype) mean |difference|",
              "D' max |difference|", "D' mean |difference|",
              "LD blocks (quardplot / PLINK / identical)"),
  value = c(sprintf("%d / %d", nrow(g$info), nrow(bim)), format(r2_g$n_pairs, big.mark = ","),
            format(r2_g$max_abs_diff, digits = 2), format(r2_g$mean_abs_diff, digits = 2),
            format(r2_h$max_abs_diff, digits = 2), format(r2_h$mean_abs_diff, digits = 2),
            format(dp$max_abs_diff, digits = 2), format(dp$mean_abs_diff, digits = 2),
            sprintf("%d / %d / %d", nrow(ours), nrow(theirs), sum(key(ours) %in% key(theirs)))))
print(out, right = FALSE, row.names = FALSE)
utils::write.csv(out, "plink_comparison.csv", row.names = FALSE)
pv <- tryCatch(system2(plink, "--version", stdout = TRUE, stderr = TRUE)[1], error = function(e) NA_character_)
writeLines(c(sprintf("Region: %s:%.0f-%.0f of %s", chr, start, end, basename(vcf)),
             paste("PLINK:", pv), "", utils::capture.output(utils::sessionInfo())), "session_info.txt")
