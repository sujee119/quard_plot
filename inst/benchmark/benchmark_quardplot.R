# Benchmark of quard_plot(): run time and peak memory
#
# Simulates a rice-like data set (GWAS results with 1,000,000 SNPs on 12
# chromosomes; genotypes for 100,000 variants on chromosome 9 with
# haplotype-block structure; gene models for chromosome 9), then runs
# quard_plot() in a fresh R process for every setting and records the
# elapsed time, the peak memory (resident set size) and the size of the PDF.
#
#   1. number of samples: 100, 500 and 1,000, with a bgzip-compressed and
#      tabix-indexed VCF and with a plain (uncompressed, unindexed) VCF;
#   2. number of variants in the plotted region (500 samples, indexed VCF):
#      windows of 50, 200, 500 and 1,000 kb.
#
# Usage (from a terminal, in an empty folder):
#   Rscript benchmark_quardplot.R            # all settings, 3 repetitions
#   Rscript benchmark_quardplot.R 5          # 5 repetitions
# Needs: quardplot, and bgzip + tabix (htslib) on the PATH for the indexed
# runs. Peak memory is read from /proc (Linux); it is NA on other systems.
# Results: bench_results.csv (one row per run), bench_summary.csv (median
# per setting) and bench_system.txt (hardware, R and package versions).

args <- commandArgs(trailingOnly = TRUE)
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])

peak_rss_mb <- function() {
  f <- "/proc/self/status"
  if (!file.exists(f)) return(NA_real_)
  x <- grep("^VmHWM:", readLines(f), value = TRUE)
  if (!length(x)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", x)) / 1024
}

# ---- one run in a fresh process -------------------------------------------------
if (length(args) && args[1] == "one") {
  suppressPackageStartupMessages(library(quardplot))
  vcf <- args[2]; gwas <- args[3]; gff <- args[4]; out <- args[5]; win <- as.numeric(args[6])
  pos <- as.numeric(args[7])
  set_species("rice", verbose = FALSE)
  t0 <- proc.time()[["elapsed"]]
  res <- if (win > 0) {
    quard_plot(gwas, gff, vcf, chr = 9, pos = pos, start = pos - win / 2, end = pos + win / 2,
               output = out, verbose = FALSE)
  } else {
    quard_plot(gwas, gff, vcf, chr = 9, pos = pos, output = out, verbose = FALSE)
  }
  el <- proc.time()[["elapsed"]] - t0
  n_loc <- if (is.null(res$genotypes)) 0L else nrow(res$genotypes$info)
  cat(sprintf("RESULT,%.3f,%.1f,%d,%d,%.0f\n", el, peak_rss_mb(), n_loc,
              if (is.null(res$ld)) 0L else nrow(res$ld), file.size(out) / 1024))
  quit(save = "no")
}

# ---- simulated data --------------------------------------------------------------
make_data <- function(dir, n_samples, n_var = 100000L, seed = 1L) {
  set.seed(seed)
  len9 <- 23012720
  pos <- sort(sample.int(len9 - 1000L, n_var)) + 500L
  # haplotype blocks: 6 founder haplotypes per block, inbred lines,
  # 2% heterozygous calls and 2% missing calls
  bsize <- pmax(5L, round(stats::rexp(n_var / 25, 1 / 40)))
  blk <- rep(seq_along(bsize), bsize)[seq_len(n_var)]
  smp <- sprintf("Acc%04d", seq_len(n_samples))
  vcf <- file.path(dir, sprintf("bench_n%d.vcf", n_samples))
  con <- file(vcf, "w")
  writeLines(c("##fileformat=VCFv4.2", "##source=quardplot benchmark (simulated)",
               sprintf("##contig=<ID=Chr9,length=%d>", len9),
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", smp),
                     collapse = "\t")), con)
  y <- NULL
  causal <- which.min(abs(pos - len9 / 2))
  pvals <- numeric(n_var)
  gt_lab <- c("0/0", "0/1", "1/1")
  bases <- c("A", "C", "G", "T")
  chunks <- split(seq_len(n_var), ceiling(seq_len(n_var) / 4000L))
  # founder haplotypes and founder choice per block
  ub <- unique(blk)
  found <- lapply(ub, function(b) {
    m <- sum(blk == b)
    list(h = matrix(stats::runif(6 * m) < stats::runif(m, 0.1, 0.9)[rep(seq_len(m), each = 6)], 6, m),
         who = sample.int(6L, n_samples, replace = TRUE))
  })
  names(found) <- ub
  bstart <- match(ub, blk)
  geno_of <- function(idx) {
    G <- matrix(0L, length(idx), n_samples)
    for (b in unique(blk[idx])) {
      k <- idx[blk[idx] == b]
      first <- bstart[b]
      cols <- k - first + 1L
      f <- found[[as.character(b)]]
      G[match(k, idx), ] <- t(2L * f$h[f$who, cols, drop = FALSE])
    }
    het <- matrix(stats::runif(length(G)) < 0.02, nrow(G))
    G[het] <- 1L
    G[matrix(stats::runif(length(G)) < 0.02, nrow(G))] <- NA_integer_
    G
  }
  gc_ <- geno_of(causal)[1, ]
  gc_[is.na(gc_)] <- 0L
  y <- 0.9 * gc_ / 2 + stats::rnorm(n_samples)
  for (idx in chunks) {
    G <- geno_of(idx)
    # association P values (correlation test with the simulated trait)
    Gm <- G
    mu <- rowMeans(Gm, na.rm = TRUE)
    Gm[is.na(Gm)] <- matrix(mu, nrow(Gm), ncol(Gm))[is.na(Gm)]
    yc <- y - mean(y)
    Gc <- Gm - mu
    r <- as.vector(Gc %*% yc) / sqrt(rowSums(Gc^2) * sum(yc^2))
    r[!is.finite(r)] <- 0
    tt <- r * sqrt((n_samples - 2) / pmax(1e-12, 1 - r^2))
    pvals[idx] <- 2 * stats::pt(-abs(tt), n_samples - 2)
    ref <- sample(bases, length(idx), replace = TRUE)
    alt <- vapply(ref, function(z) sample(setdiff(bases, z), 1), "")
    lab <- matrix(gt_lab[G + 1L], nrow(G))
    lab[is.na(G)] <- "./."
    fixed <- paste("Chr9", pos[idx], sprintf("S9_%d", pos[idx]), ref, alt, ".", "PASS", ".", "GT", sep = "\t")
    writeLines(do.call(paste, c(list(fixed), as.data.frame(lab, stringsAsFactors = FALSE), sep = "\t")), con)
  }
  close(con)
  # GWAS results: chromosome 9 from the genotypes + 900,000 background SNPs
  sizes <- chrom_sizes("rice")
  other <- sizes[sizes$chr != "9", ]
  nb <- 1000000L - n_var
  chr_b <- sample(other$chr, nb, replace = TRUE, prob = other$length)
  pos_b <- round(stats::runif(nb, 1, other$length[match(chr_b, other$chr)]))
  gw <- data.frame(SNP = c(sprintf("S%s_%d", chr_b, pos_b), sprintf("S9_%d", pos)),
                   Chr = c(chr_b, rep("9", n_var)), Pos = c(pos_b, pos),
                   P.value = signif(c(stats::runif(nb), pvals), 4))
  gw <- gw[!duplicated(paste(gw$Chr, gw$Pos)), ]
  f_gw <- file.path(dir, sprintf("bench_gwas_n%d.csv", n_samples))
  data.table::fwrite(gw, f_gw)
  list(vcf = vcf, gwas = f_gw, causal = pos[causal])
}

# ---- driver -------------------------------------------------------------------------
reps <- if (length(args)) as.integer(args[1]) else 3L
dir <- "bench_data"
dir.create(dir, showWarnings = FALSE)
dir.create("bench_out", showWarnings = FALSE)
suppressPackageStartupMessages(library(quardplot))
f_gff <- file.path(dir, "bench_genes.gff3")
set.seed(9)
writeLines(quardplot:::.sim_gff_lines("Chr9", 1, 23012720), f_gff)
has_tabix <- nzchar(Sys.which("bgzip")) && nzchar(Sys.which("tabix"))
sets <- list()
for (n in c(100L, 500L, 1000L)) {
  t0 <- proc.time()[["elapsed"]]
  d <- make_data(dir, n)
  message(sprintf("Simulated %d samples in %.0f s", n, proc.time()[["elapsed"]] - t0))
  sets[[as.character(n)]] <- d
  if (has_tabix) {
    system2("bgzip", c("-f", "-k", d$vcf))
    system2("tabix", c("-f", "-p", "vcf", paste0(d$vcf, ".gz")))
  }
}
runs <- list()
add_run <- function(setting, n, vcf, win) {
  d <- sets[[as.character(n)]]
  for (r in seq_len(reps)) {
    out <- file.path("bench_out", sprintf("%s_rep%d.pdf", gsub("[^A-Za-z0-9]+", "_", setting), r))
    o <- system2(file.path(R.home("bin"), "Rscript"),
                 c(this_file, "one", vcf, d$gwas, f_gff, out, win, d$causal), stdout = TRUE, stderr = FALSE)
    v <- strsplit(sub("^RESULT,", "", grep("^RESULT,", o, value = TRUE)), ",")[[1]]
    runs[[length(runs) + 1L]] <<- data.frame(setting = setting, samples = n,
                                             vcf = if (grepl("\\.gz$", vcf)) "bgzip + tabix" else "plain",
                                             window_kb = if (win > 0) win / 1000 else NA, rep = r,
                                             seconds = as.numeric(v[1]), peak_mb = as.numeric(v[2]),
                                             region_variants = as.integer(v[3]),
                                             heatmap_variants = as.integer(v[4]), pdf_kb = as.numeric(v[5]))
    message(sprintf("%-28s rep %d: %6.1f s, %6.0f MB", setting, r, as.numeric(v[1]), as.numeric(v[2])))
  }
}
for (n in c(100L, 500L, 1000L)) {
  d <- sets[[as.character(n)]]
  if (has_tabix) add_run(sprintf("samples %d, indexed", n), n, paste0(d$vcf, ".gz"), 0)
  add_run(sprintf("samples %d, plain", n), n, d$vcf, 0)
}
if (has_tabix) {
  for (w in c(50e3, 200e3, 500e3, 1000e3)) {
    add_run(sprintf("window %d kb", w / 1000), 500L, paste0(sets[["500"]]$vcf, ".gz"), w)
  }
}
res <- do.call(rbind, runs)
utils::write.csv(res, "bench_results.csv", row.names = FALSE)
agg <- stats::aggregate(cbind(seconds, peak_mb, region_variants, heatmap_variants, pdf_kb) ~
                          setting + samples + vcf, data = res, FUN = stats::median)
utils::write.csv(agg, "bench_summary.csv", row.names = FALSE)
info <- c(paste("Date:", format(Sys.time())),
          paste("R:", R.version.string), paste("Platform:", R.version$platform),
          paste("OS:", utils::sessionInfo()$running),
          paste("CPU:", tryCatch(trimws(sub(".*:", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1])),
                                 error = function(e) NA)),
          paste("Cores:", parallel::detectCores()),
          paste("Memory (GB):", tryCatch(round(as.numeric(gsub("[^0-9]", "", grep("MemTotal", readLines("/proc/meminfo"),
                                                                                value = TRUE))) / 1024^2, 1),
                                         error = function(e) NA)),
          paste("Packages:", paste(sprintf("%s %s", c("quardplot", "data.table", "ggplot2", "patchwork"),
                                         vapply(c("quardplot", "data.table", "ggplot2", "patchwork"),
                                                function(p) as.character(utils::packageVersion(p)), "")),
                                 collapse = ", ")),
          paste("VCF sizes (MB):", paste(sprintf("%d samples: %.0f plain / %.0f bgzip", c(100, 500, 1000),
                                                 vapply(sets, function(d) file.size(d$vcf) / 1e6, 0),
                                                 vapply(sets, function(d) {
                                                   g <- paste0(d$vcf, ".gz")
                                                   if (file.exists(g)) file.size(g) / 1e6 else NA
                                                 }, 0)), collapse = "; ")))
writeLines(info, "bench_system.txt")
print(agg)
cat(info, sep = "\n")
