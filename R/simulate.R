# Small simulated example data ---------------------------------------------------

# Simulated GFF3 gene models for a region (genes, mRNAs, exons, CDS, UTRs).
.sim_gff_lines <- function(seqid, r0, r1, prefix = "Sim09g") {
  out <- "##gff-version 3"
  notes <- c("expressed protein", "protein kinase domain containing protein",
             "MYB family transcription factor", "transporter family protein",
             "hypothetical protein", "NB-ARC domain containing protein")
  s <- r0 + 3000
  k <- 0L
  while (s < r1 - 8000) {
    k <- k + 1L
    e <- s + sample(1500:7000, 1)
    strand <- sample(c("+", "-"), 1)
    gid <- sprintf("%s%05d", prefix, 12400L + k * 10L)
    note <- sample(notes, 1)
    attr <- if (k %% 4L == 0L) {
      sprintf("ID=%s;Note=%s;Name=%s", gid, note, gid)
    } else {
      sprintf("ID=%s;Name=%s;Note=%s", gid, gid, note)
    }
    out <- c(out, paste(seqid, "sim", "gene", s, e, ".", strand, ".", attr, sep = "\t"))
    nex <- sample(2:6, 1)
    br <- sort(sample((s + 100):(e - 100), 2L * (nex - 1L)))
    ex_s <- c(s, br[seq(2L, length(br), 2L)])
    ex_e <- c(br[seq(1L, length(br), 2L)], e)
    ntx <- if (stats::runif(1) < 0.3 && nex > 2L) 2L else 1L
    for (t in seq_len(ntx)) {
      tid <- sprintf("%s.%d", gid, t)
      es <- ex_s
      ee <- ex_e
      if (t == 2L) {
        drop <- sample(2:(nex - 1L), 1)
        es <- es[-drop]
        ee <- ee[-drop]
      }
      out <- c(out, paste(seqid, "sim", "mRNA", min(es), max(ee), ".", strand, ".",
                          sprintf("ID=%s;Parent=%s;Name=%s", tid, gid, tid), sep = "\t"))
      u5 <- sample(80:250, 1)
      u3 <- sample(120:350, 1)
      lo <- if (strand == "+") min(es) + u5 else min(es) + u3
      hi <- if (strand == "+") max(ee) - u3 else max(ee) - u5
      for (x in seq_along(es)) {
        out <- c(out, paste(seqid, "sim", "exon", es[x], ee[x], ".", strand, ".",
                            sprintf("ID=%s.exon%d;Parent=%s", tid, x, tid), sep = "\t"))
      }
      cs <- pmax(es, lo)
      ce <- pmin(ee, hi)
      ok <- cs <= ce
      cs <- cs[ok]
      ce <- ce[ok]
      ord <- if (strand == "+") seq_along(cs) else rev(seq_along(cs))
      len_before <- c(0, cumsum((ce - cs + 1)[ord]))[seq_along(ord)]
      phase <- numeric(length(cs))
      phase[ord] <- (3 - (len_before %% 3)) %% 3
      for (x in seq_along(cs)) {
        out <- c(out, paste(seqid, "sim", "CDS", cs[x], ce[x], ".", strand, phase[x],
                            sprintf("ID=%s.cds;Parent=%s", tid, tid), sep = "\t"))
      }
      left_type <- if (strand == "+") "five_prime_UTR" else "three_prime_UTR"
      right_type <- if (strand == "+") "three_prime_UTR" else "five_prime_UTR"
      for (x in seq_along(es)) {
        if (es[x] < lo) {
          out <- c(out, paste(seqid, "sim", left_type, es[x], min(ee[x], lo - 1), ".", strand, ".",
                              sprintf("Parent=%s", tid), sep = "\t"))
        }
        if (ee[x] > hi) {
          out <- c(out, paste(seqid, "sim", right_type, max(es[x], hi + 1), ee[x], ".", strand, ".",
                              sprintf("Parent=%s", tid), sep = "\t"))
        }
      }
    }
    s <- e + sample(2500:15000, 1)
  }
  out
}

#' Write a small simulated example data set
#'
#' Creates a complete, internally consistent test data set (a few hundred
#' kilobytes) so the functions can be tried without large real files. All
#' data are simulated: a rice-like genome (IRGSP-1.0 chromosome lengths),
#' genome-wide GWAS P values with three association peaks, and for a 300-kb
#' region on chromosome 9: genotypes of inbred accessions with haplotype-block
#' structure (VCF with `GT:DP` fields and the same data in HapMap format), a
#' phenotype with one causal variant, GWAS P values computed from these
#' genotypes, GFF3 gene models and the LD blocks found by [find_ld_blocks()]
#' (PLINK `.blocks.det` format). The VCF uses `Chr9`, the GWAS `9` and the
#' GFF3 `Chr9`, to show that chromosome names are matched automatically.
#'
#' @param dir Output directory.
#' @param n_samples Number of accessions.
#' @param seed Random seed.
#' @return Invisibly, a named list of file paths (`gwas`, `gff`, `vcf`,
#'   `hapmap`, `phenotype`, `blocks`, `genome_hapmap` = whole-genome HapMap
#'   with SNPs and indels from four simulated subpopulations,
#'   `complete_hapmap` = the whole-genome HapMap plus the dense chr 9 region,
#'   `groups` = subpopulation of each accession) plus the target `chr` and
#'   `pos` (the causal variant).
#' @examples
#' \donttest{
#' ex <- simulate_example_data(tempfile("qp_example"))
#' str(ex)
#' }
#' @export
simulate_example_data <- function(dir = "quardplot_example", n_samples = 120L, seed = 1L) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  sizes <- chrom_sizes("rice")
  chr_t <- "9"
  r0 <- 12400000
  r1 <- 12700000
  n_snp <- 240L
  pos <- sort(sample(r0:r1, n_snp))
  n_blk <- 7L
  cuts <- sort(sample(12:(n_snp - 12), n_blk - 1L))
  blk <- findInterval(seq_len(n_snp), c(1L, cuts))
  n_f <- 6L
  founders <- matrix(0L, n_f, n_snp)
  for (b in seq_len(n_blk)) {
    # within a block the founders that carry an allele are nested sets (no
    # historical recombination: D' = 1), so each block is one LD block
    perm <- sample(n_f)
    parts <- lapply(2:4, function(z) perm[seq_len(z)])
    for (k in which(blk == b)) {
      u <- stats::runif(1)
      carriers <- if (u < 0.45) parts[[1]] else if (u < 0.7) parts[[2]] else if (u < 0.85) parts[[3]] else
        perm[seq_len(sample(1:(n_f - 1L), 1))]
      founders[carriers, k] <- 1L
    }
  }
  hap <- matrix(0L, n_samples, n_snp)
  for (s in seq_len(n_samples)) {
    f <- sample(n_f, 1)
    for (b in seq_len(n_blk)) {
      if (b > 1L && stats::runif(1) < 0.45) f <- sample(n_f, 1)
      idx <- which(blk == b)
      hap[s, idx] <- founders[f, idx]
    }
  }
  nn <- n_samples * n_snp
  flip <- matrix(stats::runif(nn) < 0.005, n_samples)
  hap[flip] <- 1L - hap[flip]
  geno <- 2L * hap
  geno[matrix(stats::runif(nn) < 0.01, n_samples)] <- 1L
  geno[matrix(stats::runif(nn) < 0.02, n_samples)] <- NA
  rare <- sample(seq_len(n_snp), 6)
  for (k in rare) {
    geno[, k] <- 0L
    geno[sample(n_samples, 2), k] <- 2L
  }
  ids_reg <- sprintf("S%s_%d", chr_t, pos)
  samples <- sprintf("Acc%03d", seq_len(n_samples))

  # genotypes: VCF (GT:DP, one multi-allelic site) and the same calls as HapMap
  bases <- c("A", "C", "G", "T")
  ref <- sample(bases, n_snp, replace = TRUE)
  alt <- vapply(ref, function(b) sample(setdiff(bases, b), 1), "")
  multi <- setdiff(seq_len(n_snp), rare)[3]
  gt <- ifelse(is.na(geno), "./.", ifelse(geno == 0L, "0/0", ifelse(geno == 1L, "0/1", "1/1")))
  dp <- matrix(sample(5:30, nn, replace = TRUE), n_samples)
  dp[is.na(geno)] <- 0L
  calls <- matrix(paste0(gt, ":", dp), n_samples)
  alt_vcf <- alt
  alt_vcf[multi] <- paste(setdiff(bases, ref[multi])[1:2], collapse = ",")
  calls[1, multi] <- "1/2:15"
  hdr <- c("##fileformat=VCFv4.2",
           "##source=quardplot::simulate_example_data (simulated data)",
           sprintf("##contig=<ID=Chr%s,length=%.0f>", sizes$chr, sizes$length),
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
                 collapse = "\t"))
  body <- vapply(seq_len(n_snp), function(k) {
    paste(c(paste0("Chr", chr_t), pos[k], ids_reg[k], ref[k], alt_vcf[k], ".", "PASS", ".", "GT:DP",
            calls[, k]), collapse = "\t")
  }, "")
  f_vcf <- file.path(dir, "example_region.vcf.gz")
  con <- gzfile(f_vcf, "w")
  writeLines(c(hdr, body), con)
  close(con)
  hm <- vapply(seq_len(n_snp), function(k) {
    a <- ref[k]
    b <- alt[k]
    g <- geno[, k]
    cl <- ifelse(is.na(g), "NN", ifelse(g == 0L, paste0(a, a), ifelse(g == 1L, paste0(a, b), paste0(b, b))))
    paste(c(ids_reg[k], paste0(a, "/", b), chr_t, pos[k], "+", "NA", "NA", "NA", "NA", "NA", "NA", cl),
          collapse = "\t")
  }, "")
  f_hmp <- file.path(dir, "example_region.hmp.txt")
  writeLines(c(paste(c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID",
                       "assayLSID", "panelLSID", "QCcode", samples), collapse = "\t"), hm), f_hmp)

  # gene models
  f_gff <- file.path(dir, "example_genes.gff3")
  writeLines(.sim_gff_lines("Chr9", r0, r1), f_gff)

  # LD blocks from the genotypes (PLINK .blocks.det format)
  g <- genotype_qc(read_vcf(f_vcf, chr = chr_t, verbose = FALSE), verbose = FALSE)
  bl <- find_ld_blocks(g, verbose = FALSE)
  f_blk <- file.path(dir, "example_blocks.det")
  write_blocks(bl, f_blk)

  # causal variant: common variant inside a multi-SNP block near the centre
  af <- colMeans(geno, na.rm = TRUE) / 2
  ok_af <- af > 0.25 & af < 0.75
  ok_af[c(rare, multi)] <- FALSE
  big <- bl[bl$nsnps >= 5, , drop = FALSE]
  causal <- NA_integer_
  if (nrow(big)) {
    big <- big[order(abs((big$start + big$end) / 2 - 12.55e6)), , drop = FALSE]
    for (b in seq_len(nrow(big))) {
      cand <- which(pos >= big$start[b] & pos <= big$end[b] & ok_af)
      if (length(cand)) {
        r2b <- suppressWarnings(stats::cor(geno[, cand, drop = FALSE], use = "pairwise.complete.obs")^2)
        causal <- cand[which.max(colSums(r2b, na.rm = TRUE))]
        break
      }
    }
  }
  if (is.na(causal)) causal <- which(ok_af)[which.min(abs(pos[ok_af] - 12.55e6))]

  # phenotype and association P values (simple linear regression)
  gc <- geno[, causal]
  gc[is.na(gc)] <- round(mean(gc, na.rm = TRUE))
  y <- 1.3 * gc / 2 + stats::rnorm(n_samples)
  pval <- vapply(seq_len(n_snp), function(k) {
    gk <- geno[, k]
    ok <- !is.na(gk)
    if (stats::var(gk[ok]) == 0) return(NA_real_)
    stats::cor.test(gk[ok], y[ok])$p.value
  }, 0)
  f_phe <- file.path(dir, "example_phenotype.csv")
  utils::write.csv(data.frame(Taxa = samples, Trait = round(y, 4)), f_phe, row.names = FALSE, quote = FALSE)

  # genome-wide GWAS results: background SNPs, two other peaks, and the region
  n_bg <- 15000L
  chr_bg <- sample(sizes$chr, n_bg, replace = TRUE, prob = sizes$length)
  pos_bg <- round(stats::runif(n_bg, 1, sizes$length[match(chr_bg, sizes$chr)]))
  p_bg <- stats::runif(n_bg)
  peaks <- data.frame(chr = c("3", "11"), pos = c(20e6, 8e6), h = c(5.5, 6.5))
  for (k in seq_len(nrow(peaks))) {
    on <- chr_bg == peaks$chr[k] & abs(pos_bg - peaks$pos[k]) < 1.5e6
    lp <- peaks$h[k] * exp(-abs(pos_bg[on] - peaks$pos[k]) / 3e5) * stats::runif(sum(on), 0.5, 1)
    p_bg[on] <- pmin(p_bg[on], 10^(-lp))
  }
  drop <- chr_bg == chr_t & pos_bg >= r0 - 1000 & pos_bg <= r1 + 1000
  okp <- !is.na(pval)
  gwas <- data.frame(SNP = c(sprintf("S%s_%d", chr_bg[!drop], pos_bg[!drop]), ids_reg[okp]),
                     CHR = c(chr_bg[!drop], rep(chr_t, sum(okp))),
                     BP = c(pos_bg[!drop], pos[okp]),
                     P = signif(c(p_bg[!drop], pval[okp]), 4), stringsAsFactors = FALSE)
  gwas <- gwas[order(as.numeric(gwas$CHR), gwas$BP), , drop = FALSE]
  gwas <- gwas[!duplicated(paste(gwas$CHR, gwas$BP)), , drop = FALSE]
  f_gwas <- file.path(dir, "example_gwas.csv")
  utils::write.csv(gwas, f_gwas, row.names = FALSE, quote = FALSE)

  # whole-genome HapMap (SNPs and +/- indels) with four subpopulations
  gw <- .sim_genome_hapmap(sizes, samples)
  f_ghm <- file.path(dir, "example_genome.hmp.txt")
  writeLines(gw$lines, f_ghm)
  f_grp <- file.path(dir, "example_groups.csv")
  utils::write.csv(data.frame(Taxa = samples, Group = gw$groups), f_grp, row.names = FALSE, quote = FALSE)

  # complete HapMap = whole genome + the dense chr 9 region, like a real
  # genotype file that is used for both the locus and genome-wide plots
  gb <- gw$lines[-1L]
  g_chr <- .field(gb, 3L, "\t")
  g_pos <- as.numeric(.field(gb, 4L, "\t"))
  keep <- !(g_chr == chr_t & g_pos >= r0 & g_pos <= r1)
  a_body <- c(gb[keep], hm)
  o <- order(as.numeric(c(g_chr[keep], rep(chr_t, length(hm)))), c(g_pos[keep], pos))
  f_all <- file.path(dir, "example_complete.hmp.txt")
  writeLines(c(gw$lines[1L], a_body[o]), f_all)

  out <- list(gwas = f_gwas, gff = f_gff, vcf = f_vcf, hapmap = f_hmp, phenotype = f_phe,
              blocks = f_blk, genome_hapmap = f_ghm, complete_hapmap = f_all, groups = f_grp,
              chr = chr_t, pos = pos[causal])
  message("Simulated example data written to '", dir, "' (target: chr ", chr_t, ":", .fmt(pos[causal]), ").")
  invisible(out)
}

# Whole-genome HapMap lines: inbred accessions from four subpopulations
# (Balding-Nichols allele frequencies, Fst ~ 0.25), SNP density higher towards
# chromosome ends, plus +/- indels.
.sim_genome_hapmap <- function(sizes, samples, n_snp = 3000L, n_indel = 400L, fst = 0.25) {
  n <- length(samples)
  lv <- c("Group A", "Group B", "Group C", "Group D")
  grp <- sample(rep(lv, length.out = n))
  m <- n_snp + n_indel
  chr <- sample(sizes$chr, m, replace = TRUE, prob = sizes$length)
  pos <- pmax(1, round(sizes$length[match(chr, sizes$chr)] * stats::rbeta(m, 0.75, 0.75)))
  d <- data.frame(chr = chr, pos = pos, indel = c(rep(FALSE, n_snp), rep(TRUE, n_indel)))
  d <- d[!duplicated(paste(d$chr, d$pos)), , drop = FALSE]
  d <- d[order(as.numeric(d$chr), d$pos), , drop = FALSE]
  m <- nrow(d)
  p_anc <- stats::runif(m, 0.05, 0.95)
  a <- p_anc * (1 - fst) / fst
  b <- (1 - p_anc) * (1 - fst) / fst
  P <- sapply(lv, function(z) stats::rbeta(m, a, b))
  gi <- match(grp, lv)
  alt <- matrix(stats::runif(m * n), m) < P[, gi]
  dos <- 2L * alt
  dos[matrix(stats::runif(m * n) < 0.015, m)] <- 1L
  dos[matrix(stats::runif(m * n) < 0.03, m)] <- NA
  bases <- c("A", "C", "G", "T")
  ref <- sample(bases, m, replace = TRUE)
  altb <- vapply(ref, function(z) sample(setdiff(bases, z), 1), "")
  ref[d$indel] <- ifelse(stats::runif(sum(d$indel)) < 0.5, "+", "-")
  altb[d$indel] <- ifelse(ref[d$indel] == "+", "-", "+")
  calls <- ifelse(is.na(dos), "NN", ifelse(dos == 0L, paste0(ref, ref), ifelse(dos == 1L, paste0(ref, altb),
                                                                               paste0(altb, altb))))
  calls <- matrix(calls, m)
  ids <- sprintf("%s_%s_%d", ifelse(d$indel, "I", "S"), d$chr, d$pos)
  meta <- paste(ids, paste0(ref, "/", altb), d$chr, d$pos, "+", "NA", "NA", "NA", "NA", "NA", "NA", sep = "\t")
  body <- do.call(paste, c(list(meta), as.data.frame(calls, stringsAsFactors = FALSE), sep = "\t"))
  hdr <- paste(c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID",
                 "panelLSID", "QCcode", samples), collapse = "\t")
  list(lines = c(hdr, body), groups = grp)
}
