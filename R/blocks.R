# LD blocks --------------------------------------------------------------------

.empty_blocks <- function() {
  data.frame(chr = character(), start = numeric(), end = numeric(), kb = numeric(),
             nsnps = integer(), snps = character(), stringsAsFactors = FALSE)
}

#' Read an LD-block file
#'
#' Accepts the `.blocks.det` file written by PLINK 1.9 `--blocks`
#' (columns `CHR BP1 BP2 KB NSNPS SNPS`; PLINK 2.0 has no `--blocks`), a
#' file written by [write_blocks()], or any table with columns `CHR`, `START`
#' and `END` (the format of the original script). Chromosome names are
#' normalised.
#'
#' @param file Path to the block file, or a data frame.
#' @param chr_map Optional name mapping passed to [normalize_chr()].
#' @return Data frame with `chr`, `start`, `end` and, when available,
#'   `nsnps` and `snps`.
#' @export
read_blocks <- function(file, chr_map = NULL) {
  if (is.data.frame(file)) {
    d <- file
  } else {
    .check_file(file, "LD block")
    first <- trimws(readLines(file, n = 1L, warn = FALSE))
    if (startsWith(first, "*")) {
      .stop("'", file, "' is a PLINK .blocks file (marker names only). Use the .blocks.det file, ",
            "which contains the block coordinates.")
    }
    d <- .fread_df(file, header = TRUE)
  }
  nm <- names(d)
  cc <- .find_col(nm, c("CHR", "chrom", "chromosome", "chr_id"))
  sc <- .find_col(nm, c("BP1", "START", "block_start", "bp_start", "from"))
  ec <- .find_col(nm, c("BP2", "END", "block_end", "bp_end", "to"))
  if (anyNA(c(cc, sc, ec))) {
    .stop("LD-block file needs columns CHR, BP1, BP2 (PLINK 1.9 --blocks .blocks.det) or ",
          "CHR, START, END. Columns found: ", paste(nm, collapse = ", "), ".")
  }
  out <- data.frame(chr = normalize_chr(d[[cc]], chr_map), start = as.numeric(d[[sc]]),
                    end = as.numeric(d[[ec]]), stringsAsFactors = FALSE)
  nc <- .find_col(nm, c("NSNPS", "nsnp", "n_snps"))
  if (!is.na(nc)) out$nsnps <- as.integer(d[[nc]])
  snc <- .find_col(nm, c("SNPS", "snp_list", "markers"))
  if (!is.na(snc)) out$snps <- as.character(d[[snc]])
  out <- out[!is.na(out$chr) & is.finite(out$start) & is.finite(out$end), , drop = FALSE]
  out <- out[order(match(out$chr, sort_chr(out$chr)), out$start), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Write LD blocks in PLINK `.blocks.det` format
#'
#' @param blocks Output of [find_ld_blocks()] or [read_blocks()].
#' @param file Output path.
#' @return The file path (invisibly).
#' @export
write_blocks <- function(blocks, file) {
  out <- data.frame(CHR = blocks$chr, BP1 = blocks$start, BP2 = blocks$end,
                    KB = round((blocks$end - blocks$start + 1) / 1000, 3),
                    NSNPS = if (!is.null(blocks$nsnps)) blocks$nsnps else NA_integer_,
                    SNPS = if (!is.null(blocks$snps)) blocks$snps else NA_character_)
  .ensure_dir(file)
  utils::write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(file)
}

.prefix2d <- function(S) {
  m <- nrow(S)
  C <- matrix(0, m + 1L, m + 1L)
  C[-1, -1] <- t(apply(apply(S, 2, cumsum), 1, cumsum))
  C
}

# 90% confidence bounds of D' per pair from the likelihood surface over
# D' = 0, 0.01, ..., 1 (Gabriel et al. 2002; Haploview / PLINK 1.9 rule).
.dprime_ci <- function(k, h, chunk = 20000L) {
  P <- length(k$k11)
  lo <- hi <- rep(NA_real_, P)
  pA <- h$p11 + h$p10
  pB <- h$p11 + h$p01
  D <- h$p11 - pA * pB
  neg <- !is.na(D) & D < 0
  k11 <- ifelse(neg, k$k10, k$k11)
  k10 <- ifelse(neg, k$k11, k$k10)
  k01 <- ifelse(neg, k$k00, k$k01)
  k00 <- ifelse(neg, k$k01, k$k00)
  pBf <- ifelse(neg, 1 - pB, pB)
  dmax <- pmin(pA * (1 - pBf), (1 - pA) * pBf)
  ok <- is.finite(pA) & is.finite(pBf) & pA > 0 & pA < 1 & pBf > 0 & pBf < 1 & is.finite(dmax) & dmax > 0
  todo <- which(ok)
  grid <- seq(0, 1, by = 0.01)
  for (s in seq(1L, length(todo), by = chunk)) {
    idx <- todo[s:min(length(todo), s + chunk - 1L)]
    a <- pA[idx]
    b <- pBf[idx]
    dm <- dmax[idx]
    LL <- matrix(0, length(idx), length(grid))
    for (gi in seq_along(grid)) {
      Dg <- grid[gi] * dm
      q11 <- pmax(a * b + Dg, 1e-10)
      q10 <- pmax(a * (1 - b) - Dg, 1e-10)
      q01 <- pmax((1 - a) * b - Dg, 1e-10)
      q00 <- pmax((1 - a) * (1 - b) + Dg, 1e-10)
      LL[, gi] <- k11[idx] * log(q11) + k10[idx] * log(q10) + k01[idx] * log(q01) +
        k00[idx] * log(q00) + k$dh[idx] * log(q11 * q00 + q10 * q01)
    }
    mx <- LL[cbind(seq_along(idx), max.col(LL, ties.method = "first"))]
    L <- exp(LL - mx)
    tot <- rowSums(L)
    cum <- L
    for (gi in 2:101) cum[, gi] <- cum[, gi - 1L] + L[, gi]
    rcum <- L
    for (gi in 100:1) rcum[, gi] <- rcum[, gi + 1L] + L[, gi]
    first <- max.col(cum > 0.05 * tot, ties.method = "first")
    last <- 102L - max.col((rcum > 0.05 * tot)[, 101:1, drop = FALSE], ties.method = "first")
    lo[idx] <- pmax(first - 2L, 0L) / 100
    hi[idx] <- pmin(last, 100L) / 100
  }
  list(lo = lo, hi = hi)
}

#' Identify haplotype blocks from genotypes
#'
#' Implements the confidence-interval method of Gabriel et al. (2002), as
#' used by Haploview and by PLINK 1.9 `--blocks`, so an external block file
#' is not needed. For each pair of variants within `max_kb`, haplotype
#' frequencies are estimated by EM and a 90% confidence interval for D' is
#' obtained from the likelihood over D' = 0, 0.01, ..., 1. Pairs are in
#' "strong LD" when the lower bound is at least `strong_lowci` and the upper
#' bound at least `strong_highci`, and show "strong evidence of recombination"
#' when the upper bound is below `recomb_highci`. A block is a run of variants in which
#' strong-LD pairs make up more than `inform_frac` of the informative pairs;
#' blocks are built from the longest strong-LD pairs first. As in Haploview,
#' two- and three-variant blocks may span at most 20 and 30 kb, and small
#' blocks use the lower-bound cut-offs 0.80 (two variants) and 0.50 (three
#' or four variants).
#'
#' Requires diploid genotypes coded 0/1/2 (or phased haplotypes). Apply
#' [genotype_qc()] first (PLINK uses MAF >= 0.05 by default).
#'
#' @param geno A `qp_geno` object (one region; a few thousand variants at most).
#' @param max_kb Maximum distance between variants considered (PLINK default 200).
#' @param strong_lowci,strong_highci,recomb_highci,inform_frac Gabriel
#'   method thresholds (Haploview / PLINK defaults).
#' @param small_max_span Apply the 20/30 kb span limits to 2- and 3-variant blocks.
#' @param verbose Print a summary.
#' @return Data frame of blocks: `chr`, `start`, `end`, `kb`, `nsnps`,
#'   `snps` (IDs separated by `|`), `first` and `last` (row indices in the
#'   genotype object).
#' @references Gabriel SB et al. (2002) The structure of haplotype blocks in
#'   the human genome. Science 296:2225-2229. Barrett JC et al. (2005)
#'   Haploview: analysis and visualization of LD and haplotype maps.
#'   Bioinformatics 21:263-265.
#' @export
find_ld_blocks <- function(geno, max_kb = 200, strong_lowci = 0.70, strong_highci = 0.98,
                           recomb_highci = 0.90, inform_frac = 0.95, small_max_span = TRUE,
                           verbose = TRUE) {
  if (!inherits(geno, "qp_geno")) .stop("`geno` must be a qp_geno object.")
  chrs <- unique(geno$info$chr)
  if (length(chrs) > 1L) {
    res <- lapply(chrs, function(cc) {
      find_ld_blocks(.geno_subset(geno, which(geno$info$chr == cc)), max_kb, strong_lowci,
                     strong_highci, recomb_highci, inform_frac, small_max_span, verbose = FALSE)
    })
    return(do.call(rbind, res))
  }
  if (nrow(geno$info) < 2L) return(.empty_blocks())
  g <- if (is.unsorted(geno$info$pos)) .geno_subset(geno, order(geno$info$pos)) else geno
  pos <- g$info$pos
  ids <- g$info$id
  m <- length(pos)
  if (m > 4000L) .warn("Block detection on ", .fmt(m), " variants needs a lot of memory; use a smaller region.")
  t0 <- proc.time()[["elapsed"]]
  last <- findInterval(pos + max_kb * 1000, pos)
  nj <- pmax(last - seq_len(m), 0L)
  pi <- rep.int(seq_len(m), nj)
  pj <- sequence(nj, from = seq_len(m) + 1L)
  P <- length(pi)
  if (!P) return(.empty_blocks())
  phased <- !is.null(g$hap1)
  if (phased) {
    H <- rbind(t(g$hap1), t(g$hap2))
  } else {
    .check_diploid(g$dosage)
    X <- t(g$dosage)
  }
  k <- list(k11 = numeric(P), k10 = numeric(P), k01 = numeric(P), k00 = numeric(P), dh = numeric(P))
  B <- 250L
  for (b0 in seq(1L, m, by = B)) {
    b1 <- min(m, b0 + B - 1L)
    sel <- which(pi >= b0 & pi <= b1)
    if (!length(sel)) next
    c0 <- min(pj[sel])
    c1 <- max(pj[sel])
    kk <- if (phased) {
      .hap_counts_from_phased(H[, b0:b1, drop = FALSE], H[, c0:c1, drop = FALSE])
    } else {
      .hap_counts_from_geno(X[, b0:b1, drop = FALSE], X[, c0:c1, drop = FALSE])
    }
    rc <- cbind(pi[sel] - b0 + 1L, pj[sel] - c0 + 1L)
    for (nm in names(k)) k[[nm]][sel] <- kk[[nm]][rc]
  }
  h <- .em_haplo(k)
  ci <- .dprime_ci(k, h)
  lo <- ci$lo
  hi <- ci$hi
  cand <- which(!is.na(lo) & lo >= strong_lowci & hi >= strong_highci)
  if (!length(cand)) {
    .msg(verbose, "LD blocks: none found among ", .fmt(m), " variants.")
    return(.empty_blocks())
  }
  strong_in <- !is.na(lo) & lo > strong_lowci & hi >= strong_highci
  rec <- !is.na(hi) & hi < recomb_highci
  S <- matrix(0, m, m)
  R <- matrix(0, m, m)
  S[cbind(pi[strong_in], pj[strong_in])] <- 1
  R[cbind(pi[rec], pj[rec])] <- 1
  CS <- .prefix2d(S)
  CR <- .prefix2d(R)
  rm(S, R)
  box <- function(C, i, j) C[j + 1L, j + 1L] - C[i, j + 1L] - C[j + 1L, i] + C[i, i]
  keys <- pi * (m + 1) + pj
  sep <- pos[pj[cand]] - pos[pi[cand]]
  cand <- cand[order(-sep, -cand)]
  # Haploview arrays indexed by the number of markers (0-based in Java)
  cut_low_small <- c(0, 0, 0.80, 0.50, 0.50)
  max_span_small <- c(0, 0, 20000, 30000)
  used <- logical(m)
  bi <- integer()
  bj <- integer()
  for (cc in cand) {
    i <- pi[cc]
    j <- pj[cc]
    if (used[i] || used[j]) next
    n_in <- j - i + 1L
    if (small_max_span && n_in < 4L && (pos[j] - pos[i]) > max_span_small[n_in + 1L]) next
    if (n_in < 5L) {
      xy <- utils::combn(i:j, 2)
      kk <- match(xy[1, ] * (m + 1) + xy[2, ], keys)
      l <- lo[kk]
      u <- hi[kk]
      okp <- !is.na(kk) & !is.na(l)
      n_strong <- sum(okp & l > cut_low_small[n_in + 1L] & u >= strong_highci)
      n_rec <- sum(okp & u < recomb_highci)
    } else {
      n_strong <- box(CS, i, j)
      n_rec <- box(CR, i, j)
    }
    need <- if (n_in > 3L) 6 else if (n_in > 2L) 3 else 1
    if (n_strong + n_rec < need) next
    if (n_strong / (n_strong + n_rec) > inform_frac) {
      bi <- c(bi, i)
      bj <- c(bj, j)
      used[i:j] <- TRUE
    }
  }
  if (!length(bi)) {
    .msg(verbose, "LD blocks: none found among ", .fmt(m), " variants.")
    return(.empty_blocks())
  }
  o <- order(bi)
  bi <- bi[o]
  bj <- bj[o]
  out <- data.frame(chr = g$info$chr[1], start = pos[bi], end = pos[bj],
                    kb = (pos[bj] - pos[bi] + 1) / 1000, nsnps = bj - bi + 1L,
                    snps = vapply(seq_along(bi), function(z) paste(ids[bi[z]:bj[z]], collapse = "|"), ""),
                    first = bi, last = bj, stringsAsFactors = FALSE)
  .msg(verbose, "LD blocks: ", nrow(out), " block(s) among ", .fmt(m), " variants (Gabriel method, ",
       .fmt(P), " pairs within ", max_kb, " kb) in ", round(proc.time()[["elapsed"]] - t0, 1), " s.")
  out
}
