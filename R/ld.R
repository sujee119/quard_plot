# Linkage disequilibrium -----------------------------------------------------

# Squared Pearson correlation of allele dosages for all pairs, using for
# each pair only the samples called at both variants (pairwise-complete).
# X: samples x variants. Exactly equal to cor(X, use = "pairwise")^2.
.ld_r2_genotype <- function(X) {
  M <- !is.na(X)
  storage.mode(M) <- "double"
  X0 <- X
  X0[is.na(X0)] <- 0
  N <- crossprod(M)
  Sx <- crossprod(X0, M)
  Sxx <- crossprod(X0^2, M)
  Sxy <- crossprod(X0)
  Sy <- t(Sx)
  Syy <- t(Sxx)
  cv <- Sxy - Sx * Sy / N
  vx <- Sxx - Sx^2 / N
  vy <- Syy - Sy^2 / N
  r2 <- cv^2 / (vx * vy)
  tiny <- 1e-9 * pmax(N, 1)
  r2[!is.finite(r2) | vx <= tiny | vy <= tiny | N < 3] <- NA_real_
  r2[!is.na(r2) & r2 > 1] <- 1
  diag(r2) <- 1
  r2
}

# Unambiguous haplotype counts and double heterozygotes for pairs (i, j)
# from diploid genotypes coded 0/1/2 (ALT copies). Allele "1" = ALT.
.hap_counts_from_geno <- function(Xi, Xj) {
  cls <- function(X, a) {
    m <- (X == a)
    m[is.na(m)] <- FALSE
    storage.mode(m) <- "double"
    m
  }
  I <- list(cls(Xi, 0), cls(Xi, 1), cls(Xi, 2))
  J <- list(cls(Xj, 0), cls(Xj, 1), cls(Xj, 2))
  n <- function(a, b) crossprod(I[[a + 1]], J[[b + 1]])
  n00 <- n(0, 0); n01 <- n(0, 1); n02 <- n(0, 2)
  n10 <- n(1, 0); n11 <- n(1, 1); n12 <- n(1, 2)
  n20 <- n(2, 0); n21 <- n(2, 1); n22 <- n(2, 2)
  list(k11 = 2 * n22 + n21 + n12, k10 = 2 * n20 + n21 + n10,
       k01 = 2 * n02 + n12 + n01, k00 = 2 * n00 + n01 + n10, dh = n11)
}

# Exact haplotype counts from phased haplotypes (stacked haplotype matrices).
.hap_counts_from_phased <- function(Hi, Hj) {
  ind <- function(H, a) {
    m <- (H == a)
    m[is.na(m)] <- FALSE
    storage.mode(m) <- "double"
    m
  }
  A1 <- ind(Hi, 1); A0 <- ind(Hi, 0)
  B1 <- ind(Hj, 1); B0 <- ind(Hj, 0)
  k11 <- crossprod(A1, B1)
  list(k11 = k11, k10 = crossprod(A1, B0), k01 = crossprod(A0, B1), k00 = crossprod(A0, B0),
       dh = k11 * 0)
}

# EM estimate of two-locus haplotype frequencies (Excoffier & Slatkin 1995),
# vectorised over pairs. Double heterozygotes are split between the
# coupling (11/00) and repulsion (10/01) phases.
.em_haplo <- function(k, tol = 1e-10, max_iter = 1000L) {
  dims <- dim(k$k11)
  K11 <- as.vector(k$k11)
  K10 <- as.vector(k$k10)
  K01 <- as.vector(k$k01)
  K00 <- as.vector(k$k00)
  DH <- as.vector(k$dh)
  N <- K11 + K10 + K01 + K00 + 2 * DH
  N[N == 0] <- NA_real_
  p11 <- (K11 + DH / 2) / N
  p10 <- (K10 + DH / 2) / N
  p01 <- (K01 + DH / 2) / N
  p00 <- (K00 + DH / 2) / N
  # iterate only the pairs with double heterozygotes that have not converged
  act <- which(DH > 0 & !is.na(N))
  it <- 0L
  while (length(act) && it < max_iter) {
    it <- it + 1L
    a <- p11[act] * p00[act]
    b <- p10[act] * p01[act]
    th <- a / (a + b)
    th[!is.finite(th)] <- 0.5
    n11 <- (K11[act] + DH[act] * th) / N[act]
    n10 <- (K10[act] + DH[act] * (1 - th)) / N[act]
    n01 <- (K01[act] + DH[act] * (1 - th)) / N[act]
    n00 <- (K00[act] + DH[act] * th) / N[act]
    d <- pmax(abs(n11 - p11[act]), abs(n10 - p10[act]))
    p11[act] <- n11
    p10[act] <- n10
    p01[act] <- n01
    p00[act] <- n00
    act <- act[!(is.na(d) | d < tol)]
  }
  shape <- function(v) if (is.null(dims)) v else matrix(v, dims[1], dims[2])
  list(p11 = shape(p11), p10 = shape(p10), p01 = shape(p01), p00 = shape(p00), N = shape(N))
}

# r2 and |D'| from haplotype frequencies.
.ld_from_freqs <- function(h) {
  pA <- h$p11 + h$p10
  pB <- h$p11 + h$p01
  D <- h$p11 - pA * pB
  den <- pA * (1 - pA) * pB * (1 - pB)
  r2 <- D^2 / den
  dmax <- ifelse(D >= 0, pmin(pA * (1 - pB), (1 - pA) * pB), pmin(pA * pB, (1 - pA) * (1 - pB)))
  dp <- abs(D) / dmax
  bad <- !is.finite(r2) | den <= 1e-12
  r2[bad] <- NA_real_
  dp[bad | !is.finite(dp)] <- NA_real_
  r2[!is.na(r2) & r2 > 1] <- 1
  dp[!is.na(dp) & dp > 1] <- 1
  list(r2 = r2, dprime = dp, D = D, pA = pA, pB = pB)
}

.check_diploid <- function(D) {
  v <- D[!is.na(D)]
  if (length(v) && any(abs(v - round(v)) > 1e-8)) {
    .stop("Haplotype-based LD (and LD blocks) need diploid genotypes coded 0/1/2. ",
          "Use method = 'genotype' for haploid or polyploid data.")
  }
  invisible(TRUE)
}

#' Pairwise linkage disequilibrium
#'
#' Computes the LD matrix for the variants of a genotype object.
#'
#' * `method = "genotype"` (default): \eqn{r^2} is the squared Pearson
#'   correlation between allele dosages (0, 1, 2 copies of ALT), using for
#'   each pair the samples genotyped at both variants (missing calls are
#'   excluded, never imputed). For inbred or fully homozygous samples, such
#'   as rice germplasm, this equals the haplotype-based \eqn{r^2}; for
#'   outbred diploids it is the "genotypic" (composite) \eqn{r^2}.
#' * `method = "haplotype"`: two-locus haplotype frequencies are estimated by
#'   the EM algorithm from unphased diploid genotypes (counted directly for
#'   phased data), and \eqn{r^2 = D^2 / (p_A(1-p_A)p_B(1-p_B))} and
#'   \eqn{D' = D / D_{max}} are derived from them. `stat = "dprime"` always
#'   uses this method.
#'
#' Values are kept within \[0, 1\]; pairs with no variation among the shared
#' samples are `NA`. Remove rare variants first ([genotype_qc()]) because
#' \eqn{r^2} is unstable when allele frequencies are close to 0.
#'
#' @param geno A `qp_geno` object (ideally after [genotype_qc()]).
#' @param stat `"r2"` or `"dprime"`.
#' @param method `"genotype"` or `"haplotype"` (see Details).
#' @return A symmetric matrix with variant IDs as dimnames and attributes
#'   `pos`, `chr`, `stat` and `method`.
#' @export
calc_ld <- function(geno, stat = c("r2", "dprime"), method = c("genotype", "haplotype")) {
  stat <- match.arg(stat)
  method <- match.arg(method)
  if (!inherits(geno, "qp_geno")) .stop("`geno` must be a qp_geno object.")
  m <- nrow(geno$info)
  if (m < 2L) .stop("At least two variants are needed to compute LD (found ", m, ").")
  if (stat == "dprime") method <- "haplotype"
  if (method == "genotype") {
    out <- .ld_r2_genotype(t(geno$dosage))
  } else {
    if (!is.null(geno$hap1)) {
      H <- rbind(t(geno$hap1), t(geno$hap2))
      k <- .hap_counts_from_phased(H, H)
    } else {
      .check_diploid(geno$dosage)
      X <- t(geno$dosage)
      k <- .hap_counts_from_geno(X, X)
    }
    ld <- .ld_from_freqs(.em_haplo(k))
    out <- if (stat == "r2") ld$r2 else ld$dprime
    diag(out) <- 1
  }
  dimnames(out) <- list(geno$info$id, geno$info$id)
  attr(out, "pos") <- geno$info$pos
  attr(out, "chr") <- geno$info$chr[1]
  attr(out, "stat") <- stat
  attr(out, "method") <- method
  out
}

#' LD of every variant with a lead variant
#'
#' @param geno A `qp_geno` object.
#' @param lead Lead variant: an ID in `geno$info$id`, a row number, or a
#'   position on the chromosome.
#' @return Named numeric vector of genotypic \eqn{r^2} with the lead variant
#'   (names are variant IDs; the lead has 1).
#' @export
ld_with_lead <- function(geno, lead) {
  info <- geno$info
  i <- if (is.character(lead)) {
    match(lead, info$id)
  } else if (lead %in% info$pos) {
    match(lead, info$pos)
  } else if (lead >= 1 && lead <= nrow(info) && lead == round(lead)) {
    as.integer(lead)
  } else {
    NA_integer_
  }
  if (is.na(i)) .stop("Lead variant '", lead, "' is not in the genotype data.")
  X <- t(geno$dosage)
  x <- X[, i]
  ml <- !is.na(x)
  x0 <- ifelse(ml, x, 0)
  M <- !is.na(X)
  X0 <- X
  X0[!M] <- 0
  N <- colSums(M * ml)
  sx <- colSums(M * x0)
  sxx <- colSums(M * x0^2)
  sy <- colSums(X0 * ml)
  syy <- colSums(X0^2 * ml)
  sxy <- colSums(X0 * x0)
  cv <- sxy - sx * sy / N
  vx <- sxx - sx^2 / N
  vy <- syy - sy^2 / N
  r2 <- cv^2 / (vx * vy)
  r2[!is.finite(r2) | vx <= 1e-9 * N | vy <= 1e-9 * N | N < 3] <- NA_real_
  r2[!is.na(r2) & r2 > 1] <- 1
  r2[i] <- 1
  stats::setNames(r2, info$id)
}

#' Read LD values computed by PLINK
#'
#' Loads `plink --r2 square` matrices (with variant IDs from `ids` or a
#' `.bim` file) or the default `plink --r2` table (`CHR_A BP_A SNP_A CHR_B
#' BP_B SNP_B R2`, optionally with a `DP` column) into a matrix that can be
#' passed to [plot_ld()] or compared with [compare_ld()].
#'
#' @param file PLINK `.ld` output file.
#' @param ids Variant IDs in matrix order (square format).
#' @param bim Alternatively, the `.bim` file of the data set used (square
#'   format); its second column gives the IDs.
#' @param stat For table input: `"R2"`, `"R"` (squared) or `"DP"` (D').
#' @return A symmetric matrix with IDs as dimnames (and `pos` attribute
#'   when positions are available).
#' @export
read_plink_ld <- function(file, ids = NULL, bim = NULL, stat = c("R2", "R", "DP")) {
  stat <- match.arg(stat)
  .check_file(file, "PLINK LD")
  first <- readLines(file, n = 1L, warn = FALSE)
  if (grepl("SNP_A", first, fixed = TRUE)) {
    d <- .fread_df(file, header = TRUE)
    col <- if (stat == "DP") "DP" else if ("R2" %in% names(d)) "R2" else "R"
    if (!col %in% names(d)) .stop("Column '", col, "' not found in '", file, "'.")
    v <- as.numeric(d[[col]])
    if (col == "R") v <- v^2
    snps <- unique(c(d$SNP_A, d$SNP_B))
    bp <- c(stats::setNames(d$BP_A, d$SNP_A), stats::setNames(d$BP_B, d$SNP_B))
    pos <- as.numeric(bp[snps])
    o <- order(pos)
    snps <- snps[o]
    pos <- pos[o]
    m <- matrix(NA_real_, length(snps), length(snps), dimnames = list(snps, snps))
    ia <- match(d$SNP_A, snps)
    ib <- match(d$SNP_B, snps)
    m[cbind(ia, ib)] <- v
    m[cbind(ib, ia)] <- v
    diag(m) <- 1
    attr(m, "pos") <- pos
    return(m)
  }
  m <- as.matrix(.fread_df(file, header = FALSE, na.strings = c("nan", "NaN", "NA", "-nan")))
  storage.mode(m) <- "double"
  if (nrow(m) != ncol(m)) .stop("'", file, "' is not a square matrix (use plink --r2 square).")
  if (is.null(ids) && !is.null(bim)) {
    b <- .fread_df(bim, header = FALSE)
    ids <- b[[2]]
    attr(m, "pos") <- as.numeric(b[[4]])
  }
  if (!is.null(ids)) {
    if (length(ids) != nrow(m)) .stop("`ids` has ", length(ids), " entries but the matrix has ", nrow(m), " rows.")
    dimnames(m) <- list(ids, ids)
  }
  m
}

#' Compare two LD matrices
#'
#' Numerical agreement between two LD matrices on their shared variants, for
#' example `calc_ld()` against `plink --r2 square` on the same samples and
#' variants.
#'
#' @param a,b LD matrices with variant IDs as dimnames.
#' @return A one-row data frame: number of shared variants and pairs, maximum
#'   and mean absolute difference, root-mean-square difference and Pearson
#'   correlation of the off-diagonal values.
#' @export
compare_ld <- function(a, b) {
  ids <- intersect(rownames(a), rownames(b))
  if (length(ids) < 2L) .stop("The two matrices share fewer than two variant IDs.")
  a <- a[ids, ids]
  b <- b[ids, ids]
  ut <- upper.tri(a)
  x <- a[ut]
  y <- b[ut]
  ok <- !is.na(x) & !is.na(y)
  d <- abs(x[ok] - y[ok])
  data.frame(n_variants = length(ids), n_pairs = sum(ok), max_abs_diff = max(d),
             mean_abs_diff = mean(d), rmse = sqrt(mean(d^2)),
             correlation = if (sum(ok) > 2) stats::cor(x[ok], y[ok]) else NA_real_)
}

#' Write an LD matrix to a CSV or Excel file
#'
#' Exports LD values (e.g. from [calc_ld()], [read_plink_ld()] or the
#' `ld_all` element returned by [quard_plot()]) for spreadsheets or other
#' software. The file type follows the file name: `.csv` (opens in Excel),
#' `.xlsx` (Excel workbook; needs the `writexl` package) or `.txt`
#' (tab-separated). With `genes`, every variant is annotated with its
#' closest gene and its location (exon / intron / upstream promoter /
#' downstream / intergenic; see [annotate_variants()]).
#'
#' * `format = "matrix"`: one row per variant: `SNP`, `CHR`, `POS`, (gene
#'   annotation columns), then one column per variant with the LD values;
#'   the diagonal is 1.
#' * `format = "pairs"`: one row per pair of variants (`SNP_A`, `POS_A`,
#'   `SNP_B`, `POS_B`, `DIST_BP`, the LD value and, with `genes`, the gene and
#'   location of both variants), optionally only pairs with a value of at
#'   least `min_value`.
#'
#' @param ld LD matrix with variant IDs as dimnames (and `pos` attribute).
#' @param file Output file (`.csv`, `.xlsx` or `.txt`).
#' @param format `"matrix"` or `"pairs"`.
#' @param digits Decimal places.
#' @param min_value For `"pairs"`: keep only pairs with LD >= this value.
#' @param genes Optional GFF3/GTF file or [read_gff()] output for the
#'   gene annotation columns.
#' @param upstream_bp,downstream_bp Promoter (upstream) and downstream
#'   windows for the annotation (default 3 kb and 1 kb).
#' @return The file path (invisibly).
#' @export
write_ld <- function(ld, file, format = c("matrix", "pairs"), digits = 4, min_value = 0,
                     genes = NULL, upstream_bp = 3000, downstream_bp = 1000) {
  format <- match.arg(format)
  out <- .ld_table(ld, format = format, digits = digits, min_value = min_value, genes = genes,
                   upstream_bp = upstream_bp, downstream_bp = downstream_bp)
  write_table(out, file, sheet = if (format == "matrix") "LD" else "LD_pairs")
  invisible(file)
}

# LD matrix as a data frame ("matrix" or "pairs" layout), optionally with the
# closest gene and location of every variant.
.ld_table <- function(ld, format = c("matrix", "pairs"), digits = 4, min_value = 0,
                      genes = NULL, upstream_bp = 3000, downstream_bp = 1000) {
  format <- match.arg(format)
  ids <- rownames(ld) %||% paste0("v", seq_len(nrow(ld)))
  pos <- attr(ld, "pos") %||% rep(NA_real_, nrow(ld))
  chr <- attr(ld, "chr") %||% NA_character_
  tag <- if (identical(attr(ld, "stat"), "dprime")) "Dprime" else "R2"
  ann <- NULL
  if (!is.null(genes)) {
    ann <- annotate_variants(data.frame(snp = ids, chr = chr, pos = pos), genes,
                             upstream_bp = upstream_bp, downstream_bp = downstream_bp)
  }
  if (format == "matrix") {
    m <- round(unclass(as.matrix(ld)), digits)
    attributes(m) <- list(dim = dim(m))
    out <- data.frame(SNP = ids, CHR = chr, POS = pos, stringsAsFactors = FALSE)
    if (!is.null(ann)) {
      out$GENE <- ann$gene_id
      out$LOCATION <- ann$location
      out$FEATURE <- ann$feature
      out$DISTANCE_BP <- ann$distance_bp
      out$GENE_STRAND <- ann$strand
      out$GENE_DESCRIPTION <- ann$description
    }
    lead_cols <- names(out)
    out <- cbind(out, as.data.frame(m, stringsAsFactors = FALSE))
    names(out) <- c(lead_cols, ids)
  } else {
    ut <- which(upper.tri(ld), arr.ind = TRUE)
    v <- round(as.matrix(ld)[ut], digits)
    a <- ut[, 1]
    b <- ut[, 2]
    out <- data.frame(SNP_A = ids[a], POS_A = pos[a], SNP_B = ids[b], POS_B = pos[b],
                      DIST_BP = abs(pos[b] - pos[a]), value = v, stringsAsFactors = FALSE)
    names(out)[6] <- tag
    if (!is.null(ann)) {
      out$GENE_A <- ann$gene_id[a]
      out$LOCATION_A <- ann$location[a]
      out$GENE_B <- ann$gene_id[b]
      out$LOCATION_B <- ann$location[b]
    }
    out <- out[is.na(out[[6]]) | out[[6]] >= min_value, , drop = FALSE]
    out <- out[order(out$POS_A, out$POS_B), , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}
