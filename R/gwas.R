# GWAS summary statistics --------------------------------------------------

.gwas_cols <- list(
  snp = c("SNP", "rs", "rsid", "rs_id", "marker", "markername", "marker_id", "snp_id",
          "variant_id", "variant", "id", "taxa"),
  chr = c("CHR", "chrom", "chromosome", "chr_id", "chr_name", "#chrom", "seqnames", "seqid"),
  pos = c("BP", "pos", "position", "ps", "base_pair_location", "bp_position",
          "physical_position", "coordinate"),
  p   = c("P", "p_value", "pvalue", "pval", "p.value", "p_wald", "p_lrt", "p_score",
          "p_bolt_lmm", "p_bolt_lmm_inf", "p_val")
)

#' Read GWAS summary statistics
#'
#' Reads comma-, tab- or space-separated GWAS results (optionally gzip, bzip2
#' or xz compressed). Columns are found by common names (PLINK `SNP CHR BP P`,
#' GAPIT `SNP Chromosome Position P.value`, GEMMA `rs chr ps p_wald`, TASSEL
#' `Marker Chr Pos p`, ...) or set explicitly. If no column name is
#' recognised, the first four columns are taken as SNP, chromosome, position
#' and P value (the behaviour of the original script), with a message.
#'
#' rMVP results, whose P value column is named after the trait and model
#' (e.g. `Height.MLM`), are recognised, and for PLINK `--linear`,
#' `--logistic` or `--glm` output with a `TEST` column only the `ADD` rows
#' (the SNP effect) are kept.
#'
#' Quality checks: rows with missing chromosome, position or P value are
#' removed; P values must lie in (0, 1] (values of exactly 0 are replaced by
#' the smallest non-zero P value with a warning; values outside \[0, 1\] are
#' removed with a warning); positions must be positive integers (1-based);
#' duplicated chromosome-position pairs keep the smallest P value; missing SNP
#' identifiers are replaced by `chr:pos`.
#'
#' @param file Path to the GWAS file, or a data frame.
#' @param snp_col,chr_col,pos_col,p_col Optional column names (or numbers) to
#'   use instead of automatic detection.
#' @param sep Field separator passed to [data.table::fread()]; `"auto"` detects it.
#' @param chr_map Optional name mapping passed to [normalize_chr()].
#' @param verbose Print progress messages.
#' @return A data frame of class `qp_gwas` with columns `snp`, `chr`
#'   (normalised), `pos`, `p` and `logp` (-log10 P), sorted by chromosome and
#'   position.
#' @export
read_gwas <- function(file, snp_col = NULL, chr_col = NULL, pos_col = NULL, p_col = NULL,
                      sep = "auto", chr_map = NULL, verbose = TRUE) {
  if (is.data.frame(file)) {
    raw <- as.data.frame(file, stringsAsFactors = FALSE)
    src <- "data frame"
  } else {
    .check_file(file, "GWAS")
    src <- basename(file)
    first <- readLines(file, n = 1L, warn = FALSE)
    if (grepl("^[[:space:]]*\".*[,;\t].*\"[[:space:]]*$", first)) {
      # Whole rows stored as quoted strings ("SNP1,Chr1,1250432,1e-5"), as
      # handled by the original script: remove the outer quotes and split.
      txt <- readLines(file, warn = FALSE)
      txt <- sub("\"[[:space:]]*$", "", sub("^[[:space:]]*\"", "", txt))
      raw <- as.data.frame(data.table::fread(text = txt[nzchar(txt)], header = TRUE,
                                             colClasses = "character", data.table = FALSE,
                                             showProgress = FALSE),
                           stringsAsFactors = FALSE)
    } else {
      raw <- .fread_df(file, sep = sep, header = TRUE, colClasses = "character")
    }
  }
  if (!nrow(raw)) .stop("GWAS input '", src, "' has no rows.")
  nm <- names(raw)
  pick <- function(user, key) {
    if (!is.null(user)) {
      if (is.numeric(user)) return(nm[user])
      if (!user %in% nm) .stop("Column '", user, "' not found in GWAS file. Columns: ",
                               paste(nm, collapse = ", "))
      return(user)
    }
    .find_col(nm, .gwas_cols[[key]])
  }
  cols <- list(snp = pick(snp_col, "snp"), chr = pick(chr_col, "chr"),
               pos = pick(pos_col, "pos"), p = pick(p_col, "p"))
  if (is.na(cols$p) && is.null(p_col)) {
    # rMVP names the P value column after the trait and model, e.g. "Height.MLM"
    rmvp <- grep("\\.(GLM|MLM|FarmCPU|BLINK)$", nm, ignore.case = TRUE, value = TRUE)
    if (length(rmvp)) {
      cols$p <- rmvp[1]
      .msg(verbose, "GWAS: using column '", rmvp[1], "' as the P value (rMVP output)",
           if (length(rmvp) > 1L) paste0("; other P value columns ignored: ", paste(rmvp[-1], collapse = ", ")),
           ".")
    }
  }
  test_col <- .find_col(nm, "TEST")
  if (!is.na(test_col)) {
    # PLINK --linear/--logistic/--glm write one row per term; keep the SNP effect
    tv <- toupper(trimws(raw[[test_col]]))
    if (any(tv == "ADD") && !all(tv == "ADD")) {
      .msg(verbose, "GWAS: kept the ", .fmt(sum(tv == "ADD")), " row(s) with TEST = ADD ",
           "(rows for covariates and other terms removed).")
      raw <- raw[tv == "ADD", , drop = FALSE]
    }
  }
  found <- !vapply(cols[c("chr", "pos", "p")], is.na, NA)
  if (!all(found)) {
    if (!any(found) && is.na(cols$snp) && ncol(raw) >= 4L) {
      cols <- list(snp = nm[1], chr = nm[2], pos = nm[3], p = nm[4])
      .msg(verbose, "GWAS column names not recognised; using columns 1-4 as SNP, ",
           "chromosome, position and P value (", paste(nm[1:4], collapse = ", "), ").")
    } else {
      .stop("Could not identify the ", paste(c("chromosome", "position", "P value")[!found],
                                              collapse = ", "),
            " column(s) in GWAS input '", src, "'. Columns found: ", paste(nm, collapse = ", "),
            ". Set them with chr_col, pos_col and p_col.")
    }
  }
  n0 <- nrow(raw)
  chr <- normalize_chr(raw[[cols$chr]], chr_map)
  pos <- suppressWarnings(as.numeric(raw[[cols$pos]]))
  p <- suppressWarnings(as.numeric(raw[[cols$p]]))
  snp <- if (!is.na(cols$snp)) as.character(raw[[cols$snp]]) else rep(NA_character_, n0)

  bad <- is.na(chr) | !nzchar(chr) | is.na(pos) | is.na(p)
  if (any(bad)) .msg(verbose, "GWAS: removed ", .fmt(sum(bad)), " row(s) with missing chromosome, position or P value.")
  badpos <- !bad & (pos < 1 | pos != round(pos))
  if (any(badpos)) .warn("GWAS: removed ", .fmt(sum(badpos)), " row(s) with non-positive or non-integer positions.")
  badp <- !bad & (p < 0 | p > 1)
  if (any(badp)) .warn("GWAS: removed ", .fmt(sum(badp)), " row(s) with P values outside [0, 1]",
                       " (are these already -log10 values?).")
  keep <- !(bad | badpos | badp)
  if (!any(keep)) .stop("No valid GWAS rows remain in '", src, "' after quality checks.")
  d <- data.frame(snp = snp[keep], chr = chr[keep], pos = pos[keep], p = p[keep],
                  stringsAsFactors = FALSE)
  zero <- d$p == 0
  if (any(zero)) {
    pmin_nz <- if (any(!zero)) min(d$p[!zero]) else .Machine$double.xmin
    d$p[zero] <- pmin_nz
    .warn("GWAS: ", .fmt(sum(zero)), " P value(s) of exactly 0 were set to the smallest non-zero P value (",
          format(pmin_nz, digits = 3), ").")
  }
  miss_id <- is.na(d$snp) | !nzchar(d$snp) | d$snp == "."
  d$snp[miss_id] <- paste0(d$chr[miss_id], ":", d$pos[miss_id])
  d <- d[order(match(d$chr, sort_chr(d$chr)), d$pos, d$p), , drop = FALSE]
  dup <- duplicated(paste(d$chr, d$pos))
  if (any(dup)) {
    .msg(verbose, "GWAS: ", .fmt(sum(dup)), " duplicated chromosome-position row(s); kept the smallest P value.")
    d <- d[!dup, , drop = FALSE]
  }
  dup_id <- duplicated(d$snp)
  if (any(dup_id)) {
    .warn("GWAS: ", .fmt(sum(dup_id)), " SNP identifier(s) occur at more than one position; ",
          "they were made unique by appending chr:pos.")
    d$snp[dup_id] <- paste0(d$snp[dup_id], "_", d$chr[dup_id], ":", d$pos[dup_id])
  }
  d$logp <- -log10(d$p)
  rownames(d) <- NULL
  .msg(verbose, "GWAS: ", .fmt(nrow(d)), " SNPs on ", length(unique(d$chr)),
       " chromosome(s) read from ", src, ".")
  attr(d, "n_input") <- n0
  attr(d, "columns") <- unlist(cols)
  class(d) <- c("qp_gwas", "data.frame")
  d
}

.as_gwas <- function(gwas, verbose = FALSE) {
  if (inherits(gwas, "qp_gwas")) return(gwas)
  read_gwas(gwas, verbose = verbose)
}

#' Genome-wide significance threshold
#'
#' @param gwas GWAS data ([read_gwas()]); used for the number of tests.
#' @param threshold `"bonferroni"` (alpha / number of SNPs), a P value
#'   (e.g. `5e-8`), a value on the -log10 scale (any number >= 1, e.g. `7`
#'   for P = 1e-7), or `NULL`/`"none"` for no line.
#' @param alpha Family-wise error rate for the Bonferroni threshold.
#' @return `NULL` or a list with `p`, `logp` and `label`.
#' @export
gwas_threshold <- function(gwas, threshold = "bonferroni", alpha = 0.05) {
  if (is.list(threshold) && !is.null(threshold$p)) return(threshold)
  if (is.null(threshold) || length(threshold) != 1L || is.na(threshold) ||
      isFALSE(threshold) || identical(tolower(as.character(threshold)), "none")) {
    return(NULL)
  }
  if (is.character(threshold)) {
    if (tolower(threshold) != "bonferroni") {
      num <- suppressWarnings(as.numeric(threshold))
      if (is.na(num)) .stop("`threshold` must be 'bonferroni', 'none' or a number.")
      return(gwas_threshold(gwas, num, alpha))
    }
    n <- nrow(gwas)
    p <- alpha / n
    return(list(p = p, logp = -log10(p),
                label = sprintf("Bonferroni (%s / %s tests)", format(alpha), .fmt(n))))
  }
  if (!is.numeric(threshold) || threshold <= 0) .stop("`threshold` must be positive.")
  p <- if (threshold < 1) threshold else 10^(-threshold)
  list(p = p, logp = -log10(p), label = sprintf("P = %s", format(p, digits = 3)))
}
