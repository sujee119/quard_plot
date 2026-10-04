# Genotypes: VCF, TASSEL HapMap, quality control ---------------------------

# k-th tab-separated field of each line (fast for k = 1).
.field <- function(x, k, sep = "\t") {
  if (k == 1L) {
    t1 <- regexpr(sep, x, fixed = TRUE)
    out <- ifelse(t1 > 0, substr(x, 1L, t1 - 1L), x)
  } else {
    s <- if (sep == "\t") "\\t" else sep
    m <- regexpr(sprintf("^(?:[^%s]*%s){%d}[^%s]*", s, s, k - 1L, s), x, perl = TRUE)
    out <- rep(NA_character_, length(x))
    if (any(m > 0)) out[m > 0] <- sub(sprintf("^.*%s", s), "", regmatches(x, m), perl = TRUE)
  }
  if (sep != "\t") out <- gsub("\"", "", out, fixed = TRUE)
  out
}

# Layout of a HapMap file: delimiter (tab, comma or semicolon), the column of
# "rs#" (a leading row-name column, as written by write.csv, is allowed) and
# the first sample column.
.hapmap_layout <- function(file) {
  hdr <- readLines(file, n = 1L, warn = FALSE)
  sep <- if (grepl("\t", hdr, fixed = TRUE)) "\t" else if (grepl(",", hdr, fixed = TRUE)) "," else
    if (grepl(";", hdr, fixed = TRUE)) ";" else " "
  cols <- trimws(strsplit(hdr, sep, fixed = TRUE)[[1]])
  cols <- gsub("^\"|\"$", "", cols)
  r <- which(grepl("^rs#?$", tolower(cols)))[1]
  if (is.na(r) || length(cols) < r + 11L) {
    .stop("'", file, "' does not look like a HapMap file. Expected a header with rs#, alleles, chrom, pos, ",
          "strand, assembly#, center, protLSID, assayLSID, panelLSID, QCcode and then the samples ",
          "(tab- or comma-separated).")
  }
  qc <- match("qccode", tolower(cols))
  fs <- if (!is.na(qc) && qc > r) qc + 1L else r + 11L
  list(sep = sep, header = hdr, cols = cols, r = r, first_sample = fs, samples = cols[fs:length(cols)])
}

.is_hapmap_file <- function(file) {
  first <- readLines(file, n = 1L, warn = FALSE)
  f <- tolower(gsub("\"", "", strsplit(first, "[\t,;]")[[1]]))
  any(f[seq_len(min(2L, length(f)))] %in% c("rs#", "rs"))
}

# Stream a large text file in chunks and keep only the lines on one
# chromosome inside [start, end]. Memory use is bounded by `chunk_size`
# lines; reading stops once the region has been passed in a sorted file.
.stream_lines <- function(file, n_skip, chr_field, pos_field, tchr, start, end, chr_map,
                          chunk_size = 10000L, every = 1L, sep = "\t") {
  con <- file(file, "r")
  on.exit(close(con))
  if (n_skip > 0L) readLines(con, n = n_skip, warn = FALSE)
  out <- list()
  k <- 0L
  cache <- character()
  seen <- FALSE
  sorted <- TRUE
  last_pos <- -Inf
  n_pass <- 0
  thin <- function(v) {
    if (every <= 1L || !length(v)) return(v)
    idx <- n_pass + seq_along(v)
    n_pass <<- n_pass + length(v)
    v[(idx - 1) %% every == 0]
  }
  repeat {
    x <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(x)) break
    x <- x[nzchar(x)]
    if (is.null(tchr)) {
      k <- k + 1L
      out[[k]] <- thin(x)
      next
    }
    lab <- .field(x, chr_field, sep)
    new <- setdiff(unique(lab), names(cache))
    if (length(new)) cache[new] <- normalize_chr(new, chr_map)
    on <- unname(cache[lab]) == tchr
    on[is.na(on)] <- FALSE
    if (any(on)) {
      seen <- TRUE
      xs <- x[on]
      p <- suppressWarnings(as.numeric(.field(xs, pos_field, sep)))
      if (sorted && (is.unsorted(p, na.rm = TRUE) || isTRUE(p[1] < last_pos))) sorted <- FALSE
      last_pos <- p[length(p)]
      sel <- !is.na(p) & p >= start & p <= end
      if (any(sel)) {
        k <- k + 1L
        out[[k]] <- thin(xs[sel])
      }
      if (sorted && isTRUE(last_pos > end)) break
    } else if (seen && sorted) {
      break
    }
  }
  unlist(out, use.names = FALSE)
}

# Extract the GT sub-field from VCF sample columns (FORMAT may be GT:AD:DP:...).
.extract_gt <- function(G, fmt) {
  if (all(fmt == "GT")) return(G)
  out <- G
  gt_first <- grepl("^GT(:|$)", fmt)
  if (any(gt_first)) out[gt_first, ] <- sub(":.*$", "", G[gt_first, , drop = FALSE])
  for (i in which(!gt_first)) {
    f <- strsplit(fmt[i], ":", fixed = TRUE)[[1]]
    k <- match("GT", f)
    if (is.na(k)) .stop("A VCF record has no GT field in FORMAT ('", fmt[i], "').")
    parts <- strsplit(G[i, ], ":", fixed = TRUE)
    out[i, ] <- vapply(parts, function(z) if (length(z) >= k) z[k] else ".", "")
  }
  out
}

# Convert GT strings to allele dosages on the 0-2 scale.
#  * missing calls ('.', './.', '.|.', partial calls) -> NA (never 0)
#  * dosage = 2 * (number of non-reference alleles) / ploidy, so haploid
#    calls ('0', '1') give 0/2 and polyploid calls are scaled to 0-2
#  * for multi-allelic sites every non-reference allele counts as ALT
#  * phased diploid data also return the two haplotypes (0 = REF, 1 = ALT)
.gt_to_dosage <- function(gt) {
  u <- unique(as.vector(gt))
  parts <- strsplit(u, "[/|]")
  n_all <- lengths(parts)
  miss <- vapply(parts, function(a) !length(a) || any(a %in% c(".", "")), NA)
  alt <- vapply(parts, function(a) sum(a != "0"), 0)
  dos <- ifelse(miss, NA_real_, 2 * alt / pmax(n_all, 1L))
  idx <- match(gt, u)
  D <- matrix(dos[idx], nrow = nrow(gt), ncol = ncol(gt))
  phased <- grepl("|", u, fixed = TRUE) & n_all == 2L & !miss
  called <- !miss[idx]
  all_phased <- any(called) && all(phased[idx][called])
  H1 <- H2 <- NULL
  if (all_phased) {
    a1 <- vapply(parts, function(a) if (length(a) >= 1L) as.numeric(a[1] != "0") else NA_real_, 0)
    a2 <- vapply(parts, function(a) if (length(a) >= 2L) as.numeric(a[2] != "0") else NA_real_, 0)
    a1[miss] <- NA_real_
    a2[miss] <- NA_real_
    H1 <- matrix(a1[idx], nrow = nrow(gt))
    H2 <- matrix(a2[idx], nrow = nrow(gt))
  }
  list(dosage = D, hap1 = H1, hap2 = H2)
}

.new_geno <- function(info, dosage, samples, hap1 = NULL, hap2 = NULL, source = NA_character_,
                      how = NA_character_) {
  rownames(info) <- NULL
  dimnames(dosage) <- NULL
  structure(list(info = info, dosage = dosage, samples = samples, hap1 = hap1, hap2 = hap2,
                 source = source, how = how, qc = NULL), class = "qp_geno")
}

.empty_geno <- function(samples = character(), source = NA_character_) {
  info <- data.frame(id = character(), chr = character(), pos = numeric(), ref = character(),
                     alt = character(), multiallelic = logical(), is_snp = logical(),
                     stringsAsFactors = FALSE)
  .new_geno(info, matrix(numeric(), 0, length(samples)), samples, source = source)
}

.geno_subset <- function(g, i) {
  g$info <- g$info[i, , drop = FALSE]
  rownames(g$info) <- NULL
  g$dosage <- g$dosage[i, , drop = FALSE]
  if (!is.null(g$hap1)) {
    g$hap1 <- g$hap1[i, , drop = FALSE]
    g$hap2 <- g$hap2[i, , drop = FALSE]
  }
  g
}

#' @export
print.qp_geno <- function(x, ...) {
  cat(sprintf("<qp_geno> %s variant(s) x %s sample(s)%s\n", .fmt(nrow(x$info)), .fmt(length(x$samples)),
              if (nrow(x$info)) sprintf(" on chr %s:%s-%s", x$info$chr[1], .fmt(min(x$info$pos)),
                                        .fmt(max(x$info$pos))) else ""))
  if (!is.null(x$qc)) print(x$qc, row.names = FALSE)
  invisible(x)
}

.region_args <- function(chr, start, end, chr_map) {
  list(tchr = if (!is.null(chr)) normalize_chr(chr, chr_map) else NULL,
       start = if (is.null(start)) -Inf else as.numeric(start),
       end = if (is.null(end)) Inf else as.numeric(end))
}

#' Read genotypes for a region from a VCF file
#'
#' Only the requested region is kept in memory. If the VCF is
#' bgzip-compressed with a tabix index (`.tbi`/`.csi`) and the `tabix`
#' program (htslib) is on the PATH, the region is retrieved directly from the
#' index; otherwise the file is streamed in chunks and only lines inside the
#' region are kept. Plain and gzip-compressed VCFs both work.
#'
#' Genotype handling: the `GT` sub-field is taken from the FORMAT column
#' (FORMAT such as `GT:AD:DP:GQ` is supported); phased (`0|1`) and unphased
#' (`0/1`) calls are equivalent for dosages; missing calls (`.`, `./.`)
#' become `NA`; haploid calls are coded 0/2 and polyploid calls are scaled to
#' the 0-2 range; at multi-allelic sites all non-reference alleles count as
#' ALT and the site is flagged (removed by default in [genotype_qc()]).
#'
#' @param file Path to a VCF (`.vcf`, `.vcf.gz`).
#' @param chr,start,end Region to read (any chromosome naming style). If
#'   `chr` is `NULL` the whole file is read.
#' @param samples Optional sample names to keep.
#' @param chr_map Optional name mapping passed to [normalize_chr()].
#' @param method `"auto"` (tabix when possible, otherwise streaming),
#'   `"tabix"` or `"stream"`.
#' @param chunk_size Lines read per chunk when streaming.
#' @param every Keep only every `every`-th variant (thinning, e.g. for
#'   genome-wide distance trees); 1 keeps all.
#' @param verbose Print progress messages.
#' @return An object of class `qp_geno`: a list with `info` (variant table:
#'   `id`, `chr`, `pos`, `ref`, `alt`, `multiallelic`, `is_snp`), `dosage`
#'   (variants x samples matrix, `NA` = missing), `samples`, and for fully
#'   phased diploid data the haplotype matrices `hap1` and `hap2`.
#' @export
read_vcf <- function(file, chr = NULL, start = NULL, end = NULL, samples = NULL, chr_map = NULL,
                     method = c("auto", "tabix", "stream"), chunk_size = 10000L, every = 1L,
                     verbose = TRUE) {
  .check_file(file, "VCF")
  method <- match.arg(method)
  hdr <- .vcf_header(file)
  if (!length(hdr) || !startsWith(hdr[length(hdr)], "#CHROM")) {
    .stop("No '#CHROM' header line found in VCF '", file, "'.")
  }
  cols <- strsplit(hdr[length(hdr)], "\t", fixed = TRUE)[[1]]
  if (length(cols) < 10L) .stop("VCF '", file, "' has no sample (genotype) columns.")
  ra <- .region_args(chr, start, end, chr_map)
  t0 <- proc.time()[["elapsed"]]
  tabix <- Sys.which("tabix")
  has_index <- file.exists(paste0(file, ".tbi")) || file.exists(paste0(file, ".csi"))
  use_tabix <- method != "stream" && !is.null(ra$tchr) && has_index && nzchar(tabix) && every <= 1L
  if (method == "tabix" && !use_tabix) {
    .stop("method = 'tabix' needs a region, a .tbi/.csi index next to the bgzip-compressed VCF ",
          "and the 'tabix' program (htslib) on the PATH.")
  }
  if (use_tabix) {
    seqs <- system2(tabix, c("-l", shQuote(file)), stdout = TRUE)
    lab <- seqs[normalize_chr(seqs, chr_map) == ra$tchr]
    if (!length(lab)) {
      .stop("Chromosome '", chr, "' is not in the VCF index. VCF chromosomes: ",
            paste(utils::head(seqs, 20), collapse = ", "), ".")
    }
    reg <- sprintf("%s:%.0f-%.0f", lab[1], max(1, ra$start), min(ra$end, 2^31 - 1))
    body <- system2(tabix, c(shQuote(file), shQuote(reg)), stdout = TRUE)
    how <- "tabix index"
  } else {
    body <- .stream_lines(file, length(hdr), 1L, 2L, ra$tchr, ra$start, ra$end, chr_map, chunk_size, every)
    how <- if (every > 1L) sprintf("streaming (every %d)", as.integer(every)) else "streaming"
  }
  smp <- cols[-(1:9)]
  keep_s <- rep(TRUE, length(smp))
  if (!is.null(samples)) {
    miss <- setdiff(samples, smp)
    if (length(miss)) .warn(length(miss), " requested sample(s) not in the VCF: ", paste(utils::head(miss, 5), collapse = ", "))
    keep_s <- smp %in% samples
    if (!any(keep_s)) .stop("None of the requested samples are in the VCF.")
  }
  where <- if (is.null(ra$tchr)) "" else sprintf(" in chr %s:%s-%s", ra$tchr, .fmt(max(1, ra$start)),
                                                  if (is.finite(ra$end)) .fmt(ra$end) else "end")
  if (!length(body)) {
    .msg(verbose, "VCF: no variants", where, ".")
    return(.empty_geno(smp[keep_s], file))
  }
  dt <- data.table::fread(text = c(hdr[length(hdr)], body), sep = "\t", header = TRUE,
                          colClasses = "character", quote = "", data.table = FALSE,
                          showProgress = FALSE, check.names = FALSE)
  G <- as.matrix(dt[, 9L + which(keep_s), drop = FALSE])
  gt <- .extract_gt(G, dt[[9]])
  conv <- .gt_to_dosage(gt)
  info <- data.frame(id = dt[[3]], chr = normalize_chr(dt[[1]], chr_map),
                     pos = as.numeric(dt[[2]]), ref = dt[[4]], alt = dt[[5]],
                     stringsAsFactors = FALSE)
  noid <- info$id %in% c(".", "")
  info$id[noid] <- paste0(info$chr[noid], ":", info$pos[noid])
  info$multiallelic <- grepl(",", info$alt, fixed = TRUE)
  info$is_snp <- nchar(info$ref) == 1L & grepl("^[ACGTacgt](,[ACGTacgt])*$", info$alt)
  g <- .new_geno(info, conv$dosage, smp[keep_s], conv$hap1, conv$hap2, source = file, how = how)
  o <- order(g$info$pos)
  if (is.unsorted(g$info$pos)) g <- .geno_subset(g, o)
  .msg(verbose, "VCF: ", .fmt(nrow(info)), " variant(s) x ", .fmt(sum(keep_s)), " sample(s)", where,
       " read by ", how, " in ", round(proc.time()[["elapsed"]] - t0, 1), " s.")
  g
}

# Dosage of the second (ALT) HapMap allele for a set of calls.
.hmp_dosage <- function(calls, alleles) {
  al <- strsplit(alleles, "/", fixed = TRUE)[[1]]
  out <- rep(NA_real_, length(calls))
  if (length(al) < 2L) return(out)
  ref <- al[1]
  alt <- al[2]
  iupac <- c(R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC")
  one <- nchar(calls) == 1L
  if (any(one)) {
    c1 <- calls[one]
    d1 <- ifelse(c1 == ref, 0, ifelse(c1 == alt, 2, NA_real_))
    het <- c1 %in% names(iupac)
    if (any(het)) {
      d1[het] <- ifelse(iupac[c1[het]] == paste(sort(c(ref, alt)), collapse = ""), 1, NA_real_)
    }
    if (setequal(c(ref, alt), c("+", "-"))) d1[c1 == "0"] <- 1
    out[one] <- d1
  }
  two <- nchar(calls) == 2L
  if (any(two)) {
    a <- substr(calls[two], 1L, 1L)
    b <- substr(calls[two], 2L, 2L)
    ok <- a %in% c(ref, alt) & b %in% c(ref, alt)
    out[two] <- ifelse(ok, (a == alt) + (b == alt), NA_real_)
  }
  out
}

#' Read genotypes for a region from a TASSEL HapMap file
#'
#' Reads HapMap (`.hmp.txt`) genotype files as written by TASSEL and used by
#' GAPIT. Calls may be single IUPAC letters (`A`, `G`, `R` = A/G, `N` =
#' missing) or two letters (`AA`, `AG`, `NN`). The first allele in the
#' `alleles` column is taken as REF and the second as ALT.
#'
#' @inheritParams read_vcf
#' @param file Path to a HapMap file (optionally gzip-compressed).
#' @return A `qp_geno` object (see [read_vcf()]).
#' @export
read_hapmap <- function(file, chr = NULL, start = NULL, end = NULL, samples = NULL,
                        chr_map = NULL, chunk_size = 10000L, every = 1L, verbose = TRUE) {
  .check_file(file, "HapMap")
  lay <- .hapmap_layout(file)
  r <- lay$r
  ra <- .region_args(chr, start, end, chr_map)
  t0 <- proc.time()[["elapsed"]]
  body <- .stream_lines(file, 1L, r + 2L, r + 3L, ra$tchr, ra$start, ra$end, chr_map, chunk_size, every,
                        sep = lay$sep)
  smp <- lay$samples
  keep_s <- if (is.null(samples)) rep(TRUE, length(smp)) else smp %in% samples
  if (!any(keep_s)) .stop("None of the requested samples are in the HapMap file.")
  if (!length(body)) {
    .msg(verbose, "HapMap: no variants in the requested region.")
    return(.empty_geno(smp[keep_s], file))
  }
  dt <- data.table::fread(text = c(lay$header, body), sep = lay$sep, header = TRUE, colClasses = "character",
                          quote = if (lay$sep == "\t") "" else "\"", data.table = FALSE,
                          showProgress = FALSE, check.names = FALSE)
  G <- toupper(as.matrix(dt[, (lay$first_sample - 1L) + which(keep_s), drop = FALSE]))
  G[is.na(G)] <- "N"
  alleles <- toupper(dt[[r + 1L]])
  D <- matrix(NA_real_, nrow(G), ncol(G))
  for (al in unique(alleles)) {
    rows <- which(alleles == al)
    sub <- G[rows, , drop = FALSE]
    u <- unique(as.vector(sub))
    D[rows, ] <- .hmp_dosage(u, al)[match(sub, u)]
  }
  al_list <- strsplit(alleles, "/", fixed = TRUE)
  info <- data.frame(id = dt[[r]], chr = normalize_chr(dt[[r + 2L]], chr_map), pos = as.numeric(dt[[r + 3L]]),
                     ref = vapply(al_list, function(a) a[1], ""),
                     alt = vapply(al_list, function(a) if (length(a) > 1L) a[2] else NA_character_, ""),
                     stringsAsFactors = FALSE)
  info$multiallelic <- lengths(al_list) > 2L
  info$is_snp <- vapply(al_list, function(a) all(a %in% c("A", "C", "G", "T")), NA)
  g <- .new_geno(info, D, smp[keep_s], source = file, how = "streaming")
  if (is.unsorted(g$info$pos)) g <- .geno_subset(g, order(g$info$pos))
  .msg(verbose, "HapMap: ", .fmt(nrow(info)), " variant(s) x ", .fmt(sum(keep_s)), " sample(s) read in ",
       round(proc.time()[["elapsed"]] - t0, 1), " s.")
  g
}

#' Read genotypes from a VCF or HapMap file
#'
#' Chooses [read_vcf()] or [read_hapmap()] from the file name or content.
#' A `qp_geno` object is returned unchanged.
#' @param file Path to a VCF or HapMap file, or a `qp_geno` object.
#' @param ... Passed to [read_vcf()] or [read_hapmap()].
#' @return A `qp_geno` object.
#' @export
read_genotypes <- function(file, ...) {
  if (inherits(file, "qp_geno")) return(file)
  .check_file(file, "Genotype")
  low <- tolower(file)
  if (grepl("\\.vcf(\\.b?gz|\\.bz2|\\.xz)?$", low)) return(read_vcf(file, ...))
  if (grepl("\\.hmp(\\.txt|\\.csv)?(\\.gz)?$", low)) return(read_hapmap(file, ...))
  first <- readLines(file, n = 1L, warn = FALSE)
  if (startsWith(first, "##fileformat=VCF")) return(read_vcf(file, ...))
  if (.is_hapmap_file(file)) return(read_hapmap(file, ...))
  .stop("Cannot tell the genotype format of '", file, "'. Use a VCF (.vcf, .vcf.gz) or a TASSEL ",
        "HapMap file (tab- or comma-separated, header starting with rs#). For numeric (0/1/2) ",
        "genotype tables use as_genotypes().")
}

#' Create a genotype object from a dosage matrix
#'
#' For genotypes that are already numeric (for example GAPIT numeric format),
#' coded as 0/1/2 copies of the alternative allele.
#'
#' @param dosage Numeric matrix, variants in rows and samples in columns
#'   (0, 1, 2, `NA` for missing).
#' @param chr,pos Chromosome and position of each variant.
#' @param id Variant identifiers (default `chr:pos`).
#' @param ref,alt Optional allele labels.
#' @param samples Sample names (default column names).
#' @return A `qp_geno` object.
#' @export
as_genotypes <- function(dosage, chr, pos, id = NULL, ref = NA_character_, alt = NA_character_,
                         samples = colnames(dosage)) {
  dosage <- as.matrix(dosage)
  storage.mode(dosage) <- "double"
  n <- nrow(dosage)
  if (length(pos) != n) .stop("`pos` must have one value per row of `dosage`.")
  chr <- normalize_chr(rep_len(chr, n))
  if (is.null(id)) id <- paste0(chr, ":", pos)
  if (is.null(samples)) samples <- paste0("sample", seq_len(ncol(dosage)))
  info <- data.frame(id = as.character(id), chr = chr, pos = as.numeric(pos),
                     ref = rep_len(ref, n), alt = rep_len(alt, n), multiallelic = FALSE,
                     is_snp = TRUE, stringsAsFactors = FALSE)
  g <- .new_geno(info, dosage, samples, source = "matrix", how = "matrix")
  if (is.unsorted(g$info$pos)) g <- .geno_subset(g, order(g$info$pos))
  g
}

#' Genotype quality control
#'
#' Applies documented marker filters before LD calculation and reports how
#' many variants each step removed (the original script used an undefined
#' "variance" filter and treated missing calls as homozygous reference).
#'
#' Steps, in order: duplicated positions (first kept), multi-allelic sites
#' (if `biallelic_only`), indels (if `snps_only`), call rate (missing
#' fraction above `max_missing`), minor allele frequency below `min_maf`
#' (monomorphic sites are always removed).
#'
#' @param geno A `qp_geno` object.
#' @param min_maf Minimum minor allele frequency (default 0.05).
#' @param max_missing Maximum fraction of missing calls per variant.
#' @param biallelic_only Remove multi-allelic sites.
#' @param snps_only Remove indels and other non-SNP variants.
#' @param verbose Print the filter summary.
#' @return The filtered `qp_geno` object with `maf` and `missing` columns
#'   added to `info` and a `qc` table (step, removed, remaining).
#' @export
genotype_qc <- function(geno, min_maf = 0.05, max_missing = 0.2, biallelic_only = TRUE,
                        snps_only = FALSE, verbose = TRUE) {
  if (!inherits(geno, "qp_geno")) .stop("`geno` must be a qp_geno object (read_vcf(), read_hapmap()).")
  info <- geno$info
  D <- geno$dosage
  n <- ncol(D)
  ncall <- rowSums(!is.na(D))
  miss <- if (n) 1 - ncall / n else rep(1, nrow(D))
  af <- rowSums(D, na.rm = TRUE) / (2 * pmax(ncall, 1))
  maf <- pmin(af, 1 - af)
  keep <- rep(TRUE, nrow(info))
  log <- data.frame(step = "input", removed = 0L, remaining = nrow(info), stringsAsFactors = FALSE)
  step <- function(fail, label) {
    rem <- sum(keep & fail)
    keep <<- keep & !fail
    log <<- rbind(log, data.frame(step = label, removed = rem, remaining = sum(keep)))
  }
  step(duplicated(paste(info$chr, info$pos)), "duplicated position")
  if (biallelic_only) step(info$multiallelic, "multi-allelic")
  if (snps_only) step(!info$is_snp, "not a SNP")
  step(miss > max_missing, sprintf("missing rate > %s", format(max_missing)))
  step(ncall == 0 | maf < min_maf | maf <= 0, sprintf("MAF < %s or monomorphic", format(min_maf)))
  out <- .geno_subset(geno, keep)
  out$info$maf <- maf[keep]
  out$info$missing <- miss[keep]
  out$qc <- log
  .msg(verbose, "Genotype QC: ", .fmt(sum(keep)), " of ", .fmt(length(keep)), " variant(s) kept (",
       paste(sprintf("%s: -%s", log$step[-1], log$removed[-1]), collapse = "; "), ").")
  out
}
