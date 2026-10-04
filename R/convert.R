# HapMap -> VCF conversion and reference FASTA access ---------------------------

.iupac_pairs <- list(R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"), W = c("A", "T"),
                     K = c("G", "T"), M = c("A", "C"))

# VCF GT strings for HapMap calls, given the allele order (REF first).
.hmp_calls_to_gt <- function(calls, alleles) {
  out <- rep("./.", length(calls))
  ix <- function(b) match(b, alleles) - 1L
  one <- nchar(calls) == 1L
  if (any(one)) {
    c1 <- calls[one]
    r <- rep("./.", length(c1))
    i <- ix(c1)
    hom <- !is.na(i)
    r[hom] <- paste0(i[hom], "/", i[hom])
    het <- !hom & c1 %in% names(.iupac_pairs)
    if (any(het)) {
      pr <- .iupac_pairs[c1[het]]
      a <- vapply(pr, function(p) ix(p[1]), 0L)
      b <- vapply(pr, function(p) ix(p[2]), 0L)
      ok <- !is.na(a) & !is.na(b)
      rh <- r[het]
      rh[ok] <- paste0(pmin(a, b)[ok], "/", pmax(a, b)[ok])
      r[het] <- rh
    }
    if (length(alleles) >= 2L && setequal(alleles[1:2], c("+", "-"))) r[c1 == "0"] <- "0/1"
    out[one] <- r
  }
  two <- nchar(calls) == 2L
  if (any(two)) {
    a <- ix(substr(calls[two], 1L, 1L))
    b <- ix(substr(calls[two], 2L, 2L))
    ok <- !is.na(a) & !is.na(b)
    r <- rep("./.", sum(two))
    r[ok] <- paste0(pmin(a, b)[ok], "/", pmax(a, b)[ok])
    out[two] <- r
  }
  out
}

# Index of a FASTA file (read from .fai, or built by streaming the file).
.fasta_index <- function(fasta) {
  fai <- paste0(fasta, ".fai")
  if (file.exists(fai)) {
    d <- utils::read.table(fai, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "",
                           comment.char = "")
    d <- d[, 1:5]
    names(d) <- c("name", "length", "offset", "line_bases", "line_bytes")
    return(d)
  }
  if (.is_compressed(fasta)) {
    .stop("A compressed FASTA needs a .fai index (samtools faidx); or use the uncompressed FASTA.")
  }
  con <- file(fasta, "rb")
  head_raw <- readBin(con, "raw", 1e5)
  close(con)
  eol <- if (any(head_raw == as.raw(13))) 2 else 1
  con <- file(fasta, "rb")
  on.exit(close(con))
  nm <- character()
  len <- off <- lb <- lby <- numeric()
  cur <- 0L
  offset <- 0
  repeat {
    x <- readLines(con, n = 200000L, warn = FALSE)
    if (!length(x)) break
    nb <- nchar(x, type = "bytes") + eol
    st <- offset + cumsum(c(0, nb[-length(nb)]))
    h <- which(startsWith(x, ">"))
    seg_a <- c(1L, h)
    seg_b <- c(h - 1L, length(x))
    for (s in seq_along(seg_a)) {
      a <- seg_a[s]
      b <- seg_b[s]
      if (a > b) next
      if (startsWith(x[a], ">")) {
        cur <- cur + 1L
        nm[cur] <- sub("^>([^[:space:]]+).*$", "\\1", x[a])
        off[cur] <- st[a] + nb[a]
        len[cur] <- 0
        lb[cur] <- NA_real_
        lby[cur] <- NA_real_
        a <- a + 1L
        if (a > b) next
      }
      if (cur == 0L) next
      sl <- nchar(x[a:b])
      len[cur] <- len[cur] + sum(sl)
      if (is.na(lb[cur])) {
        lb[cur] <- sl[1]
        lby[cur] <- sl[1] + eol
      }
    }
    offset <- offset + sum(nb)
  }
  d <- data.frame(name = nm, length = len, offset = off, line_bases = lb, line_bytes = lby,
                  stringsAsFactors = FALSE)
  try(utils::write.table(d, fai, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE),
      silent = TRUE)
  d
}

# Reference bases at positions on one contig (upper case; NA outside).
.fasta_bases <- function(fasta, idx, contig, pos) {
  out <- rep(NA_character_, length(pos))
  r <- match(contig, idx$name)
  if (is.na(r)) return(out)
  ok <- which(!is.na(pos) & pos >= 1 & pos <= idx$length[r])
  if (!length(ok)) return(out)
  p <- pos[ok]
  byte <- idx$offset[r] + ((p - 1) %/% idx$line_bases[r]) * idx$line_bytes[r] + ((p - 1) %% idx$line_bases[r])
  o <- order(byte)
  bs <- byte[o]
  res <- character(length(bs))
  con <- file(fasta, "rb")
  on.exit(close(con))
  i <- 1L
  while (i <= length(bs)) {
    j <- findInterval(bs[i] + 1e6 - 1, bs)
    seek(con, bs[i])
    raw <- readBin(con, "raw", bs[j] - bs[i] + 1)
    res[i:j] <- toupper(rawToChar(raw[bs[i:j] - bs[i] + 1], multiple = TRUE))
    i <- j + 1L
  }
  tmp <- character(length(ok))
  tmp[o] <- res
  out[ok] <- tmp
  out
}

#' Convert a TASSEL HapMap file to VCF
#'
#' Streams a HapMap file in chunks (so very large files can be converted)
#' and writes a VCF 4.2 file with `GT` genotypes. Heterozygous IUPAC codes
#' (`R`, `Y`, ...), two-letter calls (`AG`) and missing calls (`N`, `NN`)
#' are converted to `0/1`, `0/0`, `./.` and so on; multi-allelic sites keep
#' all alleles.
#'
#' HapMap files do not say which allele is the reference-genome allele. If
#' `ref_fasta` (the reference genome) is given, REF is set to the reference
#' base and sites where the reference base matches none of the HapMap alleles
#' are flagged `RefMismatch` in the FILTER column. Without a FASTA, REF is the
#' first allele of the `alleles` column (noted in the VCF header).
#' Indels coded `+`/`-` have no sequence, so they are dropped by default
#' (`indels = "drop"`) or written with symbolic ALT alleles `<INS>`/`<DEL>`
#' (`indels = "symbolic"`).
#'
#' If `output` ends in `.gz` and the `bgzip`/`tabix` programs (htslib) are on
#' the PATH, the file is bgzip-compressed and indexed; otherwise it is gzip
#' compressed (readable by [read_vcf()], but not indexable).
#'
#' @param hapmap HapMap file (`.hmp.txt`, optionally gzip-compressed).
#' @param output Output VCF (`.vcf` or `.vcf.gz`).
#' @param ref_fasta Optional reference genome FASTA (uncompressed, or with a
#'   `.fai` index from `samtools faidx`).
#' @param chrom_sizes Optional chromosome lengths for the `##contig` header
#'   lines (anything accepted by [chrom_sizes()]); taken from the FASTA index
#'   when `ref_fasta` is given, otherwise from the species chosen with
#'   [set_species()].
#' @param chr,start,end Optional region to convert.
#' @param chr_prefix Optional prefix for chromosome names in the VCF (e.g.
#'   `"Chr"` to write `Chr1`); by default the FASTA names (if given) or the
#'   HapMap names are used.
#' @param indels `"drop"` or `"symbolic"` for `+`/`-` indels.
#' @param index Create a tabix index when possible.
#' @param chunk_size Lines processed per chunk.
#' @param chr_map Optional chromosome name mapping.
#' @param verbose Print a summary.
#' @return Invisibly, a list with the output path and counts of variants
#'   read, written, indels dropped and reference mismatches.
#' @examples
#' \donttest{
#' ex <- simulate_example_data(tempfile("qp"))
#' hapmap_to_vcf(ex$hapmap, tempfile(fileext = ".vcf"))
#' }
#' @export
hapmap_to_vcf <- function(hapmap, output, ref_fasta = NULL, chrom_sizes = NULL, chr = NULL,
                          start = NULL, end = NULL, chr_prefix = NULL, indels = c("drop", "symbolic"),
                          index = TRUE, chunk_size = 20000L, chr_map = NULL, verbose = TRUE) {
  indels <- match.arg(indels)
  .check_file(hapmap, "HapMap")
  if (!grepl("\\.vcf(\\.gz)?$", output, ignore.case = TRUE)) .stop("`output` must end in .vcf or .vcf.gz.")
  .ensure_dir(output)
  lay <- .hapmap_layout(hapmap)
  hdr <- lay$header
  rr <- lay$r
  fs <- lay$first_sample
  samples <- lay$samples
  fai <- NULL
  fasta_norm <- NULL
  if (!is.null(ref_fasta)) {
    .check_file(ref_fasta, "Reference FASTA")
    fai <- .fasta_index(ref_fasta)
    fasta_norm <- normalize_chr(fai$name, chr_map)
  }
  label_of <- function(chr_raw, chr_n) {
    if (!is.null(chr_prefix)) return(paste0(chr_prefix, chr_n))
    if (!is.null(fai)) {
      lab <- fai$name[match(chr_n, fasta_norm)]
      return(ifelse(is.na(lab), chr_raw, lab))
    }
    chr_raw
  }
  ra <- .region_args(chr, start, end, chr_map)
  tmp <- tempfile(fileext = ".vcfbody")
  out_con <- file(tmp, "w")
  in_con <- file(hapmap, "r")
  readLines(in_con, n = 1L, warn = FALSE)
  n_read <- n_written <- n_indel <- n_mismatch <- 0
  seen <- character()
  last_chr <- NA_character_
  last_pos <- -Inf
  unsorted <- FALSE
  t0 <- proc.time()[["elapsed"]]
  repeat {
    x <- readLines(in_con, n = chunk_size, warn = FALSE)
    if (!length(x)) break
    x <- x[nzchar(x)]
    if (!is.null(ra$tchr) && length(x)) {
      cn <- normalize_chr(.field(x, rr + 2L, lay$sep), chr_map)
      pp <- suppressWarnings(as.numeric(.field(x, rr + 3L, lay$sep)))
      x <- x[cn == ra$tchr & !is.na(pp) & pp >= ra$start & pp <= ra$end]
    }
    if (!length(x)) next
    dt <- data.table::fread(text = c(hdr, x), sep = lay$sep, header = TRUE, colClasses = "character",
                            quote = if (lay$sep == "\t") "" else "\"", data.table = FALSE,
                            showProgress = FALSE, check.names = FALSE)
    n_read <- n_read + nrow(dt)
    al <- toupper(dt[[rr + 1L]])
    chr_raw <- dt[[rr + 2L]]
    chr_n <- normalize_chr(chr_raw, chr_map)
    pos <- suppressWarnings(as.numeric(dt[[rr + 3L]]))
    alist <- strsplit(al, "/", fixed = TRUE)
    plusminus <- vapply(alist, function(a) any(a %in% c("+", "-", "0")), NA)
    seqindel <- vapply(alist, function(a) any(nchar(a) > 1L), NA)
    keep <- !is.na(pos) & lengths(alist) >= 1L
    if (indels == "drop") {
      n_indel <- n_indel + sum(keep & plusminus)
      keep <- keep & !plusminus
    }
    if (!any(keep)) next
    dt <- dt[keep, , drop = FALSE]
    al <- al[keep]
    alist <- alist[keep]
    chr_raw <- chr_raw[keep]
    chr_n <- chr_n[keep]
    pos <- pos[keep]
    plusminus <- plusminus[keep]
    seqindel <- seqindel[keep]
    n <- nrow(dt)
    refb <- rep(NA_character_, n)
    if (!is.null(fai)) {
      for (cc in unique(chr_n)) {
        rows <- which(chr_n == cc)
        contig <- fai$name[match(cc, fasta_norm)]
        if (!is.na(contig)) refb[rows] <- .fasta_bases(ref_fasta, fai, contig, pos[rows])
      }
    }
    snp_like <- !plusminus & !seqindel
    use_ref <- !is.na(refb) & snp_like & vapply(seq_len(n), function(i) refb[i] %in% alist[[i]], NA)
    mism <- !is.na(refb) & snp_like & !use_ref
    n_mismatch <- n_mismatch + sum(mism)
    key <- paste(al, ifelse(use_ref, refb, ""), sep = "|")
    G <- toupper(as.matrix(dt[, fs:ncol(dt), drop = FALSE]))
    G[is.na(G)] <- "N"
    GT <- matrix("./.", n, ncol(G))
    REF <- ALT <- character(n)
    for (k in unique(key)) {
      rows <- which(key == k)
      a <- alist[[rows[1]]]
      ord <- if (use_ref[rows[1]]) c(refb[rows[1]], setdiff(a, refb[rows[1]])) else a
      sub <- G[rows, , drop = FALSE]
      u <- unique(as.vector(sub))
      GT[rows, ] <- .hmp_calls_to_gt(u, ord)[match(sub, u)]
      if (plusminus[rows[1]]) {
        REF[rows] <- ifelse(!is.na(refb[rows]), refb[rows], "N")
        alt_lab <- ifelse(ord[-1] == "+", "<INS>", ifelse(ord[-1] == "-", "<DEL>", "<INDEL>"))
        ALT[rows] <- if (length(alt_lab)) paste(alt_lab, collapse = ",") else "."
      } else {
        REF[rows] <- ord[1]
        ALT[rows] <- if (length(ord) > 1L) paste(ord[-1], collapse = ",") else "."
      }
    }
    lab <- label_of(chr_raw, chr_n)
    seen <- union(seen, unique(lab))
    same <- c(identical(chr_n[1], last_chr), chr_n[-1] == chr_n[-n])
    prev <- c(last_pos, pos[-n])
    if (any(same & pos < prev)) unsorted <- TRUE
    last_chr <- chr_n[n]
    last_pos <- pos[n]
    id <- dt[[rr]]
    id[is.na(id) | !nzchar(id)] <- "."
    fixed <- paste(lab, format(pos, scientific = FALSE, trim = TRUE), id, REF, ALT, ".",
                   ifelse(mism, "RefMismatch", "PASS"), ".", "GT", sep = "\t")
    body <- do.call(paste, c(list(fixed), as.data.frame(GT, stringsAsFactors = FALSE), sep = "\t"))
    writeLines(body, out_con)
    n_written <- n_written + n
  }
  close(in_con)
  close(out_con)
  on.exit(unlink(tmp), add = TRUE)

  cs <- NULL
  if (is.null(chrom_sizes) && is.null(fai)) chrom_sizes <- .species_sizes()
  if (!is.null(chrom_sizes)) {
    cs <- chrom_sizes(chrom_sizes, chr_map = chr_map)
    cs$label <- label_of(cs$chr, cs$chr)
  } else if (!is.null(fai)) {
    cs <- data.frame(label = fai$name, length = fai$length, stringsAsFactors = FALSE)
  }
  contig_lines <- if (!is.null(cs)) {
    cs <- cs[cs$label %in% seen | !length(seen), , drop = FALSE]
    sprintf("##contig=<ID=%s,length=%.0f>", cs$label, cs$length)
  } else {
    sprintf("##contig=<ID=%s>", seen)
  }
  header <- c("##fileformat=VCFv4.2",
              paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
              paste0("##source=quardplot::hapmap_to_vcf from ", basename(hapmap)),
              if (!is.null(ref_fasta)) paste0("##reference=file://", normalizePath(ref_fasta)) else
                "##quardplot_note=No reference genome given: REF/ALT follow the order of the HapMap alleles column and may not match the reference-genome base",
              contig_lines,
              if (!is.null(ref_fasta)) "##FILTER=<ID=RefMismatch,Description=\"Reference base matches none of the HapMap alleles\">",
              if (indels == "symbolic") c("##ALT=<ID=INS,Description=\"Insertion (HapMap '+')\">",
                                          "##ALT=<ID=DEL,Description=\"Deletion (HapMap '-')\">"),
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
              paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
                    collapse = "\t"))
  gz <- grepl("\\.gz$", output, ignore.case = TRUE)
  plain <- if (gz) sub("\\.gz$", "", output, ignore.case = TRUE) else output
  writeLines(header, plain)
  file.append(plain, tmp)
  how <- "plain text"
  if (gz) {
    bgzip <- Sys.which("bgzip")
    if (nzchar(bgzip)) {
      system2(bgzip, c("-f", shQuote(plain)))
      if (!identical(paste0(plain, ".gz"), output)) file.rename(paste0(plain, ".gz"), output)
      how <- "bgzip"
      tabix <- Sys.which("tabix")
      if (index && nzchar(tabix) && !unsorted) {
        system2(tabix, c("-f", "-p", "vcf", shQuote(output)))
        how <- "bgzip + tabix index"
      }
    } else {
      inp <- file(plain, "r")
      gzc <- gzfile(output, "w")
      repeat {
        l <- readLines(inp, n = 50000L)
        if (!length(l)) break
        writeLines(l, gzc)
      }
      close(inp)
      close(gzc)
      unlink(plain)
      how <- "gzip (install htslib's bgzip/tabix to make an indexable file)"
    }
  }
  if (unsorted) .warn("HapMap positions are not sorted; sort the VCF (e.g. bcftools sort) before indexing.")
  .msg(verbose, "HapMap -> VCF: ", .fmt(n_written), " of ", .fmt(n_read), " variant(s) written to ", output,
       " (", how, ")", if (n_indel) paste0("; ", .fmt(n_indel), " +/- indel(s) dropped") else "",
       if (!is.null(ref_fasta)) paste0("; ", .fmt(n_mismatch), " REF mismatch(es)") else "",
       " in ", round(proc.time()[["elapsed"]] - t0, 1), " s.")
  invisible(list(file = output, read = n_read, written = n_written, indels_dropped = n_indel,
                 ref_mismatch = n_mismatch))
}
