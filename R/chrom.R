# Chromosome names, sizes and genome layout --------------------------------

#' Normalise chromosome names
#'
#' Removes common prefixes (`chr`, `Chr`, `CHR`, `chromosome`, `ch`, and
#' tomato-style `SL4.0ch`), leading zeros (`chr01` -> `1`) and unifies
#' sex/organelle chromosome names, so that GWAS, GFF, VCF and block files
#' using different conventions can be matched.
#'
#' @param x Character or numeric vector of chromosome names.
#' @param chr_map Optional named vector for names that cannot be normalised
#'   automatically, e.g. `c(NC_029256.1 = "1")`. Names are the original
#'   labels, values the chromosome names to use.
#' @return Character vector of normalised names (`"1"`, `"2"`, ..., `"X"`,
#'   `"MT"`, `"Pt"`).
#' @examples
#' normalize_chr(c("Chr1", "chr01", "1", "CHR_X", "ChrM", "ChrC", "SL4.0ch01"))
#' @export
normalize_chr <- function(x, chr_map = NULL) {
  x <- trimws(as.character(x))
  if (!is.null(chr_map)) {
    if (is.null(names(chr_map))) {
      .stop("`chr_map` must be a named vector, e.g. c(NC_029256.1 = '1').")
    }
    hit <- !is.na(x) & x %in% names(chr_map)
    x[hit] <- as.character(chr_map[x[hit]])
  }
  y <- sub("^SL[0-9.]+(?=ch)", "", x, perl = TRUE, ignore.case = TRUE)        # tomato SL4.0ch01
  y <- sub("^(chromosome|chrom|chr|ch(?=[0-9]))[_ .-]?", "", y, perl = TRUE, ignore.case = TRUE)
  empty <- !is.na(y) & !nzchar(y)
  y[empty] <- x[empty]
  y <- sub("^0+(?=[0-9])", "", y, perl = TRUE)
  up <- toupper(y)
  sex <- !is.na(up) & up %in% c("X", "Y", "W", "Z")
  y[sex] <- up[sex]
  y[!is.na(up) & up %in% c("M", "MT", "MITO", "MITOCHONDRIA", "MITOCHONDRION")] <- "MT"
  y[!is.na(up) & up %in% c("C", "PT", "CP", "CHLOROPLAST", "PLASTID")] <- "Pt"
  y[is.na(x)] <- NA_character_
  y
}

#' Sort chromosome names in natural order
#'
#' Numbered chromosomes first (numerically), then X/Y/W/Z, other names
#' alphabetically, and organelles (MT, Pt) last.
#' @param x Chromosome names (normalised with [normalize_chr()]).
#' @return Unique chromosome names in natural order.
#' @export
sort_chr <- function(x) {
  u <- unique(as.character(x[!is.na(x)]))
  isnum <- grepl("^[0-9]+$", u)
  num <- rep(Inf, length(u))
  num[isnum] <- as.numeric(u[isnum])
  grp <- ifelse(isnum, 1L, ifelse(u %in% c("X", "Y", "W", "Z"), 2L,
                                  ifelse(u %in% c("MT", "Pt"), 4L, 3L)))
  u[order(grp, num, u)]
}

# Built-in chromosome lengths (bp). Verify against your own reference
# (.fai) when exact lengths matter; any file or named vector can be used.
.qp_presets <- list(
  # Oryza sativa Nipponbare IRGSP-1.0 / MSU7 pseudomolecules
  rice = c(`1` = 43270923, `2` = 35937250, `3` = 36413819, `4` = 35502694,
           `5` = 29958434, `6` = 31248787, `7` = 29697621, `8` = 28443022,
           `9` = 23012720, `10` = 23207287, `11` = 29021106, `12` = 27531856),
  # Arabidopsis thaliana TAIR10
  arabidopsis = c(`1` = 30427671, `2` = 19698289, `3` = 23459830,
                  `4` = 18585056, `5` = 26975502),
  # Solanum lycopersicum Heinz 1706 SL4.0 (ITAG4.0 / ITAG4.1); chromosome 0
  # (unanchored contigs) is not included
  tomato = c(`1` = 90863682, `2` = 53473368, `3` = 65298490, `4` = 64459972,
             `5` = 65269487, `6` = 47258699, `7` = 67883646, `8` = 63995357,
             `9` = 68513564, `10` = 64792705, `11` = 54379777, `12` = 66688036),
  # Homo sapiens GRCh38 (hg38)
  human_grch38 = c(`1` = 248956422, `2` = 242193529, `3` = 198295559,
                   `4` = 190214555, `5` = 181538259, `6` = 170805979,
                   `7` = 159345973, `8` = 145138636, `9` = 138394717,
                   `10` = 133797422, `11` = 135086622, `12` = 133275309,
                   `13` = 114364328, `14` = 107043718, `15` = 101991189,
                   `16` = 90338345, `17` = 83257441, `18` = 80373285,
                   `19` = 58617616, `20` = 64444167, `21` = 46709983,
                   `22` = 50818468, X = 156040895, Y = 57227415),
  # Homo sapiens GRCh37 (hg19)
  human_grch37 = c(`1` = 249250621, `2` = 243199373, `3` = 198022430,
                   `4` = 191154276, `5` = 180915260, `6` = 171115067,
                   `7` = 159138663, `8` = 146364022, `9` = 141213431,
                   `10` = 135534747, `11` = 135006516, `12` = 133851895,
                   `13` = 115169878, `14` = 107349540, `15` = 102531392,
                   `16` = 90354753, `17` = 81195210, `18` = 78077248,
                   `19` = 59128983, `20` = 63025520, `21` = 48129895,
                   `22` = 51304566, X = 155270560, Y = 59373566),
  # Mus musculus GRCm39 (mm39)
  mouse_grcm39 = c(`1` = 195154279, `2` = 181755017, `3` = 159745316,
                   `4` = 156860686, `5` = 151758149, `6` = 149588044,
                   `7` = 144995196, `8` = 130127694, `9` = 124359700,
                   `10` = 130530862, `11` = 121973369, `12` = 120092757,
                   `13` = 120883175, `14` = 125139656, `15` = 104073951,
                   `16` = 98008968, `17` = 95294699, `18` = 90720763,
                   `19` = 61420004, X = 169476592, Y = 91455967),
  # Mus musculus GRCm38 (mm10)
  mouse_grcm38 = c(`1` = 195471971, `2` = 182113224, `3` = 160039680,
                   `4` = 156508116, `5` = 151834684, `6` = 149736546,
                   `7` = 145441459, `8` = 129401213, `9` = 124595110,
                   `10` = 130694993, `11` = 122082543, `12` = 120129022,
                   `13` = 120421639, `14` = 124902244, `15` = 104043685,
                   `16` = 98207768, `17` = 94987271, `18` = 90702639,
                   `19` = 61431566, X = 171031299, Y = 91744698)
)

.qp_preset_alias <- c(
  rice = "rice", irgsp1 = "rice", "irgsp-1.0" = "rice", rice_irgsp1 = "rice",
  msu7 = "rice", rice_msu7 = "rice", oryza_sativa = "rice",
  arabidopsis = "arabidopsis", tair10 = "arabidopsis",
  arabidopsis_tair10 = "arabidopsis", athaliana = "arabidopsis",
  arabidopsis_thaliana = "arabidopsis",
  tomato = "tomato", sl4.0 = "tomato", sl4 = "tomato", itag4 = "tomato", tomato_sl4 = "tomato",
  solanum_lycopersicum = "tomato", slycopersicum = "tomato",
  oryza_sativa_japonica = "rice", homo_sapiens = "human_grch38", mus_musculus = "mouse_grcm39",
  human = "human_grch38", grch38 = "human_grch38", hg38 = "human_grch38",
  human_grch38 = "human_grch38", grch37 = "human_grch37", hg19 = "human_grch37",
  human_grch37 = "human_grch37", mouse = "mouse_grcm39", grcm39 = "mouse_grcm39",
  mm39 = "mouse_grcm39", mouse_grcm39 = "mouse_grcm39", grcm38 = "mouse_grcm38",
  mm10 = "mouse_grcm38", mouse_grcm38 = "mouse_grcm38"
)

#' Chromosome lengths
#'
#' Returns chromosome lengths from a built-in preset, a file or a vector, or
#' infers them from GWAS positions. This replaces the chromosome lengths that
#' had to be edited inside the original script.
#'
#' @param x One of:
#'   * a preset name: `"rice"` (IRGSP-1.0/MSU7), `"arabidopsis"` (TAIR10),
#'     `"tomato"` (SL4.0), `"human_grch38"`, `"human_grch37"`,
#'     `"mouse_grcm39"`, `"mouse_grcm38"` (aliases such as `"human"`,
#'     `"hg19"`, `"mm10"` and `"tair10"` also work; see [list_species()]);
#'   * a file: a FASTA index (`.fai`), the reference genome FASTA itself
#'     (lengths are counted, and a `.fai` index is saved next to an
#'     uncompressed FASTA), a GFF3 file with `##sequence-region` lines or
#'     `chromosome`/`region` records, a two-column chromosome-size file
#'     (name, length; header optional) or a VCF whose header has
#'     `##contig=<ID=...,length=...>` lines;
#'   * a named numeric vector (`c(Chr1 = 30427671, ...)`);
#'   * a data frame whose first two columns are name and length;
#'   * `NULL`: the lengths of the species chosen with [set_species()], or,
#'     when `gwas` is given, the largest GWAS position per chromosome.
#' @param gwas GWAS data from [read_gwas()]; used when `x` is `NULL`.
#' @param chr_map Optional name mapping passed to [normalize_chr()].
#' @return A data frame with columns `chr` and `length`, in natural order.
#' @examples
#' chrom_sizes("rice")
#' chrom_sizes(c(Chr1 = 1000, Chr2 = 800))
#' @export
chrom_sizes <- function(x = NULL, gwas = NULL, chr_map = NULL) {
  src <- NULL
  if (is.null(x)) {
    if (is.null(gwas)) {
      sp <- .species_sizes()
      if (!is.null(sp)) return(sp)
      .stop("No chromosome lengths: choose the species first (set_species(\"rice\"); see list_species()), ",
            "or give a preset name, a .fai/FASTA/size file, a VCF or a named vector.")
    }
    out <- stats::aggregate(gwas$pos, by = list(chr = gwas$chr), FUN = max)
    names(out) <- c("chr", "length")
    src <- "inferred from the largest GWAS position per chromosome"
  } else if (is.data.frame(x)) {
    if (ncol(x) < 2) .stop("A chromosome-size data frame needs name and length columns.")
    nm <- names(x)
    cc <- .find_col(nm, c("chr", "chrom", "chromosome", "name", "seqname", "seqid"))
    lc <- .find_col(nm, c("length", "len", "size", "bp"))
    if (is.na(cc)) cc <- nm[1]
    if (is.na(lc)) lc <- nm[2]
    out <- data.frame(chr = x[[cc]], length = x[[lc]], stringsAsFactors = FALSE)
    src <- "data frame"
  } else if (is.numeric(x)) {
    if (is.null(names(x))) .stop("A numeric `chrom_sizes` vector must be named by chromosome.")
    out <- data.frame(chr = names(x), length = as.numeric(x), stringsAsFactors = FALSE)
    src <- "named vector"
  } else if (is.character(x) && length(x) == 1L) {
    key <- tolower(x)
    if (key %in% names(.qp_preset_alias)) {
      v <- .qp_presets[[.qp_preset_alias[[key]]]]
      out <- data.frame(chr = names(v), length = as.numeric(v), stringsAsFactors = FALSE)
      src <- paste0("preset '", .qp_preset_alias[[key]], "'")
    } else if (file.exists(x)) {
      out <- .read_chrom_file(x)
      src <- paste0("file '", basename(x), "'")
    } else {
      .stop("`chrom_sizes` = '", x, "' is neither a file nor a built-in genome. Built-in genomes: ",
            paste(names(.qp_presets), collapse = ", "), " (see list_species()).")
    }
  } else {
    .stop("Unsupported `chrom_sizes` value.")
  }
  out$chr <- normalize_chr(out$chr, chr_map)
  out$length <- suppressWarnings(as.numeric(out$length))
  out <- out[!is.na(out$chr) & is.finite(out$length) & out$length > 0, , drop = FALSE]
  out <- out[!duplicated(out$chr), , drop = FALSE]
  if (!nrow(out)) .stop("No valid chromosome lengths found.")
  out <- out[match(sort_chr(out$chr), out$chr), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "source") <- src
  out
}

# Read a .fai / chrom.sizes / FASTA / GFF / VCF header into chr + length.
.read_chrom_file <- function(path) {
  .check_file(path, "Chromosome-size")
  first <- readLines(path, n = 1L, warn = FALSE)
  if (length(first) && startsWith(first, ">")) return(.fasta_lengths(path))
  if (length(first) && (startsWith(first, "##gff-version") || startsWith(first, "#gtf") ||
                        grepl("\\.(gff3?|gtf)(\\.gz)?$", path, ignore.case = TRUE))) {
    return(.gff_lengths(path))
  }
  if (length(first) && startsWith(first, "##")) {
    hdr <- .vcf_header(path)
    ctg <- hdr[startsWith(hdr, "##contig=<")]
    if (!length(ctg)) .stop("VCF '", path, "' has no ##contig lines with lengths.")
    id <- sub(".*[<,]ID=([^,>]+).*", "\\1", ctg)
    len <- ifelse(grepl("[<,]length=[0-9]+", ctg), sub(".*[<,]length=([0-9]+).*", "\\1", ctg), NA)
    if (all(is.na(len))) {
      .stop("The ##contig lines of '", basename(path), "' have no lengths. Choose the species with ",
            "set_species() (see list_species()), or give a .fai file or a named vector of lengths.")
    }
    return(data.frame(chr = id, length = as.numeric(len), stringsAsFactors = FALSE))
  }
  d <- .fread_df(path, header = "auto")
  if (ncol(d) < 2) .stop("Chromosome-size file '", path, "' needs at least two columns (name, length).")
  data.frame(chr = d[[1]], length = d[[2]], stringsAsFactors = FALSE)
}

# Sequence lengths of a FASTA file (from its .fai index when present).
.fasta_lengths <- function(path) {
  fai <- paste0(path, ".fai")
  if (file.exists(fai)) {
    d <- .fread_df(fai, header = FALSE)
    return(data.frame(chr = as.character(d[[1]]), length = as.numeric(d[[2]]), stringsAsFactors = FALSE))
  }
  if (!.is_compressed(path)) {
    d <- .fasta_index(path)
    return(data.frame(chr = d$name, length = d$length, stringsAsFactors = FALSE))
  }
  con <- file(path, "r")
  on.exit(close(con))
  nm <- character()
  len <- numeric()
  cur <- 0L
  repeat {
    x <- readLines(con, n = 500000L, warn = FALSE)
    if (!length(x)) break
    h <- startsWith(x, ">")
    g <- cur + cumsum(h)
    if (any(h)) {
      nm <- c(nm, sub("^>([^[:space:]]+).*$", "\\1", x[h]))
      len <- c(len, numeric(sum(h)))
      cur <- cur + sum(h)
    }
    keep <- !h & g > 0L
    if (any(keep)) {
      s <- rowsum(nchar(sub("\r$", "", x[keep]), type = "bytes"), g[keep])
      k <- as.integer(rownames(s))
      len[k] <- len[k] + s[, 1]
    }
  }
  if (!length(nm)) .stop("No sequences ('>' lines) found in FASTA file '", basename(path), "'.")
  data.frame(chr = nm, length = len, stringsAsFactors = FALSE)
}

# Sequence lengths from a GFF3/GTF file: ##sequence-region lines, otherwise
# records of type chromosome / region / contig / scaffold.
.gff_lengths <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  sr <- character()
  hit <- character()
  repeat {
    x <- readLines(con, n = 200000L, warn = FALSE)
    if (!length(x)) break
    sr <- c(sr, x[startsWith(x, "##sequence-region")])
    if (!length(sr)) {
      hit <- c(hit, x[grepl("^[^#][^\t]*\t[^\t]*\t(chromosome|region|contig|scaffold|supercontig)\t", x)])
    } else if (!any(startsWith(x, "#"))) {
      break
    }
  }
  if (length(sr)) {
    p <- strsplit(trimws(sub("^##sequence-region", "", sr)), "[[:space:]]+")
    d <- data.frame(chr = vapply(p, `[`, "", 1L), length = suppressWarnings(as.numeric(vapply(p, `[`, "", 3L))),
                    stringsAsFactors = FALSE)
    return(d[!duplicated(d$chr), , drop = FALSE])
  }
  if (length(hit)) {
    f <- strsplit(hit, "\t", fixed = FALSE)
    d <- data.frame(chr = vapply(f, `[`, "", 1L), start = as.numeric(vapply(f, `[`, "", 4L)),
                    length = as.numeric(vapply(f, `[`, "", 5L)), stringsAsFactors = FALSE)
    d <- d[d$start == 1, , drop = FALSE]
    d <- d[order(-d$length), , drop = FALSE]
    d <- d[!duplicated(d$chr), c("chr", "length"), drop = FALSE]
    if (nrow(d)) return(d)
  }
  .stop("The GFF file '", basename(path), "' does not contain chromosome lengths (no ##sequence-region ",
        "lines or chromosome records). Use a species preset (list_species()), the genome FASTA or its ",
        ".fai index, or a named vector of lengths.")
}

# Header lines of a VCF (plain or gzip/bgzip), up to and including #CHROM.
.vcf_header <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  hdr <- character()
  repeat {
    l <- readLines(con, n = 2000L, warn = FALSE)
    if (!length(l)) break
    k <- which(startsWith(l, "#CHROM"))
    if (length(k)) {
      hdr <- c(hdr, l[seq_len(k[1])])
      break
    }
    hdr <- c(hdr, l)
    if (!all(startsWith(l, "#"))) break
  }
  hdr
}

#' Genome layout for genome-wide plots
#'
#' Computes cumulative offsets of chromosomes from their lengths, replacing
#' the hard-coded offset vector of the original script.
#'
#' @param sizes Output of [chrom_sizes()].
#' @param gap Gap between chromosomes as a fraction of the genome length.
#' @return Data frame with `chr`, `length`, `offset`, `center`, `end_cum`
#'   and `index`.
#' @export
genome_layout <- function(sizes, gap = 0.004) {
  if (!all(c("chr", "length") %in% names(sizes))) sizes <- chrom_sizes(sizes)
  g <- gap * sum(sizes$length)
  off <- c(0, cumsum(sizes$length + g)[-nrow(sizes)])
  data.frame(chr = sizes$chr, length = sizes$length, offset = off,
             center = off + sizes$length / 2, end_cum = off + sizes$length,
             index = seq_len(nrow(sizes)), stringsAsFactors = FALSE)
}
