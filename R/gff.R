# Gene annotation (GFF3 / GTF) ---------------------------------------------

.gene_types_default <- c("gene", "pseudogene", "ncRNA_gene", "transposable_element_gene",
                         "protein_coding_gene", "rRNA_gene", "tRNA_gene", "snRNA_gene",
                         "snoRNA_gene", "miRNA_gene")
.sub_types <- c("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR", "5UTR", "3UTR",
                "five_prime_utr", "three_prime_utr", "UTR5", "UTR3", "noncoding_exon",
                "pseudogenic_exon")
.skip_types <- c("start_codon", "stop_codon", "intron", "chromosome", "region", "contig",
                 "scaffold", "supercontig", "biological_region", "Selenocysteine",
                 "repeat_region", "match", "match_part", "polypeptide", "protein")

.sub_class <- function(type) {
  out <- rep("exon", length(type))
  out[type == "CDS"] <- "CDS"
  out[type %in% c("five_prime_UTR", "5UTR", "five_prime_utr", "UTR5")] <- "UTR5"
  out[type %in% c("three_prime_UTR", "3UTR", "three_prime_utr", "UTR3")] <- "UTR3"
  out[type == "UTR"] <- "UTR"
  out
}

# Value of one GFF3 attribute (key=value;...), URL-decoded. The key must
# start the attribute string or follow a ';', so "ID" never matches
# "Locus_ID", and a key that is the last attribute (no trailing ';') works.
.attr_gff3 <- function(a, key) {
  out <- rep(NA_character_, length(a))
  m <- regexpr(paste0("(^|;)[[:space:]]*", key, "=[^;]*"), a, perl = TRUE)
  hit <- m > 0
  if (any(hit)) {
    v <- regmatches(a, m)
    out[hit] <- trimws(sub(paste0("^;?[[:space:]]*", key, "="), "", v, perl = TRUE))
  }
  enc <- !is.na(out) & grepl("%[0-9A-Fa-f]{2}", out)
  if (any(enc)) out[enc] <- vapply(out[enc], utils::URLdecode, "", USE.NAMES = FALSE)
  out[!is.na(out) & !nzchar(out)] <- NA_character_
  out
}

# Value of one GTF attribute (key "value"; ...).
.attr_gtf <- function(a, key) {
  out <- rep(NA_character_, length(a))
  m <- regexpr(paste0("(^|;)[[:space:]]*", key, "[[:space:]]+\"?[^\";]*\"?"), a, perl = TRUE)
  hit <- m > 0
  if (any(hit)) {
    v <- regmatches(a, m)
    v <- sub(paste0("^;?[[:space:]]*", key, "[[:space:]]+\"?"), "", v, perl = TRUE)
    out[hit] <- sub("\"$", "", v)
  }
  out
}

.merge_intervals <- function(s, e) {
  o <- order(s, e)
  s <- s[o]
  e <- e[o]
  rs <- s[1]
  re <- e[1]
  out_s <- numeric()
  out_e <- numeric()
  if (length(s) > 1L) {
    for (k in 2:length(s)) {
      if (s[k] <= re + 1) {
        re <- max(re, e[k])
      } else {
        out_s <- c(out_s, rs)
        out_e <- c(out_e, re)
        rs <- s[k]
        re <- e[k]
      }
    }
  }
  data.frame(start = c(out_s, rs), end = c(out_e, re))
}

# Collapse rows sharing `key` into one range (first chr/seqid/strand).
.agg_ranges <- function(df, key) {
  k <- df[[key]]
  first <- df[!duplicated(k), , drop = FALSE]
  data.frame(id = first[[key]], chr = first$chr, seqid = first$seqid, strand = first$strand,
             start = as.numeric(tapply(df$start, k, min)[first[[key]]]),
             end = as.numeric(tapply(df$end, k, max)[first[[key]]]),
             stringsAsFactors = FALSE)
}

#' Read a GFF3 or GTF gene annotation
#'
#' Parses gene models by following the `ID`/`Parent` links
#' (gene -> mRNA/transcript -> exon/CDS/UTR). Works with MSU, RAP-DB,
#' TAIR, Ensembl and NCBI style files, with GTF (`gene_id`/`transcript_id`),
#' with genes that have several transcripts, genes without transcript
#' features, transcripts without exon features (exons are then derived from
#' CDS and UTR) and overlapping genes. Attribute values are URL-decoded.
#'
#' For each gene a representative ("canonical") transcript is chosen: the one
#' with the longest coding sequence, then the longest exonic length, then the
#' first in the file. The original lines of the file are kept so that the
#' records of the genes in a region can be written out unchanged with
#' [write_region_genes()].
#'
#' @param file Path to a GFF3 or GTF file (optionally compressed).
#' @param chr Optional chromosome to keep (any naming style); reading only one
#'   chromosome is much faster for large annotations.
#' @param start,end Optional region on `chr`; genes overlapping it are kept.
#' @param chr_map Optional name mapping passed to [normalize_chr()].
#' @param gene_types Feature types treated as genes.
#' @param verbose Print progress messages.
#' @return An object of class `qp_genes`: a list with data frames `genes`,
#'   `transcripts`, `features` (exon/CDS/UTR parts) and `raw` (original
#'   lines), plus the file format.
#' @export
read_gff <- function(file, chr = NULL, start = NULL, end = NULL, chr_map = NULL,
                     gene_types = .gene_types_default, verbose = TRUE) {
  .check_file(file, "GFF")
  lines <- .gff_feature_lines(file)
  .parse_gff_lines(lines, file, chr = chr, start = start, end = end, chr_map = chr_map,
                   gene_types = gene_types, verbose = verbose)
}

# Feature lines of a GFF/GTF file (comments, blank lines and FASTA removed).
.gff_feature_lines <- function(file) {
  lines <- readLines(file, warn = FALSE)
  fa <- which(startsWith(lines, "##FASTA"))
  if (length(fa)) lines <- lines[seq_len(fa[1] - 1L)]
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) .stop("No feature lines found in GFF file '", file, "'.")
  lines
}

.parse_gff_lines <- function(lines, file, chr = NULL, start = NULL, end = NULL, chr_map = NULL,
                             gene_types = .gene_types_default, verbose = TRUE) {
  tab <- regexpr("\t", lines, fixed = TRUE)
  if (any(tab < 0)) {
    .warn("GFF: ignored ", sum(tab < 0), " line(s) that are not tab-separated.")
    lines <- lines[tab > 0]
    tab <- tab[tab > 0]
  }
  seqid <- substr(lines, 1L, tab - 1L)
  useq <- unique(seqid)
  nseq <- normalize_chr(useq, chr_map)
  chr_all <- nseq[match(seqid, useq)]
  all_chr <- sort_chr(nseq)
  if (!is.null(chr)) {
    tchr <- normalize_chr(chr, chr_map)
    keep <- chr_all == tchr
    if (!any(keep)) {
      .warn("GFF: no features on chromosome '", chr, "'. Chromosomes in '", basename(file), "': ",
            paste(utils::head(all_chr, 20), collapse = ", "),
            if (length(all_chr) > 20) ", ..." else "", ". Check chromosome naming or `chr_map`.")
      return(.empty_genes(file, all_chr))
    }
    lines <- lines[keep]
    chr_all <- chr_all[keep]
    seqid <- seqid[keep]
  }
  f <- data.table::fread(text = lines, sep = "\t", header = FALSE, quote = "", skip = 0,
                         colClasses = "character", fill = TRUE, data.table = FALSE,
                         showProgress = FALSE)
  if (ncol(f) < 9L) .stop("GFF file '", file, "' must have 9 tab-separated columns.")
  a9 <- f[[9]]
  d <- data.frame(row = seq_along(lines), chr = chr_all, seqid = seqid, type = f[[3]],
                  start = suppressWarnings(as.numeric(f[[4]])),
                  end = suppressWarnings(as.numeric(f[[5]])), strand = f[[7]],
                  stringsAsFactors = FALSE)
  is_gtf <- mean(grepl("^[[:space:]]*[A-Za-z_]+ \"", utils::head(a9, 500))) > 0.5
  if (is_gtf) {
    gid <- .attr_gtf(a9, "gene_id")
    tid <- .attr_gtf(a9, "transcript_id")
    d$id <- ifelse(d$type == "gene", gid, ifelse(d$type %in% c("transcript", "mRNA"), tid, NA))
    d$parent <- ifelse(d$type %in% c("transcript", "mRNA"), gid,
                       ifelse(d$type %in% .sub_types, tid, NA))
    d$name <- .attr_gtf(a9, "gene_name")
    d$note <- .attr_gtf(a9, "gene_biotype")
    d$gene_hint <- gid
  } else {
    d$id <- .attr_gff3(a9, "ID")
    d$parent <- .attr_gff3(a9, "Parent")
    d$name <- .attr_gff3(a9, "Name")
    gn <- .attr_gff3(a9, "gene_name")
    d$name[is.na(d$name)] <- gn[is.na(d$name)]
    note <- .attr_gff3(a9, "Note")
    de <- .attr_gff3(a9, "description")
    pr <- .attr_gff3(a9, "product")
    d$note <- ifelse(!is.na(note), note, ifelse(!is.na(de), de, pr))
    d$gene_hint <- NA_character_
  }
  badc <- is.na(d$start) | is.na(d$end)
  if (any(badc)) {
    .warn("GFF: ignored ", sum(badc), " line(s) with non-numeric coordinates.")
    d <- d[!badc, , drop = FALSE]
  }
  mod <- .build_gene_models(d, gene_types)
  if (!is.null(start) || !is.null(end)) {
    s0 <- if (is.null(start)) -Inf else start
    e0 <- if (is.null(end)) Inf else end
    keep_g <- mod$genes$gene_id[mod$genes$start <= e0 & mod$genes$end >= s0]
    mod <- .subset_models(mod, keep_g)
  }
  rows <- unique(rbind(
    data.frame(row = mod$genes$row, gene_id = mod$genes$gene_id, level = "gene", stringsAsFactors = FALSE),
    data.frame(row = mod$transcripts$row, gene_id = mod$transcripts$gene_id, level = "transcript",
               stringsAsFactors = FALSE),
    data.frame(row = mod$features$row, gene_id = mod$features$gene_id, level = "feature",
               stringsAsFactors = FALSE)))
  rows <- rows[!is.na(rows$row), , drop = FALSE]
  rows$line <- lines[rows$row]
  out <- structure(list(genes = mod$genes, transcripts = mod$transcripts, features = mod$features,
                        raw = rows, format = if (is_gtf) "gtf" else "gff3", file = file,
                        seqids = all_chr), class = "qp_genes")
  .msg(verbose, "Genes: ", .fmt(nrow(out$genes)), " gene(s) and ", .fmt(nrow(out$transcripts)),
       " transcript(s) read from ", basename(file),
       if (!is.null(chr)) paste0(" (chromosome ", normalize_chr(chr, chr_map), ")") else "", ".")
  out
}

.empty_genes <- function(file, seqids = character()) {
  structure(list(
    genes = data.frame(gene_id = character(), chr = character(), seqid = character(),
                       start = numeric(), end = numeric(), strand = character(),
                       name = character(), note = character(), type = character(),
                       row = integer(), canonical = character(), stringsAsFactors = FALSE),
    transcripts = data.frame(transcript_id = character(), gene_id = character(),
                             chr = character(), start = numeric(), end = numeric(),
                             strand = character(), type = character(), exon_len = numeric(),
                             cds_len = numeric(), canonical = logical(), row = integer(),
                             stringsAsFactors = FALSE),
    features = data.frame(transcript_id = character(), gene_id = character(),
                          class = character(), start = numeric(), end = numeric(),
                          row = integer(), stringsAsFactors = FALSE),
    raw = data.frame(row = integer(), gene_id = character(), level = character(),
                     line = character(), stringsAsFactors = FALSE),
    format = "gff3", file = file, seqids = seqids), class = "qp_genes")
}

.build_gene_models <- function(d, gene_types) {
  is_gene <- d$type %in% gene_types
  is_sub <- d$type %in% .sub_types
  g <- d[is_gene, , drop = FALSE]
  g$id[is.na(g$id)] <- paste0("gene_line", g$row[is.na(g$id)])
  g <- g[!duplicated(g$id), , drop = FALSE]

  fp <- sub(",.*$", "", d$parent)
  sub_par <- unique(unlist(strsplit(d$parent[is_sub & !is.na(d$parent)], ",", fixed = TRUE)))
  cand <- !is_gene & !is_sub & !(d$type %in% .skip_types)
  is_tx <- cand & ((!is.na(fp) & fp %in% g$id) | (!is.na(d$id) & d$id %in% sub_par) |
                     (is.na(d$parent) & d$type %in% c("mRNA", "transcript")))
  tx <- d[is_tx, , drop = FALSE]
  tx$gene_id <- fp[is_tx]
  tx$id[is.na(tx$id)] <- paste0("transcript_line", tx$row[is.na(tx$id)])
  tx <- tx[!duplicated(tx$id), , drop = FALSE]
  nog <- is.na(tx$gene_id) | !(tx$gene_id %in% g$id)
  if (any(nog)) {
    tx$gene_id[nog] <- ifelse(!is.na(tx$gene_id[nog]), tx$gene_id[nog],
                              ifelse(!is.na(tx$gene_hint[nog]), tx$gene_hint[nog], tx$id[nog]))
  }

  s <- d[is_sub, , drop = FALSE]
  if (nrow(s)) {
    pl <- strsplit(s$parent, ",", fixed = TRUE)
    pl[vapply(pl, length, 1L) == 0L] <- list(NA_character_)
    s <- s[rep(seq_len(nrow(s)), lengths(pl)), , drop = FALSE]
    s$parent <- unlist(pl, use.names = FALSE)
  }
  s$class <- .sub_class(s$type)
  s$transcript_id <- s$parent
  s$gene_id <- rep(NA_character_, nrow(s))
  in_tx <- !is.na(s$parent) & s$parent %in% tx$id
  s$gene_id[in_tx] <- tx$gene_id[match(s$parent[in_tx], tx$id)]
  on_gene <- !in_tx & !is.na(s$parent) & s$parent %in% g$id
  s$transcript_id[on_gene] <- paste0(s$parent[on_gene], ".model")
  s$gene_id[on_gene] <- s$parent[on_gene]
  orphan <- !in_tx & !on_gene
  if (any(orphan)) {
    nop <- orphan & is.na(s$parent)
    s$transcript_id[nop] <- paste0("feature_line", s$row[nop])
    s$gene_id[orphan] <- ifelse(!is.na(s$gene_hint[orphan]), s$gene_hint[orphan],
                                s$transcript_id[orphan])
  }

  # implicit transcripts for sub-features without a transcript feature
  need_tx <- setdiff(unique(s$transcript_id), tx$id)
  if (length(need_tx)) {
    ss <- s[s$transcript_id %in% need_tx, , drop = FALSE]
    a <- .agg_ranges(ss, "transcript_id")
    a$gene_id <- ss$gene_id[match(a$id, ss$transcript_id)]
    add <- data.frame(row = NA_integer_, chr = a$chr, seqid = a$seqid, type = "transcript",
                      start = a$start, end = a$end, strand = a$strand, id = a$id,
                      parent = a$gene_id, name = NA_character_, note = NA_character_,
                      gene_hint = NA_character_, gene_id = a$gene_id, stringsAsFactors = FALSE)
    tx <- rbind(tx, add[, names(tx)])
  }
  # implicit genes for transcripts whose gene has no gene feature
  need_g <- setdiff(unique(tx$gene_id), g$id)
  if (length(need_g)) {
    tt <- tx[tx$gene_id %in% need_g, , drop = FALSE]
    a <- .agg_ranges(tt, "gene_id")
    nm <- tt$name[match(a$id, tt$gene_id)]
    add <- data.frame(row = NA_integer_, chr = a$chr, seqid = a$seqid, type = "gene",
                      start = a$start, end = a$end, strand = a$strand, id = a$id,
                      parent = NA_character_, name = nm, note = NA_character_,
                      gene_hint = NA_character_, stringsAsFactors = FALSE)
    g <- rbind(g, add[, names(g)])
  }
  # genes without any transcript: draw the gene body as one block
  lone <- !(g$id %in% tx$gene_id)
  if (any(lone)) {
    gl <- g[lone, , drop = FALSE]
    add_tx <- data.frame(row = NA_integer_, chr = gl$chr, seqid = gl$seqid, type = "gene_only",
                         start = gl$start, end = gl$end, strand = gl$strand,
                         id = paste0(gl$id, ".gene"), parent = gl$id, name = gl$name,
                         note = NA_character_, gene_hint = NA_character_, gene_id = gl$id,
                         stringsAsFactors = FALSE)
    tx <- rbind(tx, add_tx[, names(tx)])
    add_s <- data.frame(row = NA_integer_, chr = gl$chr, seqid = gl$seqid, type = "exon",
                        start = gl$start, end = gl$end, strand = gl$strand,
                        id = NA_character_, parent = add_tx$id, name = NA_character_,
                        note = NA_character_, gene_hint = NA_character_, class = "exon",
                        transcript_id = add_tx$id, gene_id = gl$id, stringsAsFactors = FALSE)
    s <- rbind(s, add_s[, names(s)])
  }
  # exons from CDS/UTR when a transcript has no exon features
  no_ex <- setdiff(unique(s$transcript_id), unique(s$transcript_id[s$class == "exon"]))
  if (length(no_ex)) {
    parts <- s[s$transcript_id %in% no_ex, , drop = FALSE]
    derived <- lapply(split(parts, parts$transcript_id), function(p) {
      m <- .merge_intervals(p$start, p$end)
      data.frame(row = NA_integer_, chr = p$chr[1], seqid = p$seqid[1], type = "exon",
                 start = m$start, end = m$end, strand = p$strand[1], id = NA_character_,
                 parent = p$transcript_id[1], name = NA_character_, note = NA_character_,
                 gene_hint = NA_character_, class = "exon", transcript_id = p$transcript_id[1],
                 gene_id = p$gene_id[1], stringsAsFactors = FALSE)
    })
    s <- rbind(s, do.call(rbind, derived)[, names(s)])
  }

  w <- s$end - s$start + 1
  ex_len <- tapply(w[s$class == "exon"], s$transcript_id[s$class == "exon"], sum)
  cds_len <- tapply(w[s$class == "CDS"], s$transcript_id[s$class == "CDS"], sum)
  tx$exon_len <- as.numeric(ex_len[tx$id])
  tx$cds_len <- as.numeric(cds_len[tx$id])
  tx$exon_len[is.na(tx$exon_len)] <- 0
  tx$cds_len[is.na(tx$cds_len)] <- 0
  o <- order(tx$gene_id, -tx$cds_len, -tx$exon_len, ifelse(is.na(tx$row), Inf, tx$row))
  first <- !duplicated(tx$gene_id[o])
  canon <- stats::setNames(tx$id[o][first], tx$gene_id[o][first])
  tx$canonical <- tx$id %in% canon

  genes <- data.frame(gene_id = g$id, chr = g$chr, seqid = g$seqid, start = g$start, end = g$end,
                      strand = g$strand, name = ifelse(is.na(g$name), g$id, g$name),
                      note = g$note, type = g$type, row = g$row,
                      canonical = unname(canon[g$id]), stringsAsFactors = FALSE)
  genes <- genes[order(genes$start, genes$end), , drop = FALSE]
  transcripts <- data.frame(transcript_id = tx$id, gene_id = tx$gene_id, chr = tx$chr,
                            start = tx$start, end = tx$end, strand = tx$strand, type = tx$type,
                            exon_len = tx$exon_len, cds_len = tx$cds_len,
                            canonical = tx$canonical, row = tx$row, stringsAsFactors = FALSE)
  features <- data.frame(transcript_id = s$transcript_id, gene_id = s$gene_id, class = s$class,
                         start = s$start, end = s$end, row = s$row, stringsAsFactors = FALSE)
  rownames(genes) <- rownames(transcripts) <- rownames(features) <- NULL
  list(genes = genes, transcripts = transcripts, features = features)
}

.subset_models <- function(mod, gene_ids) {
  mod$genes <- mod$genes[mod$genes$gene_id %in% gene_ids, , drop = FALSE]
  mod$transcripts <- mod$transcripts[mod$transcripts$gene_id %in% gene_ids, , drop = FALSE]
  mod$features <- mod$features[mod$features$gene_id %in% gene_ids, , drop = FALSE]
  if (!is.null(mod$raw)) mod$raw <- mod$raw[mod$raw$gene_id %in% gene_ids, , drop = FALSE]
  mod
}

#' @export
print.qp_genes <- function(x, ...) {
  cat(sprintf("<qp_genes> %s gene(s), %s transcript(s) from %s (%s)\n", .fmt(nrow(x$genes)),
              .fmt(nrow(x$transcripts)), basename(x$file), toupper(x$format)))
  if (nrow(x$genes)) print(utils::head(x$genes[, c("gene_id", "chr", "start", "end", "strand", "name")]))
  invisible(x)
}

#' Genes overlapping a region
#'
#' @param genes Output of [read_gff()].
#' @param region A region (list with `chr`, `start`, `end`, e.g. from
#'   [define_region()]); alternatively give `chr`, `start` and `end`.
#' @param chr,start,end Region coordinates if `region` is not given.
#' @return A `qp_genes` object restricted to genes overlapping the region.
#' @export
genes_in_region <- function(genes, region = NULL, chr = NULL, start = NULL, end = NULL) {
  if (!inherits(genes, "qp_genes")) .stop("`genes` must come from read_gff().")
  r <- if (!is.null(region)) .as_region(region) else .as_region(list(chr = chr, start = start, end = end))
  gg <- genes$genes
  keep <- gg$chr == r$chr & gg$start <= r$end & gg$end >= r$start
  out <- .subset_models(genes, gg$gene_id[keep])
  attr(out, "region") <- r
  out
}

#' Write the GFF records of the genes in a region
#'
#' Writes the original, unmodified annotation lines of the genes overlapping
#' the region to a text file in the same format as the input (GFF3 or GTF),
#' so the candidate genes and their locations can be opened in a spreadsheet,
#' genome browser or any GFF tool.
#'
#' @param genes Output of [read_gff()].
#' @param file Output path (for example `"locus_genes.txt"`).
#' @param region Optional region; if `NULL`, all genes in `genes` are written.
#' @param features `"gene"` writes only the gene lines (one line per gene);
#'   `"all"` also writes the mRNA, exon, CDS and UTR lines of those genes.
#' @param header Add a `##gff-version 3` line and a comment with the region.
#' @return The file path (invisibly).
#' @export
write_region_genes <- function(genes, file, region = NULL, features = c("gene", "all"),
                               header = TRUE) {
  features <- match.arg(features)
  if (!inherits(genes, "qp_genes")) .stop("`genes` must come from read_gff().")
  r <- NULL
  if (!is.null(region)) {
    r <- .as_region(region)
    genes <- genes_in_region(genes, r)
  }
  raw <- genes$raw
  if (features == "gene") {
    with_gene_line <- unique(raw$gene_id[raw$level == "gene"])
    raw <- raw[raw$level == "gene" | (raw$level == "transcript" & !(raw$gene_id %in% with_gene_line)), ,
               drop = FALSE]
  }
  raw <- raw[!duplicated(raw$row), , drop = FALSE]
  raw <- raw[order(raw$row), , drop = FALSE]
  hdr <- character()
  if (header) {
    hdr <- if (genes$format == "gtf") "#gtf" else "##gff-version 3"
    where <- if (!is.null(r)) {
      sprintf(" overlapping chromosome %s:%s-%s", r$chr, .fmt(round(r$start)), .fmt(round(r$end)))
    } else ""
    hdr <- c(hdr, sprintf("# %d gene(s)%s; original records from %s", nrow(genes$genes), where,
                          basename(genes$file)))
  }
  .ensure_dir(file)
  writeLines(c(hdr, raw$line), file)
  invisible(file)
}
