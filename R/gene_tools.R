# Gene lists and single-gene structure plots ---------------------------------

# Regions from "chr:start-end" strings or a data frame (chr, start, end, name).
.parse_regions <- function(regions, chr_map = NULL) {
  if (is.character(regions)) {
    r <- gsub("[, ]", "", regions)
    m <- regmatches(r, regexec("^([^:]+):([0-9.eE+]+)-([0-9.eE+]+)$", r))
    bad <- lengths(m) != 4L
    if (any(bad)) {
      .stop("Regions must look like '9:12400000-12700000'. Not understood: ",
            paste(regions[bad], collapse = ", "), ".")
    }
    d <- data.frame(chr = vapply(m, `[`, "", 2), start = as.numeric(vapply(m, `[`, "", 3)),
                    end = as.numeric(vapply(m, `[`, "", 4)), name = regions, stringsAsFactors = FALSE)
  } else if (is.data.frame(regions)) {
    nm <- names(regions)
    cc <- .find_col(nm, c("chr", "chrom", "chromosome"))
    sc <- .find_col(nm, c("start", "from", "bp1", "begin"))
    ec <- .find_col(nm, c("end", "to", "bp2", "stop"))
    nc <- .find_col(nm, c("name", "region", "label", "id", "qtl"))
    if (anyNA(c(cc, sc, ec))) .stop("`regions` needs columns chr, start and end.")
    d <- data.frame(chr = as.character(regions[[cc]]), start = as.numeric(regions[[sc]]),
                    end = as.numeric(regions[[ec]]),
                    name = if (!is.na(nc)) as.character(regions[[nc]]) else NA_character_,
                    stringsAsFactors = FALSE)
  } else {
    .stop("`regions` must be a data frame (chr, start, end) or strings like '9:12400000-12700000'.")
  }
  d$chr <- normalize_chr(d$chr, chr_map)
  if (anyNA(d$start) || anyNA(d$end) || any(d$end < d$start)) .stop("Each region needs start <= end.")
  d$name[is.na(d$name)] <- sprintf("%s:%s-%s", d$chr[is.na(d$name)], d$start[is.na(d$name)],
                                   d$end[is.na(d$name)])
  d
}

# Original annotation lines of the given genes (gene lines only or all).
.raw_lines_for <- function(genes, gene_ids, features = "gene") {
  raw <- genes$raw[genes$raw$gene_id %in% gene_ids, , drop = FALSE]
  if (features == "gene") {
    with_gene_line <- unique(raw$gene_id[raw$level == "gene"])
    raw <- raw[raw$level == "gene" | (raw$level == "transcript" & !(raw$gene_id %in% with_gene_line)), ,
               drop = FALSE]
  }
  raw <- raw[!duplicated(raw$row), , drop = FALSE]
  raw$line[order(raw$row)]
}

#' List the genes in one or more genomic ranges
#'
#' Returns the genes overlapping a range (or several ranges, or a window
#' around one or more positions) as a table, and optionally writes it as a
#' tab-separated text file or as the original GFF records.
#'
#' @param gff GFF3/GTF file or [read_gff()] output.
#' @param chr,start,end One range (bp). Any chromosome naming style works.
#' @param pos,flank Alternative to `start`/`end`: position(s) and the number
#'   of bp to add on each side (for example `pos = 12550313, flank = 100000`).
#'   A `distance_bp` column then gives the distance of each gene to the
#'   nearest position (0 = overlapping).
#' @param regions Several ranges: a data frame with `chr`, `start`, `end`
#'   (and optional `name`), or strings such as `"9:12400000-12700000"`.
#' @param output Optional output file. For `format = "table"` the file type
#'   follows the name: `.csv` (opens in Excel), `.xlsx` (Excel workbook; needs
#'   the `writexl` package) or `.txt` (tab-separated).
#' @param format `"table"` (a table of genes) or `"gff"` (the original GFF
#'   lines, same format as the input file).
#' @param features For `format = "gff"`: `"gene"` lines only or `"all"`
#'   (gene, mRNA, exon, CDS, UTR lines).
#' @param chr_map Optional chromosome name mapping (see [normalize_chr()]).
#' @param verbose Print messages.
#' @return A data frame: `region`, `gene_id`, `name`, `chr`, `start`, `end`,
#'   `strand`, `length_bp`, `n_transcripts`, `description` (and
#'   `distance_bp` when `pos` is given).
#' @examples
#' \donttest{
#' ex <- simulate_example_data(tempfile("qp"))
#' get_genes(ex$gff, chr = 9, start = 12500000, end = 12600000)
#' get_genes(ex$gff, regions = c("9:12450000-12470000", "Chr9:12600000-12650000"))
#' }
#' @export
get_genes <- function(gff, chr = NULL, start = NULL, end = NULL, pos = NULL, flank = 0,
                      regions = NULL, output = NULL, format = c("table", "gff"),
                      features = c("gene", "all"), chr_map = NULL, verbose = TRUE) {
  format <- match.arg(format)
  features <- match.arg(features)
  if (is.null(regions)) {
    if (is.null(chr)) .stop("Give `chr` with `start` and `end` (or `pos` and `flank`), or a `regions` table.")
    if (!is.null(pos)) {
      s <- min(pos) - flank
      e <- max(pos) + flank
    } else {
      if (is.null(start) || is.null(end)) .stop("Give both `start` and `end`, or `pos` (and `flank`).")
      s <- start
      e <- end
    }
    regions <- data.frame(chr = as.character(chr), start = max(1, s), end = e, name = NA_character_)
  }
  reg <- .parse_regions(regions, chr_map)
  genes <- if (inherits(gff, "qp_genes")) {
    gff
  } else if (length(unique(reg$chr)) == 1L) {
    read_gff(gff, chr = reg$chr[1], chr_map = chr_map, verbose = verbose)
  } else {
    read_gff(gff, chr_map = chr_map, verbose = verbose)
  }
  out <- data.frame(region = character(), gene_id = character(), name = character(), chr = character(),
                    start = numeric(), end = numeric(), strand = character(), length_bp = numeric(),
                    n_transcripts = integer(), description = character(), stringsAsFactors = FALSE)
  for (k in seq_len(nrow(reg))) {
    sub <- genes_in_region(genes, chr = reg$chr[k], start = reg$start[k], end = reg$end[k])
    g <- sub$genes
    if (!nrow(g)) next
    ntx <- table(sub$transcripts$gene_id[sub$transcripts$type != "gene_only"])
    out <- rbind(out, data.frame(region = reg$name[k], gene_id = g$gene_id, name = g$name, chr = g$chr,
                                 start = g$start, end = g$end, strand = g$strand,
                                 length_bp = g$end - g$start + 1,
                                 n_transcripts = as.integer(ifelse(is.na(ntx[g$gene_id]), 0L, ntx[g$gene_id])),
                                 description = ifelse(is.na(g$note), "", g$note), stringsAsFactors = FALSE))
  }
  if (!is.null(pos) && nrow(out)) {
    out$distance_bp <- vapply(seq_len(nrow(out)), function(i) {
      d <- ifelse(pos >= out$start[i] & pos <= out$end[i], 0,
                  pmin(abs(out$start[i] - pos), abs(out$end[i] - pos)))
      min(d)
    }, 0)
  }
  rownames(out) <- NULL
  .msg(verbose, .fmt(nrow(out)), " gene(s) in ", nrow(reg), " range(s).")
  if (!is.null(output)) {
    .ensure_dir(output)
    if (format == "table") {
      output <- write_table(out, output, sheet = "Genes")
    } else {
      hdr <- c(if (genes$format == "gtf") "#gtf" else "##gff-version 3",
               paste0("# ", length(unique(out$gene_id)), " gene(s) in: ", paste(reg$name, collapse = "; ")))
      writeLines(c(hdr, .raw_lines_for(genes, unique(out$gene_id), features)), output)
    }
    .msg(verbose, "Gene list written to ", paste(output, collapse = ", "), ".")
  }
  out
}

# Find a gene (ID, name or transcript ID) and the annotation of its chromosome.
.locate_gene <- function(gff, gene, chr_map = NULL) {
  if (inherits(gff, "qp_genes")) {
    genes <- gff
  } else {
    .check_file(gff, "GFF")
    lines <- .gff_feature_lines(gff)
    hit <- grep(gene, lines, fixed = TRUE)
    if (!length(hit)) hit <- grep(gene, lines, ignore.case = TRUE)
    if (!length(hit)) .stop("Gene '", gene, "' was not found in '", basename(gff), "'.")
    tab_all <- regexpr("\t", lines, fixed = TRUE)
    seq_all <- substr(lines, 1L, tab_all - 1L)
    chrs <- unique(normalize_chr(seq_all[hit], chr_map))
    useq <- unique(seq_all)
    keep_seq <- useq[normalize_chr(useq, chr_map) %in% chrs]
    genes <- .parse_gff_lines(lines[seq_all %in% keep_seq], gff, chr_map = chr_map, verbose = FALSE)
  }
  gg <- genes$genes
  tx <- genes$transcripts
  tid <- NULL
  k <- which(gg$gene_id == gene | gg$name == gene)
  if (!length(k)) {
    t <- which(tx$transcript_id == gene)
    if (length(t)) {
      k <- which(gg$gene_id == tx$gene_id[t[1]])
      tid <- gene
    }
  }
  if (!length(k)) k <- which(toupper(gg$gene_id) == toupper(gene) | toupper(gg$name) == toupper(gene))
  if (!length(k)) {
    cand <- unique(c(gg$gene_id[grepl(gene, gg$gene_id, ignore.case = TRUE)],
                     gg$name[grepl(gene, gg$name, ignore.case = TRUE)]))
    .stop("Gene '", gene, "' was not found",
          if (length(cand)) paste0("; similar: ", paste(utils::head(cand, 8), collapse = ", ")) else "", ".")
  }
  if (length(k) > 1L) .warn("Several genes match '", gene, "'; using ", gg$gene_id[k[1]], ".")
  list(genes = genes, gene = gg[k[1], , drop = FALSE], transcript = tid)
}

# Location of variants relative to a transcript model.
.variant_location <- function(pos, gene, feats, strand) {
  out <- rep("intron", length(pos))
  ex <- feats[feats$class == "exon", , drop = FALSE]
  cds <- feats[feats$class == "CDS", , drop = FALSE]
  in_any <- function(p, s, e) vapply(p, function(z) any(z >= s & z <= e), NA)
  if (nrow(ex)) {
    in_ex <- in_any(pos, ex$start, ex$end)
    if (nrow(cds)) {
      in_cds <- in_any(pos, cds$start, cds$end)
      lo <- min(cds$start)
      hi <- max(cds$end)
      utr_left <- if (strand == "-") "3'UTR" else "5'UTR"
      utr_right <- if (strand == "-") "5'UTR" else "3'UTR"
      out[in_ex & !in_cds & pos < lo] <- utr_left
      out[in_ex & !in_cds & pos > hi] <- utr_right
      out[in_ex & !in_cds & pos >= lo & pos <= hi] <- "exon (non-coding part)"
      out[in_cds] <- "CDS"
    } else {
      out[in_ex] <- "exon (non-coding transcript)"
    }
  }
  up <- if (strand == "-") pos > gene$end else pos < gene$start
  down <- if (strand == "-") pos < gene$start else pos > gene$end
  out[up] <- "upstream"
  out[down] <- "downstream"
  out
}

#' Plot the structure of a single gene
#'
#' Two styles are available. `"lollipop"` (default, similar to geneHapR)
#' draws the representative transcript on a single gene axis: coding exons
#' (CDS), 5' UTRs and 3' UTRs as boxes in different colours on a line that
#' stands for the introns and flanking regions, an arrow at the transcription
#' start site, and each variant as a "balloon" on a stem whose top carries
#' its alleles (`A/T`; indels as `-/TTAAA`). Stems fan out so that the labels
#' do not overlap; with GWAS results the balloons are coloured by
#' -log10 P. `"classic"` draws all (or the representative) transcripts with
#' exon numbers and a separate lollipop track above the gene, with
#' heights = -log10 P. In both styles each variant is classified as
#' upstream, 5'UTR, CDS, intron, 3'UTR or downstream (returned in the
#' `variants` attribute).
#'
#' @param gff GFF3/GTF file or [read_gff()] output.
#' @param gene Gene ID, gene name or transcript ID (e.g. `"LOC_Os09g12345"`).
#' @param transcripts `"all"` transcripts or only the `"canonical"` one.
#' @param flank Bases shown on each side of the gene (and searched for variants).
#' @param variants Optional variants: a VCF or HapMap file, a `qp_geno`
#'   object, or a data frame with `pos` (and optional `type`, `id`).
#' @param gwas Optional GWAS results ([read_gwas()] or file); lollipop
#'   heights become -log10 P.
#' @param style `"lollipop"` (single gene axis with allele balloons) or
#'   `"classic"` (all transcripts, P-value lollipop track).
#' @param variant_labels Text at the top of each balloon (`"lollipop"`
#'   style): `"alleles"` (from VCF/HapMap; the variant ID when unknown),
#'   `"id"`, `"position"` or `"none"`.
#' @param exon_numbers Number the exons (`"classic"` style).
#' @param intron_style `"hat"` (angled introns) or `"line"`.
#' @param coordinates `"genomic"` (Mb) or `"relative"` (kb from the gene start).
#' @param colors Colours of the boxes: `CDS`, `UTR5` and `UTR3` (lollipop
#'   style), `UTR` (classic style) and non-coding `exon`.
#' @param variant_colors Colours for `SNP` and `indel`.
#' @param label_variants Label variants with their IDs (useful for few variants).
#' @param chr_map Optional chromosome name mapping.
#' @param base_size Base font size (pt).
#' @param title Add a title with gene ID, location, strand and description.
#' @return A ggplot (or patchwork) object; attribute `variants` holds the
#'   variant table with locations. Save with [save_pdf()].
#' @examples
#' \donttest{
#' ex <- simulate_example_data(tempfile("qp"))
#' p <- plot_gene_structure(ex$gff, "Sim09g12520", variants = ex$hapmap, gwas = ex$gwas)
#' save_pdf(p, tempfile(fileext = ".pdf"), width = 170, height = 90)
#' attr(p, "variants")
#' }
#' @export
plot_gene_structure <- function(gff, gene, transcripts = c("all", "canonical"), flank = 0,
                                variants = NULL, gwas = NULL, style = c("lollipop", "classic"),
                                variant_labels = c("alleles", "id", "position", "none"), exon_numbers = TRUE,
                                intron_style = c("hat", "line"), coordinates = c("genomic", "relative"),
                                colors = c(CDS = "#2B5C8A", UTR5 = "#E69F00", UTR3 = "#56B4E9",
                                           UTR = "#9ECAE1", exon = "#74A9CF"),
                                variant_colors = c(SNP = "#D55E00", indel = "#CC79A7"),
                                label_variants = FALSE, chr_map = NULL, base_size = 9, title = TRUE) {
  style <- match.arg(style)
  variant_labels <- match.arg(variant_labels)
  transcripts <- match.arg(transcripts)
  intron_style <- match.arg(intron_style)
  coordinates <- match.arg(coordinates)
  loc <- .locate_gene(gff, gene, chr_map)
  g <- loc$gene
  G <- loc$genes
  tx <- G$transcripts[G$transcripts$gene_id == g$gene_id, , drop = FALSE]
  if (transcripts == "canonical") tx <- tx[tx$canonical, , drop = FALSE]
  tx <- tx[order(!tx$canonical, tx$transcript_id), , drop = FALSE]
  f <- G$features[G$features$transcript_id %in% tx$transcript_id, , drop = FALSE]
  r0 <- max(1, g$start - flank)
  r1 <- g$end + flank
  if (coordinates == "relative") {
    tr <- function(x) (x - g$start) / 1e3
    xlab <- "Position relative to the gene start (kb)"
  } else {
    tr <- function(x) x / 1e6
    xlab <- sprintf("Chromosome %s (Mb)", g$chr)
  }
  n <- nrow(tx)
  tx$y <- -seq_len(n)
  f$y <- tx$y[match(f$transcript_id, tx$transcript_id)]
  coding <- unique(f$transcript_id[f$class == "CDS"])
  ex <- f[f$class == "exon", , drop = FALSE]
  ex$h <- ifelse(ex$transcript_id %in% coding, 0.16, 0.24)
  ex$fill <- ifelse(ex$transcript_id %in% coding, "UTR", "exon")
  cds <- f[f$class == "CDS", , drop = FALSE]
  w <- r1 - r0

  # introns
  intr <- do.call(rbind, lapply(split(ex, ex$transcript_id), function(e) {
    e <- e[order(e$start), , drop = FALSE]
    if (nrow(e) < 2L) return(NULL)
    a <- e$end[-nrow(e)]
    b <- e$start[-1]
    keep <- b > a + 1
    if (!any(keep)) return(NULL)
    a <- a[keep]
    b <- b[keep]
    y <- e$y[1]
    if (intron_style == "hat") {
      mid <- (a + b) / 2
      data.frame(x = c(a, mid), xend = c(mid, b), y = c(rep(y, length(a)), rep(y + 0.2, length(a))),
                 yend = c(rep(y + 0.2, length(a)), rep(y, length(a))))
    } else {
      data.frame(x = a, xend = b, y = y, yend = y)
    }
  }))
  # exon numbers in transcription order
  exn <- do.call(rbind, lapply(split(ex, ex$transcript_id), function(e) {
    s <- tx$strand[match(e$transcript_id[1], tx$transcript_id)]
    e <- e[order(e$start, decreasing = identical(s, "-")), , drop = FALSE]
    data.frame(x = (e$start + e$end) / 2, y = e$y + 0.34, lab = seq_len(nrow(e)))
  }))
  # TSS arrows
  arr <- tx[tx$strand %in% c("+", "-"), , drop = FALSE]
  if (nrow(arr)) {
    arr$tss <- ifelse(arr$strand == "+", arr$start, arr$end)
    arr$dir <- ifelse(arr$strand == "+", 1, -1)
  }
  ylim_top <- -0.4 + if (exon_numbers) 0.25 else 0
  p <- ggplot2::ggplot()
  if (!is.null(intr) && nrow(intr)) {
    p <- p + ggplot2::geom_segment(data = intr, ggplot2::aes(x = tr(.data$x), xend = tr(.data$xend), y = .data$y,
                                                             yend = .data$yend),
                                   colour = "grey35", linewidth = 0.35)
  }
  p <- p + ggplot2::geom_rect(data = ex, ggplot2::aes(xmin = tr(.data$start), xmax = tr(.data$end + 1),
                                                      ymin = .data$y - .data$h, ymax = .data$y + .data$h,
                                                      fill = .data$fill), colour = "grey20", linewidth = 0.15)
  if (nrow(cds)) {
    cds$fill <- "CDS"
    p <- p + ggplot2::geom_rect(data = cds, ggplot2::aes(xmin = tr(.data$start), xmax = tr(.data$end + 1),
                                                         ymin = .data$y - 0.3, ymax = .data$y + 0.3,
                                                         fill = .data$fill), colour = "grey15", linewidth = 0.15)
  }
  if (nrow(arr)) {
    p <- p +
      ggplot2::geom_segment(data = arr, ggplot2::aes(x = tr(.data$tss), xend = tr(.data$tss), y = .data$y + 0.3,
                                                     yend = .data$y + 0.58), linewidth = 0.4) +
      ggplot2::geom_segment(data = arr, ggplot2::aes(x = tr(.data$tss), xend = tr(.data$tss + .data$dir * 0.045 * w),
                                                     y = .data$y + 0.58, yend = .data$y + 0.58), linewidth = 0.4,
                            arrow = grid::arrow(length = grid::unit(1.3, "mm"), type = "closed", angle = 30))
    ylim_top <- max(ylim_top, -1 + 0.75)
  }
  if (exon_numbers && !is.null(exn) && nrow(exn)) {
    p <- p + ggplot2::geom_text(data = exn, ggplot2::aes(x = tr(.data$x), y = .data$y, label = .data$lab),
                                size = .txt(base_size, 0.65), colour = "grey30", vjust = 0)
  }

  vt <- NULL
  if (!is.null(variants)) {
    if (is.data.frame(variants) && !inherits(variants, "qp_geno")) {
      if (is.null(variants$pos)) .stop("A `variants` data frame needs a `pos` column.")
      vt <- data.frame(id = if (!is.null(variants$id)) as.character(variants$id) else paste0(g$chr, ":", variants$pos),
                       pos = as.numeric(variants$pos),
                       type = if (!is.null(variants$type)) as.character(variants$type) else "SNP",
                       ref = if (!is.null(variants$ref)) as.character(variants$ref) else NA_character_,
                       alt = if (!is.null(variants$alt)) as.character(variants$alt) else NA_character_,
                       stringsAsFactors = FALSE)
      if (!is.null(variants$alleles) && is.null(variants$ref)) {
        al <- strsplit(as.character(variants$alleles), "/", fixed = TRUE)
        vt$ref <- vapply(al, function(z) z[1], "")
        vt$alt <- vapply(al, function(z) paste(z[-1], collapse = ","), "")
      }
    } else {
      gv <- if (inherits(variants, "qp_geno")) variants else
        read_genotypes(variants, chr = g$chr, start = r0, end = r1, chr_map = chr_map, verbose = FALSE)
      inf <- gv$info
      inf <- inf[inf$chr == g$chr & inf$pos >= r0 & inf$pos <= r1, , drop = FALSE]
      vt <- data.frame(id = inf$id, pos = inf$pos, type = ifelse(inf$is_snp, "SNP", "indel"),
                       ref = as.character(inf$ref), alt = as.character(inf$alt), stringsAsFactors = FALSE)
    }
    vt <- vt[vt$pos >= r0 & vt$pos <= r1, , drop = FALSE]
  }
  if (!is.null(gwas)) {
    gw <- .as_gwas(gwas)
    gs <- gw[gw$chr == g$chr & gw$pos >= r0 & gw$pos <= r1, , drop = FALSE]
    if (is.null(vt)) {
      vt <- data.frame(id = gs$snp, pos = gs$pos, type = "SNP", ref = NA_character_, alt = NA_character_,
                       stringsAsFactors = FALSE)
    }
    vt$logp <- gs$logp[match(vt$pos, gs$pos)]
  }
  canon <- tx$transcript_id[tx$canonical][1]
  if (is.na(canon)) canon <- tx$transcript_id[1]
  if (style == "lollipop") {
    return(.gene_lollipop(g, f[f$transcript_id == canon, , drop = FALSE], canon, vt, r0, r1, tr, xlab,
                          colors, variant_colors, variant_labels, base_size, title))
  }
  if (!is.null(vt) && nrow(vt)) {
    vt$location <- .variant_location(vt$pos, g, f[f$transcript_id == canon, , drop = FALSE], g$strand)
    vt$type[!vt$type %in% names(variant_colors)] <- "indel"
    p <- p + ggplot2::geom_segment(data = vt, ggplot2::aes(x = tr(.data$pos), xend = tr(.data$pos), y = -n - 0.45,
                                                           yend = -0.55, colour = .data$type),
                                   linewidth = 0.25, alpha = 0.45, linetype = "22")
  }
  brk <- intersect(c("CDS", "UTR", "exon"), c(ex$fill, if (nrow(cds)) "CDS"))
  lab_map <- c(CDS = "CDS", UTR = "UTR", exon = "exon (non-coding)")
  p <- p +
    ggplot2::scale_fill_manual(values = colors, breaks = brk, labels = unname(lab_map[brk]), name = NULL) +
    ggplot2::scale_colour_manual(values = variant_colors, guide = "none") +
    ggplot2::scale_y_continuous(breaks = tx$y, labels = tx$transcript_id) +
    ggplot2::coord_cartesian(xlim = tr(c(r0, r1)), ylim = c(-n - 0.5, ylim_top + 0.15), expand = FALSE, clip = "off") +
    ggplot2::labs(x = xlab, y = NULL) +
    theme_quard(base_size) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "bottom")
  ttl <- NULL
  if (title) {
    ttl <- sprintf("%s%s  |  chr %s:%s-%s (%s strand)  |  %s bp", g$gene_id,
                   if (!is.na(g$name) && g$name != g$gene_id) paste0(" (", g$name, ")") else "",
                   g$chr, .fmt(g$start), .fmt(g$end), g$strand, .fmt(g$end - g$start + 1))
  }
  sub <- if (title && !is.na(g$note) && nzchar(g$note)) g$note else NULL
  if (!is.null(vt) && nrow(vt)) {
    has_p <- !is.null(vt$logp) && any(!is.na(vt$logp))
    vt$h <- if (has_p) vt$logp else 1
    vt$h[is.na(vt$h)] <- 0
    pv <- ggplot2::ggplot(vt, ggplot2::aes(x = tr(.data$pos))) +
      ggplot2::geom_segment(ggplot2::aes(xend = tr(.data$pos), y = 0, yend = .data$h, colour = .data$type),
                            linewidth = 0.35) +
      ggplot2::geom_point(ggplot2::aes(y = .data$h, fill = .data$type), shape = 21, colour = "grey20",
                          stroke = 0.2, size = 2) +
      ggplot2::scale_colour_manual(values = variant_colors, name = NULL) +
      ggplot2::scale_fill_manual(values = variant_colors, name = NULL) +
      ggplot2::coord_cartesian(xlim = tr(c(r0, r1)), ylim = c(0, max(vt$h, 1) * 1.15), expand = FALSE) +
      ggplot2::labs(x = NULL, y = if (has_p) expression(bold(-log[10](italic(P)))) else "Variants") +
      theme_quard(base_size) + .no_x() +
      ggplot2::theme(legend.position = "right")
    if (!has_p) pv <- pv + .no_y()
    if (label_variants) {
      pv <- pv + ggplot2::geom_text(ggplot2::aes(y = .data$h, label = .data$id), angle = 90, hjust = -0.25,
                                    size = .txt(base_size, 0.6))
    }
    out <- patchwork::wrap_plots(pv, p, ncol = 1, heights = c(1.1, max(1, 0.55 * n + 0.6)))
  } else {
    out <- p
  }
  if (title) {
    out <- out + patchwork::plot_annotation(title = ttl, subtitle = sub,
                                            theme = ggplot2::theme(plot.title = ggplot2::element_text(size = base_size + 1, face = "bold"),
                                                                   plot.subtitle = ggplot2::element_text(size = base_size - 0.5, colour = "grey30")))
  }
  if (!is.null(vt)) vt$h <- NULL
  attr(out, "variants") <- vt
  attr(out, "gene") <- g
  out
}


# Overlaps of points (variants) with intervals, joined on one key column.
# Returns the matching (point index, interval index) pairs.
.point_overlaps <- function(key_p, pos, key_i, start, end) {
  if (!length(pos) || !length(start)) return(list(p = integer(), i = integer()))
  xs <- data.table::data.table(.p = seq_along(pos), .key = as.character(key_p),
                               .s = as.numeric(pos), .e = as.numeric(pos))
  ys <- data.table::data.table(.i = seq_along(start), .key = as.character(key_i),
                               .lo = as.numeric(start), .hi = as.numeric(end))
  data.table::setkeyv(ys, c(".key", ".lo", ".hi"))
  ov <- data.table::foverlaps(xs, ys, by.x = c(".key", ".s", ".e"), by.y = c(".key", ".lo", ".hi"),
                              type = "any", nomatch = 0L)
  list(p = ov[[".p"]], i = ov[[".i"]])
}

# Exon / intron status of variants inside genes. `gid` and `strand` give the
# gene each position is tested against. Exon hits are refined to CDS, 5'UTR,
# 3'UTR (GFF UTR records, or the position relative to the gene's coding span)
# or non-coding exon. Genes without exon/CDS records give "genic".
.genic_status <- function(p, gid, strand, genes) {
  m <- length(p)
  loc <- rep("intron", m)
  feat <- rep(NA_character_, m)
  if (!m) return(list(location = loc, feature = feat))
  Fe <- genes$features[genes$features$gene_id %in% unique(gid), , drop = FALSE]
  if (nrow(Fe)) {
    ov <- .point_overlaps(gid, p, Fe$gene_id, Fe$start, Fe$end)
    hit <- function(cl) {
      z <- logical(m)
      z[unique(ov$p[Fe$class[ov$i] %in% cl])] <- TRUE
      z
    }
    in_ex <- hit(c("exon", "UTR", "UTR5", "UTR3", "CDS"))
    in_cds <- hit("CDS")
    in_u5 <- hit("UTR5")
    in_u3 <- hit("UTR3")
    in_utr <- hit("UTR")
    loc[in_ex] <- "exon"
    feat[in_ex] <- "non-coding exon"
    cds <- Fe[Fe$class == "CDS", , drop = FALSE]
    if (nrow(cds)) {
      clo <- as.numeric(tapply(cds$start, cds$gene_id, min)[gid])
      chi <- as.numeric(tapply(cds$end, cds$gene_id, max)[gid])
      minus <- strand %in% "-"
      left <- in_ex & !is.na(clo) & p < clo
      right <- in_ex & !is.na(chi) & p > chi
      feat[left] <- ifelse(minus[left], "3'UTR", "5'UTR")
      feat[right] <- ifelse(minus[right], "5'UTR", "3'UTR")
      feat[in_utr & !left & !right] <- "UTR"
    } else {
      feat[in_utr] <- "UTR"
    }
    feat[in_u3] <- "3'UTR"
    feat[in_u5] <- "5'UTR"
    feat[in_cds] <- "CDS"
  }
  tx <- genes$transcripts
  if (!is.null(tx) && nrow(tx)) {
    only <- tapply(tx$type == "gene_only", tx$gene_id, all)
    go <- gid %in% names(only)[only %in% TRUE]
    loc[go] <- "genic"
    feat[go] <- NA_character_
  }
  list(location = loc, feature = feat)
}

#' Locate variants relative to genes
#'
#' For each variant, finds the closest gene and classifies the position as
#' `exon` (with `feature` = CDS, 5'UTR, 3'UTR, UTR or non-coding exon),
#' `intron`, `upstream` (promoter region: within `upstream_bp` before the
#' transcription start site, strand-aware), `downstream` (within
#' `downstream_bp` after the gene end), or `intergenic` (then the closest
#' gene is reported with its distance). `genic` is used for genes that have
#' no exon/CDS records in the GFF. Exon/intron status uses all transcripts
#' of a gene (exon in any transcript = exon).
#'
#' When several genes qualify (overlapping genes, bidirectional promoters),
#' the reported gene is chosen by exon > intron > genic > upstream >
#' downstream, then by the shortest distance; the others are listed in
#' `other_genes` as `gene(location)`.
#'
#' `distance_bp` is 0 inside a gene, the distance to the transcription start
#' site for upstream variants, to the transcription end for downstream
#' variants, and to the nearest gene edge for intergenic variants.
#'
#' @param x Variants: a data frame with `pos` (and `chr`, optional `snp`; e.g.
#'   [read_gwas()] output, whose other columns are kept), a `qp_geno` object,
#'   an LD matrix (IDs as dimnames, `pos` attribute; e.g. `res$ld_all` from
#'   [quard_plot()]) or a numeric vector of positions (then give `chr`).
#' @param genes GFF3/GTF file or [read_gff()] output.
#' @param chr Chromosome, if `x` does not contain it.
#' @param upstream_bp Promoter / upstream window (default 3 kb).
#' @param downstream_bp Downstream window (default 1 kb).
#' @param chr_map Optional chromosome name mapping.
#' @return Data frame: `snp`, `chr`, `pos`, `gene_id`, `gene_name`,
#'   `strand`, `location`, `feature`, `distance_bp`, `description`,
#'   `other_genes` (then any other columns of a data-frame `x`).
#' @examples
#' \donttest{
#' ex <- simulate_example_data(tempfile("qp"))
#' annotate_variants(data.frame(chr = 9, pos = c(12480000, 12535820, 12600000)), ex$gff)
#' }
#' @export
annotate_variants <- function(x, genes, chr = NULL, upstream_bp = 3000, downstream_bp = 1000,
                              chr_map = NULL) {
  ann_cols <- c("snp", "chr", "pos", "gene_id", "gene_name", "strand", "location", "feature",
                "distance_bp", "description", "other_genes")
  extra <- NULL
  if (inherits(x, "qp_geno")) {
    d <- data.frame(snp = x$info$id, chr = x$info$chr, pos = x$info$pos, stringsAsFactors = FALSE)
  } else if (is.matrix(x) && !is.null(attr(x, "pos"))) {
    d <- data.frame(snp = rownames(x) %||% paste0("v", seq_len(nrow(x))),
                    chr = attr(x, "chr") %||% chr %||% NA_character_, pos = attr(x, "pos"),
                    stringsAsFactors = FALSE)
  } else if (is.data.frame(x)) {
    if (is.null(x$pos)) .stop("`x` needs a `pos` column.")
    d <- data.frame(snp = if (!is.null(x$snp)) as.character(x$snp) else paste0("v", seq_len(nrow(x))),
                    chr = if (!is.null(x$chr)) as.character(x$chr) else (chr %||% NA_character_),
                    pos = as.numeric(x$pos), stringsAsFactors = FALSE)
    keep <- setdiff(names(x), ann_cols)
    if (length(keep)) extra <- as.data.frame(x, stringsAsFactors = FALSE)[, keep, drop = FALSE]
  } else if (is.numeric(x)) {
    if (is.null(chr)) .stop("Give `chr` for a vector of positions.")
    d <- data.frame(snp = paste0(chr, ":", x), chr = chr, pos = x, stringsAsFactors = FALSE)
  } else {
    .stop("Unsupported `x`.")
  }
  d$chr <- normalize_chr(d$chr, chr_map)
  if (!inherits(genes, "qp_genes")) {
    chrs <- unique(stats::na.omit(d$chr))
    genes <- if (length(chrs) == 1L) read_gff(genes, chr = chrs, chr_map = chr_map, verbose = FALSE) else
      read_gff(genes, chr_map = chr_map, verbose = FALSE)
  }
  n <- nrow(d)
  out <- data.frame(snp = d$snp, chr = d$chr, pos = d$pos, gene_id = rep(NA_character_, n),
                    gene_name = NA_character_, strand = NA_character_, location = "intergenic",
                    feature = NA_character_, distance_bp = NA_real_, description = NA_character_,
                    other_genes = NA_character_, stringsAsFactors = FALSE)
  G <- genes$genes
  ok <- which(!is.na(d$pos) & !is.na(d$chr) & d$chr %in% G$chr)
  if (length(ok)) {
    win <- max(upstream_bp, downstream_bp, 0)
    # 1. candidate (variant, gene) pairs: variant within `win` bp of the gene
    ov <- .point_overlaps(d$chr[ok], d$pos[ok], G$chr, G$start - win, G$end + win)
    i <- ok[ov$p]
    k <- ov$i
    p <- d$pos[i]
    gs <- G$start[k]
    ge <- G$end[k]
    plus <- !(G$strand[k] %in% "-")
    inside <- p >= gs & p <= ge
    up <- !inside & ((plus & p < gs & gs - p <= upstream_bp) | (!plus & p > ge & p - ge <= upstream_bp))
    dn <- !inside & !up & ((plus & p > ge & p - ge <= downstream_bp) |
                             (!plus & p < gs & gs - p <= downstream_bp))
    loc <- rep(NA_character_, length(i))
    feat <- rep(NA_character_, length(i))
    dist <- rep(NA_real_, length(i))
    loc[up] <- "upstream"
    loc[dn] <- "downstream"
    dist[inside] <- 0
    dist[up] <- ifelse(plus[up], gs[up] - p[up], p[up] - ge[up])
    dist[dn] <- ifelse(plus[dn], p[dn] - ge[dn], gs[dn] - p[dn])
    if (any(inside)) {
      st <- .genic_status(p[inside], G$gene_id[k[inside]], G$strand[k[inside]], genes)
      loc[inside] <- st$location
      feat[inside] <- st$feature
    }
    q <- which(!is.na(loc))
    if (length(q)) {
      rk <- match(loc[q], c("exon", "intron", "genic", "upstream", "downstream"))
      q <- q[order(i[q], rk, dist[q], k[q])]
      best <- q[!duplicated(i[q])]
      v <- i[best]
      kb <- k[best]
      out$gene_id[v] <- G$gene_id[kb]
      out$gene_name[v] <- G$name[kb]
      out$strand[v] <- G$strand[kb]
      out$location[v] <- loc[best]
      out$feature[v] <- feat[best]
      out$distance_bp[v] <- dist[best]
      out$description[v] <- G$note[kb]
      oth <- setdiff(q, best)
      if (length(oth)) {
        lab <- tapply(paste0(G$gene_id[k[oth]], "(", loc[oth], ")"), i[oth], paste, collapse = ";")
        out$other_genes[as.integer(names(lab))] <- unname(lab)
      }
    }
    # 2. intergenic variants: closest gene edge on either side
    rest <- ok[is.na(out$gene_id[ok])]
    for (cc in unique(d$chr[rest])) {
      r <- rest[d$chr[rest] == cc]
      gi <- which(G$chr == cc)
      oe <- gi[order(G$end[gi])]
      os <- gi[order(G$start[gi])]
      pr <- d$pos[r]
      jl <- findInterval(pr, G$end[oe], left.open = TRUE)
      jr <- findInterval(pr, G$start[os]) + 1L
      kl <- oe[pmax(jl, 1L)]
      kr <- os[pmin(jr, length(os))]
      dl <- ifelse(jl >= 1L, pr - G$end[kl], Inf)
      dr <- ifelse(jr <= length(os), G$start[kr] - pr, Inf)
      kk <- ifelse(dl <= dr, kl, kr)
      out$gene_id[r] <- G$gene_id[kk]
      out$gene_name[r] <- G$name[kk]
      out$strand[r] <- G$strand[kk]
      out$distance_bp[r] <- pmin(dl, dr)
      out$description[r] <- G$note[kk]
    }
  }
  if (!is.null(extra)) out <- cbind(out, extra)
  rownames(out) <- NULL
  out
}

# Allele label of a variant: "A/T" for SNPs; for indels the shared padding
# base of VCF records is removed and an empty allele is written "-"
# (REF T, ALT TTAAAA -> "-/TAAAA"). Multi-allelic sites list all alleles.
.allele_label <- function(ref, alt) {
  out <- rep(NA_character_, length(ref))
  for (i in seq_along(ref)) {
    r <- ref[i]
    a <- alt[i]
    if (is.na(r) || is.na(a) || !nzchar(r) || !nzchar(a) || a == ".") next
    al <- c(r, strsplit(a, ",", fixed = TRUE)[[1]])
    if (any(nchar(al) > 1L) && all(grepl("^[ACGTNacgtn]+$", al)) &&
        length(unique(toupper(substr(al, 1L, 1L)))) == 1L) {
      al <- substring(al, 2L)
    }
    al[!nzchar(al)] <- "-"
    out[i] <- paste(al, collapse = "/")
  }
  out
}

# Spread sorted label positions so that neighbours are at least `d` apart
# inside [lo, hi]; clusters stay centred on their variants.
.spread_positions <- function(x, d, lo, hi) {
  n <- length(x)
  if (n < 2L) return(pmin(pmax(x, lo), hi))
  o <- order(x)
  y <- x[o]
  if ((n - 1) * d > hi - lo) d <- (hi - lo) / (n - 1)
  pass <- function(v) {
    for (i in 2:n) v[i] <- max(v[i], v[i - 1] + d)
    if (v[n] > hi) {
      v[n] <- hi
      for (i in (n - 1):1) v[i] <- min(v[i], v[i + 1] - d)
    }
    v
  }
  left <- pass(y)
  right <- rev(-pass(rev(-y)))   # the same from the right end
  if (right[1] < lo) {
    right[1] <- lo
    for (i in 2:n) right[i] <- max(right[i], right[i - 1] + d)
  }
  out <- numeric(n)
  out[o] <- (left + right) / 2
  out
}

# Gene model on one axis with variants as lollipops ("balloons"), in the
# style of geneHapR: CDS, 5'UTR and 3'UTR boxes on a line (introns and
# flanks), stems that fan out so that the allele labels do not overlap.
.gene_lollipop <- function(g, f, canon, vt, r0, r1, tr, xlab, colors, variant_colors, variant_labels,
                           base_size, title) {
  strand <- g$strand
  ex <- f[f$class == "exon", , drop = FALSE]
  cds <- f[f$class == "CDS", , drop = FALSE]
  boxes <- list()
  if (nrow(cds)) {
    lo <- min(cds$start)
    hi <- max(cds$end)
    boxes[[1]] <- data.frame(start = cds$start, end = cds$end, part = "CDS")
    for (k in seq_len(nrow(ex))) {
      s0 <- ex$start[k]
      e0 <- ex$end[k]
      if (s0 < lo) boxes[[length(boxes) + 1L]] <- data.frame(start = s0, end = min(e0, lo - 1),
                                                              part = if (strand == "-") "UTR3" else "UTR5")
      if (e0 > hi) boxes[[length(boxes) + 1L]] <- data.frame(start = max(s0, hi + 1), end = e0,
                                                              part = if (strand == "-") "UTR5" else "UTR3")
    }
  } else if (nrow(ex)) {
    boxes[[1]] <- data.frame(start = ex$start, end = ex$end, part = "exon")
  }
  bx <- do.call(rbind, boxes)
  if (is.null(bx)) bx <- data.frame(start = g$start, end = g$end, part = "exon")
  bx <- bx[bx$end >= bx$start, , drop = FALSE]
  bx$h <- ifelse(bx$part == "CDS", 0.2, 0.12)
  part_lab <- c(CDS = "CDS", UTR5 = "5' UTR", UTR3 = "3' UTR", exon = "exon (non-coding)")
  brk <- intersect(c("CDS", "UTR5", "UTR3", "exon"), bx$part)

  x0 <- tr(r0)
  x1 <- tr(r1)
  w <- x1 - x0
  tss <- if (strand == "-") g$end else g$start
  dir <- if (strand == "-") -1 else 1
  # Layers from back to front: gene axis, stems and TSS arrow (both start on
  # the axis), boxes (which cover the stem parts inside exons), balloons, labels.
  p <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = x0, xend = x1, y = 0, yend = 0, colour = "grey30", linewidth = 0.5)
  ytop <- 0.6
  if (!is.null(vt) && nrow(vt)) {
    vt$location <- .variant_location(vt$pos, g, f, strand)
    vt$type[!vt$type %in% names(variant_colors)] <- "indel"
    vt$alleles <- .allele_label(vt$ref, vt$alt)
    vt$label <- switch(variant_labels,
                       alleles = ifelse(is.na(vt$alleles), vt$id, vt$alleles),
                       id = vt$id,
                       position = .fmt(vt$pos),
                       none = "")
    vt <- vt[order(vt$pos), , drop = FALSE]
    vt$x <- tr(vt$pos)
    vt$xs <- .spread_positions(vt$x, d = w / 45, lo = x0 + 0.01 * w, hi = x1 - 0.01 * w)
    y1 <- 0.55
    y2 <- 1.15
    y3 <- 1.35
    stems <- rbind(data.frame(x = vt$x, xend = vt$x, y = 0, yend = y1),
                   data.frame(x = vt$x, xend = vt$xs, y = y1, yend = y2),
                   data.frame(x = vt$xs, xend = vt$xs, y = y2, yend = y3))
    p <- p + ggplot2::geom_segment(data = stems, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                                                              yend = .data$yend),
                                   colour = "grey45", linewidth = 0.3)
  }
  p <- p +
    ggplot2::annotate("segment", x = tr(tss), xend = tr(tss), y = 0, yend = -0.42, linewidth = 0.4) +
    ggplot2::annotate("segment", x = tr(tss), xend = tr(tss) + dir * 0.05 * w, y = -0.42, yend = -0.42,
                      linewidth = 0.4, arrow = grid::arrow(length = grid::unit(1.3, "mm"), type = "closed", angle = 30)) +
    ggplot2::geom_rect(data = bx, ggplot2::aes(xmin = tr(.data$start), xmax = tr(.data$end + 1),
                                               ymin = -.data$h, ymax = .data$h, fill = .data$part),
                       colour = "grey15", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = colors, breaks = brk, labels = unname(part_lab[brk]), name = NULL)
  if (!is.null(vt) && nrow(vt)) {
    has_p <- !is.null(vt$logp) && any(!is.na(vt$logp))
    p <- p + ggplot2::geom_point(data = vt, ggplot2::aes(x = .data$xs, y = y3), shape = 21, size = 2.9,
                                 fill = "grey15", colour = "grey15")
    if (has_p) {
      p <- p + ggplot2::geom_point(data = vt, ggplot2::aes(x = .data$xs, y = y3, colour = .data$logp), size = 2.3) +
        ggplot2::scale_colour_gradient(low = "#FFF5EB", high = "#A50F15", na.value = "grey80",
                                       name = expression(-log[10](italic(P))))
    } else {
      p <- p + ggplot2::geom_point(data = vt, ggplot2::aes(x = .data$xs, y = y3, colour = .data$type), size = 2.3) +
        ggplot2::scale_colour_manual(values = variant_colors, name = NULL)
    }
    if (variant_labels != "none") {
      p <- p + ggplot2::geom_text(data = vt, ggplot2::aes(x = .data$xs, y = y3 + 0.12, label = .data$label),
                                  angle = 90, hjust = 0, vjust = 0.5, size = .txt(base_size, 0.75))
      ytop <- y3 + 0.12 + 0.11 * max(nchar(vt$label), 3)
    } else {
      ytop <- y3 + 0.2
    }
  }
  p <- p +
    ggplot2::coord_cartesian(xlim = c(x0, x1), ylim = c(-0.6, ytop), expand = FALSE, clip = "off") +
    ggplot2::labs(x = xlab, y = NULL) +
    theme_quard(base_size) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(), legend.position = "bottom",
                   plot.margin = ggplot2::margin(6, 8, 4, 8))
  if (title) {
    ttl <- sprintf("%s%s  |  chr %s:%s-%s (%s strand)  |  %s bp  |  transcript %s", g$gene_id,
                   if (!is.na(g$name) && g$name != g$gene_id) paste0(" (", g$name, ")") else "",
                   g$chr, .fmt(g$start), .fmt(g$end), g$strand, .fmt(g$end - g$start + 1), canon)
    sub <- if (!is.na(g$note) && nzchar(g$note)) g$note else NULL
    p <- p + ggplot2::labs(title = ttl, subtitle = sub) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = base_size + 1, face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = base_size - 0.5, colour = "grey30"))
  }
  if (!is.null(vt)) vt$x <- vt$xs <- vt$label <- NULL
  attr(p, "variants") <- vt
  attr(p, "gene") <- g
  p
}
