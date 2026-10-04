# Karyotype-style SNP / indel density maps ------------------------------------

#' Variant positions and types
#'
#' Reads only the position and allele columns of a HapMap or VCF file (fast,
#' even for whole-genome files) and classifies each variant as `SNP` (single
#' base alleles A/C/G/T), `indel` (HapMap `+`/`-` alleles, alleles of
#' different length, or symbolic `<INS>`/`<DEL>` VCF alleles) or `other`.
#'
#' @param x HapMap or VCF file, a `qp_geno` object, or a data frame with
#'   `chr`, `pos` and `type`.
#' @param chr_map Optional chromosome name mapping.
#' @param verbose Print a summary.
#' @return Data frame with `chr`, `pos` and `type` (`SNP`, `indel`, `other`).
#' @export
variant_sites <- function(x, chr_map = NULL, verbose = TRUE) {
  lv <- c("SNP", "indel", "other")
  if (inherits(x, "qp_geno")) {
    return(data.frame(chr = x$info$chr, pos = x$info$pos,
                      type = factor(ifelse(x$info$is_snp, "SNP", "indel"), levels = lv)))
  }
  if (is.data.frame(x)) {
    if (!all(c("chr", "pos") %in% names(x))) .stop("A data frame of variants needs columns chr and pos.")
    ty <- if (!is.null(x$type)) as.character(x$type) else "SNP"
    return(data.frame(chr = normalize_chr(x$chr, chr_map), pos = as.numeric(x$pos),
                      type = factor(ty, levels = lv)))
  }
  .check_file(x, "Variant")
  first <- readLines(x, n = 1L, warn = FALSE)
  if (startsWith(first, "##") || startsWith(first, "#CHROM")) {
    d <- .fread_df(x, skip = "#CHROM", header = TRUE, select = c(1L, 2L, 4L, 5L), colClasses = "character")
    chr <- d[[1]]
    pos <- suppressWarnings(as.numeric(d[[2]]))
    ref <- toupper(d[[3]])
    alts <- strsplit(toupper(d[[4]]), ",", fixed = TRUE)
    snp <- grepl("^[ACGT]$", ref) & vapply(alts, function(a) all(grepl("^[ACGT*]$", a)) && any(a != "*"), NA)
    indel <- !snp & vapply(seq_along(alts), function(i) {
      a <- alts[[i]]
      any(a %in% c("<INS>", "<DEL>", "<INDEL>", "*")) || any(grepl("^[ACGTN]+$", a) & nchar(a) != nchar(ref[i]))
    }, NA)
  } else if (.is_hapmap_file(x)) {
    lay <- .hapmap_layout(x)
    d <- .fread_df(x, header = TRUE, select = lay$r + 1:3, colClasses = "character")
    al <- strsplit(toupper(d[[1]]), "/", fixed = TRUE)
    chr <- d[[2]]
    pos <- suppressWarnings(as.numeric(d[[3]]))
    snp <- vapply(al, function(a) length(a) >= 2L && all(a %in% c("A", "C", "G", "T")), NA)
    indel <- !snp & vapply(al, function(a) any(a %in% c("+", "-", "0")) ||
                             (all(grepl("^[ACGT]+$", a)) && length(unique(nchar(a))) > 1L), NA)
  } else {
    .stop("'", x, "' is neither a VCF nor a HapMap file.")
  }
  type <- ifelse(snp, "SNP", ifelse(indel, "indel", "other"))
  out <- data.frame(chr = normalize_chr(chr, chr_map), pos = pos, type = factor(type, levels = lv),
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$pos) & !is.na(out$chr), , drop = FALSE]
  tt <- table(out$type)
  .msg(verbose, "Variants: ", .fmt(tt[["SNP"]]), " SNPs, ", .fmt(tt[["indel"]]), " indels",
       if (tt[["other"]]) paste0(", ", .fmt(tt[["other"]]), " other") else "", " read from ", basename(x), ".")
  out
}

#' Variant density per genomic bin
#'
#' @param x Variants: anything accepted by [variant_sites()].
#' @param chrom_sizes Chromosome lengths (see [chrom_sizes()]); default:
#'   the species chosen with [set_species()], otherwise the largest variant
#'   position per chromosome.
#' @param bin_size Bin size in bp (default 1 Mb).
#' @param chr_map Optional chromosome name mapping.
#' @return Data frame with `chr`, `bin_start`, `bin_end`, `snp`, `indel`,
#'   `other` and `total` counts (attribute `sizes` holds the chromosome lengths).
#' @export
variant_density <- function(x, chrom_sizes = NULL, bin_size = 1e6, chr_map = NULL) {
  s <- if (is.data.frame(x) && all(c("chr", "pos", "type") %in% names(x))) x else
    variant_sites(x, chr_map = chr_map, verbose = FALSE)
  chrom_sizes <- chrom_sizes %||% .species_sizes()
  sizes <- if (!is.null(chrom_sizes)) {
    chrom_sizes(chrom_sizes, chr_map = chr_map)
  } else {
    chrom_sizes(NULL, gwas = data.frame(chr = s$chr, pos = s$pos))
  }
  s <- s[s$chr %in% sizes$chr, , drop = FALSE]
  res <- lapply(seq_len(nrow(sizes)), function(k) {
    L <- sizes$length[k]
    nb <- max(1L, ceiling(L / bin_size))
    sk <- s[s$chr == sizes$chr[k], , drop = FALSE]
    b <- pmin(nb, floor((sk$pos - 1) / bin_size) + 1L)
    cnt <- function(t) tabulate(b[sk$type == t], nbins = nb)
    data.frame(chr = sizes$chr[k], bin_start = (seq_len(nb) - 1) * bin_size + 1,
               bin_end = pmin(seq_len(nb) * bin_size, L), snp = cnt("SNP"), indel = cnt("indel"),
               other = cnt("other"), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res)
  out$total <- out$snp + out$indel + out$other
  attr(out, "sizes") <- sizes
  attr(out, "bin_size") <- bin_size
  out
}

.ramp_hex <- function(v, maxv, col) {
  t <- pmin(pmax(v / max(maxv, 1e-9), 0), 1)
  m <- grDevices::colorRamp(c("white", col))(t)
  grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
}

.rounded_rect <- function(x0, x1, y0, y1, rx, ry, n = 10L) {
  a <- seq(0, pi / 2, length.out = n)
  data.frame(x = c(x1 - rx + rx * cos(a), x0 + rx + rx * cos(a + pi / 2),
                   x0 + rx + rx * cos(a + pi), x1 - rx + rx * cos(a + 3 * pi / 2)),
             y = c(y1 - ry + ry * sin(a), y1 - ry + ry * sin(a + pi / 2),
                   y0 + ry + ry * sin(a + pi), y0 + ry + ry * sin(a + 3 * pi / 2)))
}

# White masks between the square corners and the rounded outline.
.corner_masks <- function(x0, x1, y0, y1, rx, ry, n = 10L) {
  a <- seq(0, pi / 2, length.out = n)
  rbind(
    data.frame(x = c(x1, x1 - rx + rx * cos(a)), y = c(y1, y1 - ry + ry * sin(a)), part = "tr"),
    data.frame(x = c(x0, x0 + rx + rx * cos(a + pi / 2)), y = c(y1, y1 - ry + ry * sin(a + pi / 2)), part = "tl"),
    data.frame(x = c(x0, x0 + rx + rx * cos(a + pi)), y = c(y0, y0 + ry + ry * sin(a + pi)), part = "bl"),
    data.frame(x = c(x1, x1 - rx + rx * cos(a + 3 * pi / 2)), y = c(y0, y0 + ry + ry * sin(a + 3 * pi / 2)), part = "br"))
}

.density_legend <- function(colors, max_snp, max_ind, unit_lab, base_size) {
  v <- seq(0, 1, length.out = 101)
  bars <- rbind(data.frame(x = v, y = 1, col = .ramp_hex(v, 1, colors[["SNP"]])),
                data.frame(x = v + 1.5, y = 1, col = .ramp_hex(v, 1, colors[["indel"]])))
  ticks <- data.frame(x = c(0, 0.5, 1, 1.5, 2, 2.5), y = 0.55,
                      lab = .fmt(round(c(0, max_snp / 2, max_snp, 0, max_ind / 2, max_ind))))
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = bars, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$col),
                       width = 0.0101, height = 0.5) +
    ggplot2::annotate("rect", xmin = c(-0.005, 1.495), xmax = c(1.005, 2.505), ymin = 0.75, ymax = 1.25,
                      fill = NA, colour = "grey40", linewidth = 0.2) +
    ggplot2::geom_text(data = ticks, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                       size = .txt(base_size, 0.8), vjust = 1) +
    ggplot2::annotate("text", x = c(0.5, 2), y = 1.5, label = c(paste("SNPs per", unit_lab), paste("Indels per", unit_lab)),
                      size = .txt(base_size, 0.9), fontface = "bold", vjust = 0) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(-0.4, 2.9), ylim = c(0, 1.9), expand = FALSE) +
    ggplot2::theme_void()
}

#' Karyotype-style SNP and indel density plot
#'
#' Draws every chromosome as an ideogram whose two halves are coloured by
#' variant density in each bin: SNPs (red, left/top half) and indels (yellow,
#' right/bottom half). Colour intensity is the number of variants per bin
#' (capped at the `cap` quantile so a few dense bins do not wash out the
#' rest). Works directly on whole-genome HapMap or VCF files.
#'
#' @param x HapMap or VCF file, `qp_geno` object or [variant_sites()] table.
#' @param chrom_sizes Chromosome lengths (see [chrom_sizes()]); default: the
#'   species chosen with [set_species()], otherwise the largest variant
#'   position per chromosome.
#' @param bin_size Bin size in bp (default 1 Mb).
#' @param colors Named colours for `SNP` and `indel`.
#' @param layout `"vertical"` (chromosomes side by side) or `"horizontal"`.
#' @param centromeres Optional data frame with `chr` and `pos` (bp) to mark
#'   centromeres.
#' @param show_counts Print SNP and indel totals for each chromosome.
#' @param cap Quantile of non-empty bins used as the top of each colour scale.
#' @param title Plot title (default describes the bin size).
#' @param chr_map Optional chromosome name mapping.
#' @param base_size Base font size (pt).
#' @param verbose Print a summary.
#' @return A patchwork object (plot + legend); attribute `density` holds the
#'   counts per bin. Save with [save_pdf()].
#' @examples
#' \donttest{
#' ex <- simulate_example_data(tempfile("qp"))
#' p <- plot_karyotype(ex$genome_hapmap, chrom_sizes = "rice", bin_size = 2e6)
#' save_pdf(p, tempfile(fileext = ".pdf"), width = 180, height = 120)
#' }
#' @export
plot_karyotype <- function(x, chrom_sizes = NULL, bin_size = 1e6,
                           colors = c(SNP = "#D7191C", indel = "#F2B705"),
                           layout = c("vertical", "horizontal"), centromeres = NULL, show_counts = TRUE,
                           cap = 0.99, title = NULL, chr_map = NULL, base_size = 8, verbose = TRUE) {
  layout <- match.arg(layout)
  s <- if (is.data.frame(x) && all(c("chr", "pos", "type") %in% names(x))) x else
    variant_sites(x, chr_map = chr_map, verbose = verbose)
  dens <- variant_density(s, chrom_sizes, bin_size, chr_map)
  sizes <- attr(dens, "sizes")
  K <- nrow(sizes)
  qmax <- function(v) {
    v <- v[v > 0]
    if (!length(v)) 1 else max(1, as.numeric(stats::quantile(v, cap)))
  }
  max_snp <- qmax(dens$snp)
  max_ind <- qmax(dens$indel)
  dens$col_snp <- .ramp_hex(dens$snp, max_snp, colors[["SNP"]])
  dens$col_ind <- .ramp_hex(dens$indel, max_ind, colors[["indel"]])
  ci <- match(dens$chr, sizes$chr)
  Lmb <- sizes$length / 1e6
  maxL <- max(Lmb)
  w <- 0.62
  unit_lab <- if (bin_size >= 1e6) paste(format(bin_size / 1e6), "Mb") else paste(format(bin_size / 1e3), "kb")
  tot <- data.frame(chr = sizes$chr, snp = vapply(sizes$chr, function(cc) sum(dens$snp[dens$chr == cc]), 0),
                    indel = vapply(sizes$chr, function(cc) sum(dens$indel[dens$chr == cc]), 0))
  if (layout == "vertical") {
    bins <- rbind(
      data.frame(xmin = ci - w / 2, xmax = ci, ymin = (dens$bin_start - 1) / 1e6, ymax = dens$bin_end / 1e6,
                 col = dens$col_snp),
      data.frame(xmin = ci, xmax = ci + w / 2, ymin = (dens$bin_start - 1) / 1e6, ymax = dens$bin_end / 1e6,
                 col = dens$col_ind))
    rx <- w / 2
    ry <- min(maxL * 0.02, min(Lmb) / 3)
    outl <- do.call(rbind, lapply(seq_len(K), function(k) {
      cbind(.rounded_rect(k - w / 2, k + w / 2, 0, Lmb[k], rx, ry), id = k)
    }))
    mask <- do.call(rbind, lapply(seq_len(K), function(k) {
      m <- .corner_masks(k - w / 2, k + w / 2, 0, Lmb[k], rx, ry)
      m$id <- paste(k, m$part)
      m
    }))
  } else {
    yk <- K - seq_len(K) + 1
    h <- 0.62
    bins <- rbind(
      data.frame(ymin = yk[ci], ymax = yk[ci] + h / 2, xmin = (dens$bin_start - 1) / 1e6, xmax = dens$bin_end / 1e6,
                 col = dens$col_snp),
      data.frame(ymin = yk[ci] - h / 2, ymax = yk[ci], xmin = (dens$bin_start - 1) / 1e6, xmax = dens$bin_end / 1e6,
                 col = dens$col_ind))
    rx <- min(maxL * 0.012, min(Lmb) / 3)
    ry <- h / 2
    outl <- do.call(rbind, lapply(seq_len(K), function(k) {
      cbind(.rounded_rect(0, Lmb[k], yk[k] - h / 2, yk[k] + h / 2, rx, ry), id = k)
    }))
    mask <- do.call(rbind, lapply(seq_len(K), function(k) {
      m <- .corner_masks(0, Lmb[k], yk[k] - h / 2, yk[k] + h / 2, rx, ry)
      m$id <- paste(k, m$part)
      m
    }))
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bins, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
                                                 ymax = .data$ymax, fill = .data$col), colour = NA) +
    ggplot2::geom_polygon(data = mask, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                          fill = "white", colour = NA) +
    ggplot2::geom_polygon(data = outl, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                          fill = NA, colour = "grey25", linewidth = 0.35) +
    ggplot2::scale_fill_identity()
  if (!is.null(centromeres)) {
    cen <- data.frame(chr = normalize_chr(centromeres[[1]], chr_map), pos = as.numeric(centromeres[[2]]) / 1e6)
    cen <- cen[cen$chr %in% sizes$chr, , drop = FALSE]
    if (nrow(cen)) {
      k <- match(cen$chr, sizes$chr)
      tri <- if (layout == "vertical") {
        do.call(rbind, lapply(seq_along(k), function(z) {
          x0 <- k[z]
          y0 <- cen$pos[z]
          d <- maxL * 0.012
          data.frame(id = rep(c(paste0(z, "l"), paste0(z, "r")), each = 3),
                     x = c(x0 - w / 2 - 0.12, x0 - w / 2 - 0.12, x0 - w / 2, x0 + w / 2 + 0.12, x0 + w / 2 + 0.12, x0 + w / 2),
                     y = c(y0 - d, y0 + d, y0, y0 - d, y0 + d, y0))
        }))
      } else {
        yk <- K - seq_len(K) + 1
        do.call(rbind, lapply(seq_along(k), function(z) {
          y0 <- yk[k[z]]
          x0 <- cen$pos[z]
          d <- maxL * 0.008
          data.frame(id = rep(c(paste0(z, "t"), paste0(z, "b")), each = 3),
                     x = c(x0 - d, x0 + d, x0, x0 - d, x0 + d, x0),
                     y = c(y0 + 0.31 + 0.12, y0 + 0.31 + 0.12, y0 + 0.31, y0 - 0.31 - 0.12, y0 - 0.31 - 0.12, y0 - 0.31))
        }))
      }
      p <- p + ggplot2::geom_polygon(data = tri, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                                     fill = "grey20", colour = NA)
    }
  }
  if (layout == "vertical") {
    if (show_counts) {
      lab <- data.frame(x = seq_len(K), y = Lmb + maxL * 0.02,
                        snp = paste0(.fmt(tot$snp)), ind = paste0(.fmt(tot$indel)))
      p <- p +
        ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$x, y = .data$y, label = .data$snp), vjust = 1,
                           size = .txt(base_size, 0.7), colour = "#9E0E10") +
        ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$x, y = .data$y + maxL * 0.035, label = .data$ind),
                           vjust = 1, size = .txt(base_size, 0.7), colour = "#8C6A00")
    }
    p <- p +
      ggplot2::scale_x_continuous(breaks = seq_len(K), labels = sizes$chr, position = "top",
                                  expand = ggplot2::expansion(add = 0.5)) +
      ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = c(0.02, if (show_counts) 0.09 else 0.02))) +
      ggplot2::labs(x = "Chromosome", y = "Position (Mb)")
  } else {
    yk <- K - seq_len(K) + 1
    if (show_counts) {
      lab <- data.frame(x = Lmb + maxL * 0.01, y = yk,
                        lab = paste0(.fmt(tot$snp), " SNPs / ", .fmt(tot$indel), " indels"))
      p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab), hjust = 0,
                                  size = .txt(base_size, 0.7), colour = "grey25")
    }
    p <- p +
      ggplot2::scale_y_continuous(breaks = yk, labels = sizes$chr, expand = ggplot2::expansion(add = 0.5)) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, if (show_counts) 0.18 else 0.02))) +
      ggplot2::labs(y = "Chromosome", x = "Position (Mb)")
  }
  p <- p + theme_quard(base_size) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_line(linewidth = 0.25))
  if (layout == "vertical") p <- p + ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
  if (layout == "horizontal") p <- p + ggplot2::theme(axis.ticks.y = ggplot2::element_blank())
  ttl <- title %||% sprintf("SNP and indel density (%s bins): %s SNPs, %s indels", unit_lab,
                            .fmt(sum(tot$snp)), .fmt(sum(tot$indel)))
  p <- p + ggplot2::ggtitle(ttl) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = base_size + 1))
  leg <- .density_legend(colors, max_snp, max_ind, unit_lab, base_size)
  out <- patchwork::wrap_plots(p, leg, ncol = 1, heights = c(1, 0.11))
  attr(out, "density") <- dens[, c("chr", "bin_start", "bin_end", "snp", "indel", "other", "total")]
  out
}
