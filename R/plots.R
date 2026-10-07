# Plot panels (ggplot2) ---------------------------------------------------------

#' ggplot2 theme used by quardplot
#'
#' @param base_size Base font size in points (8 pt suits a 170-180 mm wide
#'   journal figure).
#' @param base_family Font family.
#' @return A ggplot2 theme.
#' @export
theme_quard <- function(base_size = 8, base_family = "") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.3, colour = "grey20"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "grey20"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size * 0.9),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.title = ggplot2::element_text(face = "bold", size = base_size * 0.9),
      legend.text = ggplot2::element_text(size = base_size * 0.85),
      legend.key.size = grid::unit(3, "mm"),
      legend.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 6, 3, 3)
    )
}

.no_x <- function() {
  ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank())
}

.no_y <- function() {
  ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                 axis.line.y = ggplot2::element_blank())
}

.txt <- function(base_size, rel = 0.85) base_size * rel / ggplot2::.pt

# Remove overplotted points below `keep_above` (-log10 P) on a fine grid,
# which shrinks the PDF without changing what is visible.
.thin_points <- function(d, ymax, keep_above = 2, nx = 1500L, ny = 300L) {
  low <- d$logp < keep_above
  if (sum(low) < 20000L) return(d)
  rx <- range(d$x)
  bx <- floor((d$x[low] - rx[1]) / max(diff(rx), 1) * nx)
  by <- floor(d$logp[low] / max(ymax, 1e-9) * ny)
  key <- paste(d$grp[low], bx, by)
  rbind(d[!low, , drop = FALSE], d[low, , drop = FALSE][!duplicated(key), , drop = FALSE])
}

#' Genome-wide Manhattan plot
#'
#' @param gwas GWAS data ([read_gwas()]) or a path.
#' @param chrom_sizes Chromosome lengths: anything accepted by
#'   [chrom_sizes()] (preset name, `.fai`/size file, VCF, named vector), or
#'   `NULL`: the species chosen with [set_species()], otherwise estimated
#'   from the GWAS positions.
#' @param threshold Significance line: `"bonferroni"`, a P value, a -log10
#'   value, or `NULL` (see [gwas_threshold()]).
#' @param suggestive Optional second (dotted) line, same format.
#' @param highlight Optional region (e.g. from [define_region()]) to shade.
#' @param target Optional SNP ID to mark and label.
#' @param colors Two alternating chromosome colours.
#' @param point_size Point size (mm).
#' @param thin Remove invisible overplotted points with -log10 P < 2 to keep
#'   the PDF small.
#' @param base_size Base font size (pt).
#' @return A ggplot object.
#' @export
plot_manhattan <- function(gwas, chrom_sizes = NULL, threshold = "bonferroni", suggestive = NULL,
                           highlight = NULL, target = NULL, colors = c("#3B6E9C", "#A3A9AE"),
                           point_size = 0.5, thin = TRUE, base_size = 8) {
  gwas <- .as_gwas(gwas)
  chrom_sizes <- chrom_sizes %||% .species_sizes()
  sizes <- if (is.data.frame(chrom_sizes) && all(c("chr", "length") %in% names(chrom_sizes))) {
    chrom_sizes
  } else {
    chrom_sizes(chrom_sizes, gwas = gwas)
  }
  lay <- genome_layout(sizes)
  d <- data.frame(snp = gwas$snp, chr = gwas$chr, pos = gwas$pos, logp = gwas$logp,
                  stringsAsFactors = FALSE)
  li <- match(d$chr, lay$chr)
  if (anyNA(li)) {
    .warn(.fmt(sum(is.na(li))), " SNP(s) on chromosomes missing from `chrom_sizes` were not plotted (",
          paste(utils::head(unique(d$chr[is.na(li)]), 5), collapse = ", "), ").")
    d <- d[!is.na(li), , drop = FALSE]
    li <- li[!is.na(li)]
  }
  beyond <- d$pos > lay$length[li]
  if (any(beyond)) {
    .warn(.fmt(sum(beyond)), " SNP position(s) exceed the chromosome lengths in `chrom_sizes`: ",
          "check that the GWAS and the chromosome sizes use the same genome assembly.")
  }
  d$x <- d$pos + lay$offset[li]
  d$grp <- factor(lay$index[li] %% 2L, levels = c(1L, 0L))
  th <- gwas_threshold(gwas, threshold)
  sg <- if (is.null(suggestive)) NULL else gwas_threshold(gwas, suggestive)
  ymax <- max(c(d$logp, th$logp, sg$logp), na.rm = TRUE) * if (is.null(target)) 1.06 else 1.14
  if (thin) d <- .thin_points(d, ymax)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$logp, colour = .data$grp))
  if (!is.null(highlight)) {
    hr <- .as_region(highlight)
    off <- lay$offset[match(hr$chr, lay$chr)]
    if (!is.na(off)) {
      w <- max(hr$end - hr$start, 0.006 * max(lay$end_cum))
      mid <- (hr$start + hr$end) / 2 + off
      p <- p + ggplot2::annotate("rect", xmin = mid - w / 2, xmax = mid + w / 2, ymin = 0, ymax = ymax,
                                 fill = "#E69F00", alpha = 0.25)
    }
  }
  p <- p +
    ggplot2::geom_point(size = point_size, shape = 16, stroke = 0) +
    ggplot2::scale_colour_manual(values = c(`1` = colors[1], `0` = colors[2]), guide = "none") +
    ggplot2::scale_x_continuous(breaks = lay$center, labels = lay$chr, expand = c(0, 0),
                                limits = c(0, max(lay$end_cum))) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::labs(x = "Chromosome", y = expression(bold(-log[10](italic(P))))) +
    theme_quard(base_size)
  if (!is.null(th)) {
    p <- p + ggplot2::geom_hline(yintercept = th$logp, linetype = "dashed", colour = "#D55E00", linewidth = 0.4)
  }
  if (!is.null(sg)) {
    p <- p + ggplot2::geom_hline(yintercept = sg$logp, linetype = "dotted", colour = "grey40", linewidth = 0.4)
  }
  if (!is.null(target) && !is.na(target)) {
    k <- match(target, gwas$snp)
    if (!is.na(k) && gwas$chr[k] %in% lay$chr) {
      tp <- data.frame(x = gwas$pos[k] + lay$offset[match(gwas$chr[k], lay$chr)], y = gwas$logp[k],
                       lab = target)
      p <- p +
        ggplot2::geom_point(data = tp, ggplot2::aes(x = .data$x, y = .data$y), inherit.aes = FALSE,
                            shape = 23, fill = "#7B3294", colour = "black", stroke = 0.3,
                            size = point_size * 3.5) +
        ggplot2::geom_text(data = tp, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                           inherit.aes = FALSE, vjust = -0.9, size = .txt(base_size), fontface = "bold")
    }
  }
  attr(p, "quard_layout") <- lay
  attr(p, "quard_panel") <- "manhattan"
  p
}

#' Connector lines from the genome-wide plot to the regional panels
#'
#' @param manhattan A plot from [plot_manhattan()].
#' @param region The region shown in the panels below.
#' @param colour Line colour.
#' @return A ggplot object to place between the Manhattan plot and the
#'   regional panels (see [combine_panels()]).
#' @export
plot_zoom_connector <- function(manhattan, region, colour = "grey45") {
  lay <- attr(manhattan, "quard_layout")
  if (is.null(lay)) .stop("`manhattan` must be made by plot_manhattan().")
  r <- .as_region(region)
  total <- max(lay$end_cum)
  off <- lay$offset[match(r$chr, lay$chr)]
  if (is.na(off)) .stop("Chromosome ", r$chr, " is not in the Manhattan plot.")
  w <- max(r$end - r$start, 0.006 * total)
  mid <- (r$start + r$end) / 2 + off
  d <- data.frame(x = c(mid - w / 2, mid + w / 2) / total, xend = c(0, 1), y = 1, yend = 0)
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
                          colour = colour, linewidth = 0.35, linetype = "22") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 6, 0, 3))
  attr(p, "quard_panel") <- "connector"
  p
}

#' Regional association plot
#'
#' -log10 P of the variants in the region, coloured by \eqn{r^2} with the
#' lead variant (colour-blind-safe sequential palette; crosses = no genotype
#' data). The lead variant is drawn as a purple diamond and labelled. The
#' legend is drawn in one row above the panel.
#'
#' @param gwas GWAS data ([read_gwas()]); the threshold uses all its SNPs.
#' @param region Region list (`chr`, `start`, `end`), e.g. [define_region()].
#' @param lead Lead SNP ID (default: most significant SNP in the region).
#' @param r2 Optional named vector of \eqn{r^2} with the lead SNP, names =
#'   GWAS SNP IDs (or a two-column data frame: SNP, r2).
#' @param threshold See [gwas_threshold()].
#' @param point_size Point size (mm).
#' @param label_lead Label the lead SNP.
#' @param r2_colors Five colours for \eqn{r^2} bins 0-0.2, ..., 0.8-1.
#' @param lead_color Fill colour of the lead SNP.
#' @param base_size Base font size (pt).
#' @return A ggplot object.
#' @export
plot_regional <- function(gwas, region, lead = NULL, r2 = NULL, threshold = "bonferroni",
                          point_size = 1.6, label_lead = TRUE,
                          r2_colors = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026"),
                          lead_color = "#7B3294", base_size = 8) {
  gwas <- .as_gwas(gwas)
  r <- .as_region(region)
  th <- gwas_threshold(gwas, threshold)
  keep <- gwas$chr == r$chr & gwas$pos >= r$start & gwas$pos <= r$end
  d <- data.frame(snp = gwas$snp[keep], pos = gwas$pos[keep], p = gwas$p[keep],
                  logp = gwas$logp[keep], stringsAsFactors = FALSE)
  if (!nrow(d)) {
    .stop("No GWAS SNPs in chr ", r$chr, ":", .fmt(r$start), "-", .fmt(r$end),
          ". Check the coordinates and that GWAS and annotation use the same assembly.")
  }
  if (is.null(lead) || is.na(lead) || !(lead %in% d$snp)) lead <- d$snp[which.min(d$p)]
  d$r2 <- NA_real_
  if (!is.null(r2)) {
    if (is.data.frame(r2)) r2 <- stats::setNames(r2[[2]], r2[[1]])
    d$r2 <- unname(r2[d$snp])
  }
  bins <- c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")
  b <- as.character(cut(d$r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = bins))
  b[is.na(b)] <- "no LD data"
  d$bin <- factor(b, levels = c(rev(bins), "no LD data"))
  pal <- c(stats::setNames(rev(r2_colors), rev(bins)), `no LD data` = "grey55")
  # SNPs without LD data are crosses, so they stay distinct when printed in grey
  shp <- c(stats::setNames(rep(21L, length(bins)), rev(bins)), `no LD data` = 4L)
  d <- d[order(!is.na(d$r2), d$r2), , drop = FALSE]
  ld <- d[d$snp == lead, , drop = FALSE]
  ymax <- max(c(d$logp, th$logp), na.rm = TRUE) * 1.18
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$pos / 1e6, y = .data$logp))
  if (!is.null(th)) {
    p <- p + ggplot2::geom_hline(yintercept = th$logp, linetype = "dashed", colour = "#D55E00", linewidth = 0.4)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(fill = .data$bin, shape = .data$bin), colour = "grey30", stroke = 0.3,
                        size = point_size) +
    ggplot2::geom_point(data = ld, shape = 23, fill = lead_color, colour = "black", stroke = 0.3,
                        size = point_size * 1.7) +
    ggplot2::scale_fill_manual(values = pal, name = expression(bold(italic(r)^2)), drop = TRUE) +
    ggplot2::scale_shape_manual(values = shp, name = expression(bold(italic(r)^2)), drop = TRUE) +
    ggplot2::coord_cartesian(xlim = c(r$start, r$end) / 1e6, ylim = c(-0.02 * ymax, ymax), expand = FALSE) +
    ggplot2::labs(x = sprintf("Chromosome %s (Mb)", r$chr), y = expression(bold(-log[10](italic(P))))) +
    theme_quard(base_size) +
    # legend in one row above the panel, so it never hides SNPs
    ggplot2::theme(legend.position = "top", legend.justification = c(1, 0), legend.direction = "horizontal",
                   legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.margin = ggplot2::margin(0, 0, -4, 0),
                   legend.key.width = grid::unit(2.5, "mm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 2.2), nrow = 1),
                    shape = ggplot2::guide_legend(nrow = 1))
  if (label_lead && nrow(ld)) {
    rel <- (ld$pos[1] - r$start) / (r$end - r$start)
    hj <- if (rel < 0.15) 0 else if (rel > 0.85) 1 else 0.5
    p <- p + ggplot2::geom_text(data = ld, ggplot2::aes(label = .data$snp), vjust = -1.2, hjust = hj,
                                size = .txt(base_size), fontface = "bold")
  }
  if (all(is.na(d$r2))) p <- p + ggplot2::theme(legend.position = "none")
  attr(p, "quard_panel") <- "regional"
  attr(p, "quard_region") <- r
  p
}

#' Gene-structure track
#'
#' Draws the genes of a region with exons, coding sequence (tall boxes),
#' UTRs (short boxes) and introns (lines); colour and the arrow at the 3'
#' end give the strand. Genes are placed in non-overlapping rows, leaving
#' room for their (italic) labels.
#'
#' @param genes Output of [read_gff()].
#' @param region Region list (`chr`, `start`, `end`).
#' @param transcripts `"canonical"` (one representative transcript per gene)
#'   or `"all"`.
#' @param label Label genes (gene name, or transcript ID with
#'   `transcripts = "all"`).
#' @param label_size Label size in mm (default from `base_size`).
#' @param strand_colors Colours for `+` and `-` strand genes.
#' @param highlight Optional gene IDs or names to emphasise.
#' @param highlight_color Colour of highlighted genes.
#' @param base_size Base font size (pt).
#' @return A ggplot object (attribute `quard_rows` gives the number of rows).
#' @export
plot_genes <- function(genes, region, transcripts = c("canonical", "all"), label = TRUE,
                       label_size = NULL, strand_colors = c(`+` = "#0072B2", `-` = "#D55E00"),
                       highlight = NULL, highlight_color = "#CC79A7", base_size = 8) {
  transcripts <- match.arg(transcripts)
  if (!inherits(genes, "qp_genes")) .stop("`genes` must come from read_gff().")
  r <- .as_region(region)
  sub <- genes_in_region(genes, r)
  label_size <- label_size %||% .txt(base_size, 0.8)
  rng <- c(r$start, r$end) / 1e6
  xl <- sprintf("Chromosome %s (Mb)", r$chr)
  if (!nrow(sub$genes)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = mean(rng), y = 0, label = "No annotated genes in this region",
                        size = label_size, colour = "grey40") +
      ggplot2::coord_cartesian(xlim = rng, ylim = c(-1, 1), expand = FALSE) +
      ggplot2::labs(x = xl, y = "Genes") + theme_quard(base_size) + .no_y()
    attr(p, "quard_panel") <- "genes"
    attr(p, "quard_rows") <- 1L
    return(p)
  }
  tx <- sub$transcripts
  if (transcripts == "canonical") tx <- tx[tx$canonical, , drop = FALSE]
  gi <- match(tx$gene_id, sub$genes$gene_id)
  tx$name <- sub$genes$name[gi]
  tx$label <- if (transcripts == "all") ifelse(tx$type == "gene_only", tx$name, tx$transcript_id) else tx$name
  tx$strand[!tx$strand %in% c("+", "-")] <- "."
  tx$hl <- !is.null(highlight) & (tx$gene_id %in% highlight | tx$name %in% highlight)
  w <- r$end - r$start
  char_bp <- w * 0.0075 * (label_size / 2.3)
  tx <- tx[order(tx$start, tx$end), , drop = FALSE]
  lab_w <- if (label) nchar(tx$label) * char_bp else 0
  mid <- (pmax(tx$start, r$start) + pmin(tx$end, r$end)) / 2
  left <- pmin(tx$start, mid - lab_w / 2)
  right <- pmax(tx$end, mid + lab_w / 2)
  pad <- w * 0.012
  row_end <- numeric()
  tx$row <- 0L
  for (k in seq_len(nrow(tx))) {
    free <- which(row_end < left[k] - pad)
    rr <- if (length(free)) free[1] else length(row_end) + 1L
    row_end[rr] <- right[k]
    tx$row[k] <- rr
  }
  nrows <- max(tx$row)
  tx$y <- -tx$row
  tx$col <- ifelse(tx$hl, "hl", tx$strand)
  f <- sub$features[sub$features$transcript_id %in% tx$transcript_id, , drop = FALSE]
  ti <- match(f$transcript_id, tx$transcript_id)
  f$y <- tx$y[ti]
  f$col <- tx$col[ti]
  coding <- unique(f$transcript_id[f$class == "CDS"])
  ex <- f[f$class == "exon", , drop = FALSE]
  ex$h <- ifelse(ex$transcript_id %in% coding, 0.13, 0.22)
  cds <- f[f$class == "CDS", , drop = FALSE]
  arr_len <- w * 0.015
  arr <- tx[tx$strand %in% c("+", "-"), , drop = FALSE]
  arr$x0 <- ifelse(arr$strand == "+", arr$end, arr$start)
  arr$x1 <- ifelse(arr$strand == "+", arr$end + arr_len, arr$start - arr_len)
  tx$lx <- pmin(pmax(mid, r$start + lab_w / 2), r$end - lab_w / 2)
  pal <- c(strand_colors, `.` = "grey45", hl = highlight_color)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = tx, ggplot2::aes(x = .data$start / 1e6, xend = .data$end / 1e6,
                                                  y = .data$y, yend = .data$y, colour = .data$col),
                          linewidth = 0.35) +
    ggplot2::geom_rect(data = ex, ggplot2::aes(xmin = .data$start / 1e6, xmax = (.data$end + 1) / 1e6,
                                               ymin = .data$y - .data$h, ymax = .data$y + .data$h,
                                               fill = .data$col), colour = NA)
  if (nrow(cds)) {
    p <- p + ggplot2::geom_rect(data = cds, ggplot2::aes(xmin = .data$start / 1e6, xmax = (.data$end + 1) / 1e6,
                                                         ymin = .data$y - 0.26, ymax = .data$y + 0.26,
                                                         fill = .data$col), colour = NA)
  }
  if (nrow(arr)) {
    p <- p + ggplot2::geom_segment(data = arr, ggplot2::aes(x = .data$x0 / 1e6, xend = .data$x1 / 1e6,
                                                            y = .data$y, yend = .data$y, colour = .data$col),
                                   linewidth = 0.35,
                                   arrow = grid::arrow(length = grid::unit(1.1, "mm"), angle = 30,
                                                       type = "closed"))
  }
  if (label) {
    p <- p + ggplot2::geom_text(data = tx, ggplot2::aes(x = .data$lx / 1e6, y = .data$y - 0.47,
                                                        label = .data$label),
                                size = label_size, fontface = "italic", vjust = 1)
  }
  p <- p +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::coord_cartesian(xlim = rng, ylim = c(-nrows - 0.95, -0.55), expand = FALSE) +
    ggplot2::labs(x = xl, y = "Genes") +
    theme_quard(base_size) + .no_y()
  attr(p, "quard_panel") <- "genes"
  attr(p, "quard_rows") <- nrows
  p
}

#' Pairwise LD heatmap
#'
#' Draws the upper triangle of an LD matrix rotated by 45 degrees, with
#' connector lines from each variant's genomic position to its column, so the
#' heatmap aligns with the other regional panels.
#'
#' @param ld An LD matrix from [calc_ld()] or [read_plink_ld()] (needs the
#'   `pos` attribute), or a `qp_geno` object (LD is then computed).
#' @param region Region list; defaults to the span of the variants.
#' @param lead Optional lead variant ID to highlight.
#' @param max_snps Maximum number of variants drawn; larger sets are thinned
#'   to evenly spaced variants (the lead is always kept). About 300 is the
#'   practical upper limit for a readable vector figure.
#' @param lead_label Text of the label of `lead` (default: its ID).
#' @param label_snps `"lead"`, `"none"` or `"all"` (only sensible for up to
#'   about 60 variants).
#' @param blocks Optional block table (`start`, `end`) to outline.
#' @param colors Colour ramp from 0 to 1 (default: colour-blind-safe,
#'   lightness decreases monotonically, so it also reads in greyscale).
#' @param base_size Base font size (pt).
#' @return A ggplot object.
#' @export
plot_ld <- function(ld, region = NULL, lead = NULL, max_snps = 300, label_snps = c("lead", "none", "all"),
                    blocks = NULL, colors = c("#FFFFFF", "#FEE391", "#FE9929", "#D94701", "#7F2704"),
                    base_size = 8, lead_label = NULL) {
  label_snps <- match.arg(label_snps)
  if (inherits(ld, "qp_geno")) ld <- calc_ld(ld)
  pos <- attr(ld, "pos")
  if (is.null(pos)) .stop("The LD matrix has no 'pos' attribute (use calc_ld() or read_plink_ld(bim = ...)).")
  stat <- attr(ld, "stat") %||% "r2"
  ids <- rownames(ld) %||% paste0("v", seq_len(nrow(ld)))
  r <- if (!is.null(region)) .as_region(region) else
    .as_region(list(chr = attr(ld, "chr") %||% "", start = min(pos), end = max(pos)))
  keep <- which(pos >= r$start & pos <= r$end)
  keep <- keep[order(pos[keep])]
  n <- length(keep)
  if (n < 2L) .stop("Fewer than two variants with LD values in the region.")
  if (n > max_snps) {
    li <- if (!is.null(lead)) match(lead, ids[keep]) else NA
    sel <- sort(unique(c(round(seq(1, n, length.out = max_snps)), stats::na.omit(li))))
    message("LD heatmap: ", n, " variants thinned to ", length(sel), " evenly spaced variants (max_snps = ",
            max_snps, ").")
    keep <- keep[sel]
    n <- length(keep)
  }
  L <- ld[keep, keep, drop = FALSE]
  pos <- pos[keep]
  ids <- ids[keep]
  x0 <- r$start / 1e6
  x1 <- r$end / 1e6
  s <- (x1 - x0) / n
  xk <- x0 + (seq_len(n) - 0.5) * s
  h <- 1
  ij <- which(upper.tri(L), arr.ind = TRUE)
  i <- ij[, 1]
  j <- ij[, 2]
  cx <- (xk[i] + xk[j]) / 2
  cy <- -(j - i) * h
  poly <- data.frame(id = rep(seq_along(i), each = 4L),
                     x = as.vector(rbind(cx, cx + s / 2, cx, cx - s / 2)),
                     y = as.vector(rbind(cy + h, cy, cy - h, cy)),
                     value = rep(L[ij], each = 4L))
  lab_band <- if (label_snps == "all") max(3, 0.3 * n) else 0
  ct <- lab_band + max(2, 0.14 * n)
  conn <- data.frame(gx = pos / 1e6, xk = xk, y0 = ct, y1 = if (lab_band > 0) lab_band else -h)
  lead_i <- if (!is.null(lead)) match(lead, ids) else NA
  top_pad <- if (!is.na(lead_i) && label_snps != "none") max(1.5, 0.08 * n) else 0.5
  fill_name <- if (stat == "dprime") "D'" else expression(bold(italic(r)^2))
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = poly, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                                    fill = .data$value),
                          colour = if (n <= 80) "white" else NA, linewidth = 0.1) +
    ggplot2::annotate("segment", x = x0, xend = x1, y = ct, yend = ct, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_segment(data = conn, ggplot2::aes(x = .data$gx, xend = .data$xk, y = .data$y0,
                                                    yend = .data$y1),
                          colour = "grey60", linewidth = 0.2) +
    ggplot2::scale_fill_gradientn(colours = colors, limits = c(0, 1), na.value = "grey85",
                                  name = fill_name, breaks = c(0, 0.5, 1))
  if (label_snps == "all") {
    p <- p + ggplot2::annotate("text", x = xk, y = lab_band - 0.2, label = ids, angle = 90, hjust = 1,
                               size = .txt(base_size, 0.6))
  }
  if (!is.null(blocks) && nrow(blocks)) {
    bl <- if (!is.null(blocks$chr)) blocks[blocks$chr == r$chr, , drop = FALSE] else blocks
    seg <- lapply(seq_len(nrow(bl)), function(k) {
      a <- which(pos >= bl$start[k])[1]
      b <- utils::tail(which(pos <= bl$end[k]), 1)
      if (is.na(a) || !length(b) || b <= a) return(NULL)
      xa <- xk[a]
      xb <- xk[b]
      data.frame(x = c(xa, (xa + xb) / 2), xend = c((xa + xb) / 2, xb),
                 y = c(-h, -(b - a + 1) * h), yend = c(-(b - a + 1) * h, -h))
    })
    seg <- do.call(rbind, seg)
    if (!is.null(seg)) {
      p <- p + ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                                                              yend = .data$yend),
                                     colour = "black", linewidth = 0.4)
    }
  }
  if (!is.na(lead_i)) {
    lc <- conn[lead_i, , drop = FALSE]
    p <- p + ggplot2::geom_segment(data = lc, ggplot2::aes(x = .data$gx, xend = .data$xk, y = .data$y0,
                                                           yend = .data$y1),
                                   colour = "#7B3294", linewidth = 0.6)
    if (label_snps != "none") {
      rel <- (lc$gx - x0) / (x1 - x0)
      hj <- if (rel < 0.15) 0 else if (rel > 0.85) 1 else 0.5
      p <- p + ggplot2::annotate("text", x = lc$gx, y = ct + 0.25, label = if (is.null(lead_label)) ids[lead_i] else lead_label, vjust = 0, hjust = hj,
                                 size = .txt(base_size, 0.8), fontface = "bold", colour = "#7B3294")
    }
  }
  p <- p +
    ggplot2::coord_cartesian(xlim = c(x0, x1), ylim = c(-n * h - 0.5, ct + top_pad), expand = FALSE,
                             clip = "off") +
    ggplot2::labs(x = sprintf("Chromosome %s (Mb)", r$chr), y = NULL) +
    theme_quard(base_size) + .no_y() +
    .legend_inside(0.995, 0.03, c(1, 0)) +
    ggplot2::theme(legend.direction = "horizontal", legend.key.width = grid::unit(7, "mm"),
                   legend.key.height = grid::unit(2.2, "mm"))
  attr(p, "quard_panel") <- "ld"
  attr(p, "quard_n") <- n
  p
}

#' Custom annotation track
#'
#' Draws user intervals (for example QTL, peaks, candidate regions) on the
#' regional coordinate axis so they can be stacked with the other panels.
#'
#' @param data Data frame with `start` and `end` (bp) and optionally `chr`,
#'   `label` and `group`.
#' @param region Region list (`chr`, `start`, `end`).
#' @param title Y-axis title.
#' @param fill Fill colour (used when there is no `group` column).
#' @param base_size Base font size (pt).
#' @return A ggplot object.
#' @export
plot_track <- function(data, region, title = "Track", fill = "#009E73", base_size = 8) {
  r <- .as_region(region)
  d <- as.data.frame(data)
  if (!all(c("start", "end") %in% names(d))) .stop("`data` needs columns start and end.")
  if (!is.null(d$chr)) d <- d[normalize_chr(d$chr) == r$chr, , drop = FALSE]
  d <- d[d$end >= r$start & d$start <= r$end, , drop = FALSE]
  if (!nrow(d)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = (r$start + r$end) / 2e6, y = 0, label = "No intervals in this region",
                        size = .txt(base_size, 0.75), colour = "grey40") +
      ggplot2::coord_cartesian(xlim = c(r$start, r$end) / 1e6, ylim = c(-1, 1), expand = FALSE) +
      ggplot2::labs(x = sprintf("Chromosome %s (Mb)", r$chr), y = title) +
      theme_quard(base_size) + .no_y()
    attr(p, "quard_panel") <- "track"
    return(p)
  }
  if (is.null(d$group)) d$group <- rep("x", nrow(d))
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$start / 1e6, xmax = .data$end / 1e6, ymin = -0.3,
                                    ymax = 0.3, fill = .data$group), colour = NA) +
    ggplot2::coord_cartesian(xlim = c(r$start, r$end) / 1e6, ylim = c(-1, 1), expand = FALSE) +
    ggplot2::labs(x = sprintf("Chromosome %s (Mb)", r$chr), y = title) +
    theme_quard(base_size) + .no_y()
  p <- if (length(unique(d$group)) == 1L) {
    p + ggplot2::scale_fill_manual(values = stats::setNames(fill, unique(d$group)), guide = "none")
  } else {
    p + ggplot2::labs(fill = NULL)
  }
  if (!is.null(d$label) && nrow(d)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(x = (pmax(.data$start, r$start) + pmin(.data$end, r$end)) / 2e6,
                                             y = -0.55, label = .data$label),
                                size = .txt(base_size, 0.75), vjust = 1)
  }
  attr(p, "quard_panel") <- "track"
  p
}

#' Stack panels with a shared genomic axis
#'
#' Combines any panels from this package (or your own ggplot objects that
#' use the same region limits) into one column, aligning the plotting
#' areas. The x axis is kept only on the lowest regional panel above the LD
#' heatmap.
#'
#' @param ... ggplot objects, or one list of them, from top to bottom.
#' @param heights Relative heights.
#' @param axis_panel Index of the panel that keeps its x axis (default:
#'   chosen automatically).
#' @return A patchwork object (print it, or save it with [save_pdf()]).
#' @export
combine_panels <- function(..., heights = NULL, axis_panel = NULL) {
  plots <- list(...)
  if (length(plots) == 1L && is.list(plots[[1]]) && !inherits(plots[[1]], "gg")) plots <- plots[[1]]
  plots <- Filter(Negate(is.null), plots)
  if (!length(plots)) .stop("No panels to combine.")
  kinds <- vapply(plots, function(p) attr(p, "quard_panel") %||% "other", "")
  regional <- which(kinds %in% c("regional", "genes", "track", "ld", "other"))
  if (is.null(axis_panel) && length(regional)) {
    non_ld <- setdiff(regional, which(kinds == "ld"))
    axis_panel <- if (length(non_ld)) max(non_ld) else max(regional)
  }
  for (k in regional) if (!identical(k, axis_panel)) plots[[k]] <- plots[[k]] + .no_x()
  patchwork::wrap_plots(plots, ncol = 1, heights = heights)
}

#' Save a figure as PDF
#'
#' @param plot A ggplot or patchwork object.
#' @param file Output file; must end in `.pdf`.
#' @param width,height Size. Default: the size recommended by the plot (for
#'   trees from [plot_tree()], large enough for all sample names), otherwise
#'   180 x 200 mm (BMC full-page width is 170-180 mm, maximum height 225 mm).
#' @param units `"mm"`, `"cm"` or `"in"`.
#' @return The file path (invisibly).
#' @export
save_pdf <- function(plot, file, width = NULL, height = NULL, units = c("mm", "cm", "in")) {
  units <- match.arg(units)
  if (!grepl("\\.pdf$", file, ignore.case = TRUE)) .stop("Only PDF output is supported: `file` must end in .pdf.")
  f <- c(mm = 1 / 25.4, cm = 1 / 2.54, `in` = 1)[[units]]
  auto <- attr(plot, "size_mm")
  w_in <- if (!is.null(width)) width * f else if (!is.null(auto)) auto[1] / 25.4 else 180 / 25.4
  h_in <- if (!is.null(height)) height * f else if (!is.null(auto)) auto[2] / 25.4 else 200 / 25.4
  dir <- dirname(file)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  grDevices::pdf(file, width = w_in, height = h_in, useDingbats = FALSE, onefile = FALSE)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)
  print(plot)
  invisible(file)
}
