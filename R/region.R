# Region and target definition ---------------------------------------------

#' Define the plotting region around a target variant
#'
#' By default the region is the LD block that contains the target variant:
#'
#' 1. `start`/`end`, if given, are used as they are;
#' 2. otherwise the block in `blocks` (from a block file or
#'    [find_ld_blocks()]) that contains the target;
#' 3. if the target is in no block and genotypes are given (`geno`), its
#'    LD interval: the span of all variants with \eqn{r^2 \ge} `ld_r2`
#'    with the target (the target's own LD block);
#' 4. without genotypes, a window of `flank` bp on each side.
#'
#' For LD blocks and LD intervals a small `margin` (default 2% of the width
#' on each side, at least 200 bp) is added so the first and last variants of
#' the block are not drawn on the panel border; explicit windows are used as
#' they are.
#' Set `min_width` to widen very small blocks, or pass `genes` to extend the
#' region to genes that cross its borders.
#'
#' @param chr,pos Target chromosome and position (bp, 1-based).
#' @param blocks Optional block table (`chr`, `start`, `end`).
#' @param geno Optional genotypes around the target (a `qp_geno` object,
#'   ideally after [genotype_qc()]) used for the LD interval when the target
#'   is in no block.
#' @param ld_r2 \eqn{r^2} threshold with the target for the LD interval.
#' @param flank Half-width (bp) of the window used without block or genotypes.
#' @param start,end Optional explicit region (overrides everything else).
#' @param min_width Minimum region width (bp); 0 keeps the block as it is.
#' @param margin Fraction of the width added on each side (at least 200 bp).
#' @param genes Optional `qp_genes` object; genes crossing the region
#'   borders then extend the region.
#' @param chrom_length Optional chromosome length used to clip the region.
#' @param verbose Print the chosen region.
#' @return A `qp_region` list: `chr`, `start`, `end` (plot range), `core`
#'   (the block or interval itself), `source`, `block`, `target`.
#' @export
define_region <- function(chr, pos, blocks = NULL, geno = NULL, ld_r2 = 0.6, flank = 50000,
                          start = NULL, end = NULL, min_width = 0, margin = 0.02, genes = NULL,
                          chrom_length = NULL, verbose = TRUE) {
  chr <- normalize_chr(chr)
  .check_number(pos, "pos", lower = 1)
  blk <- NULL
  src <- NULL
  s <- e <- NA_real_
  if (!is.null(start) && !is.null(end)) {
    s <- as.numeric(start)
    e <- as.numeric(end)
    src <- "user-defined"
  } else {
    if (!is.null(blocks) && nrow(blocks)) {
      hit <- which(blocks$chr == chr & blocks$start <= pos & blocks$end >= pos)
      if (length(hit)) {
        blk <- blocks[hit[1], , drop = FALSE]
        s <- blk$start
        e <- blk$end
        src <- sprintf("LD block containing the target%s",
                       if (!is.null(blk$nsnps) && !is.na(blk$nsnps)) paste0(", ", blk$nsnps, " SNPs") else "")
      }
    }
    if (is.na(s) && !is.null(geno) && inherits(geno, "qp_geno")) {
      gi <- geno$info
      on <- which(gi$chr == chr)
      if (length(on) >= 2L) {
        if (!is.null(blocks)) {
          .msg(verbose, "Target chr ", chr, ":", .fmt(pos), " is not inside any Gabriel LD block; ",
               "using its LD interval (variants with r2 >= ", ld_r2, " with the target).")
        }
        g1 <- .geno_subset(geno, on)
        li <- which.min(abs(g1$info$pos - pos))
        if (abs(g1$info$pos[li] - pos) > 0) {
          .msg(verbose, "The target position is not genotyped; LD is taken from the nearest variant ",
               g1$info$id[li], " (", .fmt(g1$info$pos[li]), ").")
        }
        r2 <- ld_with_lead(g1, li)
        sel <- which(!is.na(r2) & r2 >= ld_r2)
        s <- min(g1$info$pos[sel])
        e <- max(g1$info$pos[sel])
        if (s == e) {
          s <- if (li > 1L) g1$info$pos[li - 1L] else s - 1000
          e <- if (li < nrow(g1$info)) g1$info$pos[li + 1L] else e + 1000
          src <- sprintf("no variant in LD (r2 >= %s) with the target; range between its neighbouring variants", ld_r2)
        } else {
          src <- sprintf("LD interval of the target (%d variants with r2 >= %s)", length(sel), ld_r2)
        }
      }
    }
    if (is.na(s)) {
      if (!is.null(blocks)) {
        .msg(verbose, "Target chr ", chr, ":", .fmt(pos), " is not inside any LD block and no genotypes ",
             "were given; using +/-", .fmt(flank), " bp around it.")
      }
      s <- pos - flank
      e <- pos + flank
      src <- "window"
    }
  }
  if (e < s) .stop("Region end must be larger than its start.")
  if (pos < s || pos > e) {
    .warn("The target position (", .fmt(pos), ") lies outside the region ", .fmt(s), "-", .fmt(e), ".")
  }
  core <- c(s, e)
  if (e - s < min_width) {
    mid <- (s + e) / 2
    s <- mid - min_width / 2
    e <- mid + min_width / 2
  }
  if (!is.null(genes) && inherits(genes, "qp_genes") && nrow(genes$genes)) {
    gg <- genes$genes[genes$genes$chr == chr & genes$genes$start <= e & genes$genes$end >= s, , drop = FALSE]
    if (nrow(gg)) {
      s <- min(s, gg$start)
      e <- max(e, gg$end)
    }
  }
  ld_based <- !(src %in% c("user-defined", "window"))
  pad <- if (ld_based && margin > 0) max(margin * (e - s), 200) else 0
  s <- max(1, floor(s - pad))
  e <- ceiling(e + pad)
  if (!is.null(chrom_length) && is.finite(chrom_length)) {
    if (pos > chrom_length) {
      .warn("Target position ", .fmt(pos), " is beyond the length of chromosome ", chr, " (",
            .fmt(chrom_length), " bp): check that all files use the same genome assembly.")
    }
    e <- min(e, chrom_length)
  }
  r <- structure(list(chr = chr, start = s, end = e, core = core, source = src, block = blk, target = pos),
                 class = "qp_region")
  .msg(verbose, sprintf("Region: chr %s:%s-%s (%s kb; %s).", chr, .fmt(core[1]), .fmt(core[2]),
                        format(round((core[2] - core[1]) / 1000, 1), nsmall = 1), src))
  r
}

# Resolve the target variant from `target` (SNP ID, "chr:pos" or "top")
# and/or chr + pos.
.resolve_target <- function(gwas, target = NULL, chr = NULL, pos = NULL, verbose = TRUE) {
  if (!is.null(target)) {
    target <- as.character(target)
    if (tolower(target) == "top") {
      d <- if (!is.null(chr)) gwas[gwas$chr == normalize_chr(chr), , drop = FALSE] else gwas
      if (!nrow(d)) .stop("No GWAS SNPs on chromosome '", chr, "'.")
      k <- which.min(d$p)
      .msg(verbose, "Target: most significant SNP ", d$snp[k], " (chr ", d$chr[k], ":", .fmt(d$pos[k]),
           ", P = ", format(d$p[k], digits = 3), ").")
      return(list(chr = d$chr[k], pos = d$pos[k], snp = d$snp[k]))
    }
    k <- match(target, gwas$snp)
    if (!is.na(k)) return(list(chr = gwas$chr[k], pos = gwas$pos[k], snp = gwas$snp[k]))
    if (grepl("^[^:]+:[0-9]+$", target)) {
      chr <- sub(":.*$", "", target)
      pos <- as.numeric(sub("^.*:", "", target))
    } else {
      near <- utils::head(gwas$snp[agrep(target, gwas$snp, max.distance = 0.1)], 5)
      .stop("Target SNP '", target, "' was not found in the GWAS file.",
            if (length(near)) paste0(" Similar IDs: ", paste(near, collapse = ", "), ".") else "")
    }
  }
  if (is.null(chr) || is.null(pos)) {
    .stop("Give the target as `target` (SNP ID, 'chr:pos' or 'top') or as `chr` and `pos`.")
  }
  chr <- normalize_chr(chr)
  pos <- as.numeric(pos)
  .check_number(pos, "pos", lower = 1)
  if (!chr %in% gwas$chr) {
    .stop("Chromosome '", chr, "' is not in the GWAS data. GWAS chromosomes: ",
          paste(utils::head(sort_chr(gwas$chr), 30), collapse = ", "), ".")
  }
  k <- which(gwas$chr == chr & gwas$pos == pos)
  if (length(k)) return(list(chr = chr, pos = pos, snp = gwas$snp[k[1]]))
  .msg(verbose, "No GWAS SNP at exactly chr ", chr, ":", .fmt(pos),
       "; the most significant SNP in the region will be used as the lead variant.")
  list(chr = chr, pos = pos, snp = NA_character_)
}
