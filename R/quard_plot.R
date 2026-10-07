# Main function -------------------------------------------------------------------

#' Locus figure (PDF) and tables (CSV or Excel)
#'
#' The main function. From your GWAS results, gene annotation (GFF) and
#' genotypes (VCF or HapMap) it finds the region around the SNP you choose,
#' calculates linkage disequilibrium (LD) and writes:
#'
#' * **the figure** (PDF): the panels you choose with `panels`, all on the
#'   same genomic axis - genome-wide Manhattan plot, regional plot coloured
#'   by LD (\eqn{r^2}) with the lead SNP, gene models and the LD heatmap;
#' * **the tables** you choose with `tables`, as CSV files (open in Excel)
#'   or as one Excel workbook (`table_format = "xlsx"`):
#'   `"genes"` = genes in the region (GFF lines, `<output>_genes.txt`),
#'   `"snps"` = GWAS SNPs in the region with P value, \eqn{r^2} with the lead
#'   SNP, closest gene and where the SNP lies (exon, intron, promoter, ...),
#'   `"ld"` = LD between every pair of SNPs as a matrix, `"ld_pairs"` = the
#'   same as a list of pairs. A README explains every column.
#'
#' The tables do not depend on the panels: the LD table is written even when
#' the LD heatmap is not drawn, and `panels = "none"` writes only the tables.
#'
#' **Which region is shown?** With `region = "auto"` (default): `start`/`end`
#' if you give them; otherwise the LD block that contains the target SNP
#' (from your `blocks` file, or found in the genotypes with
#' [find_ld_blocks()], the same method as PLINK `--blocks`); if the SNP is in
#' no block, all SNPs with \eqn{r^2 \ge} `ld_r2` with it (its LD interval);
#' without genotypes, `flank` bp on each side. A 2% margin is added on each
#' side (`margin`).
#'
#' **LD.** Genotypes in the region are filtered with [genotype_qc()]
#' (two alleles, missing rate <= `max_missing`, minor allele frequency >=
#' `min_maf`; missing calls stay missing). \eqn{r^2} is calculated with
#' [calc_ld()] (identical to PLINK 1.9 `--r2`).
#'
#' All files must use the same reference genome version, and the genotypes
#' should be from the GWAS population. Choose the species first with
#' [set_species()] so that the Manhattan plot uses the right chromosome
#' lengths.
#'
#' @param gwas GWAS results file (CSV/TSV from GAPIT, TASSEL, PLINK, GEMMA,
#'   rMVP, ...; columns are recognised automatically) or [read_gwas()] output.
#' @param gff Gene annotation, GFF3 or GTF file (or [read_gff()] output).
#'   Optional; needed for the gene panel, gene list and gene locations.
#' @param vcf Genotype file: VCF (`.vcf`, `.vcf.gz`) or HapMap (`.hmp.txt`,
#'   tab- or comma-separated). Optional; needed for LD.
#' @param blocks Optional LD-block file (PLINK `.blocks.det` or a table with
#'   CHR, START, END; see [read_blocks()]).
#' @param chr,pos Chromosome (any style: `9`, `"Chr9"`, `"chr09"`) and
#'   position (bp) of the SNP of interest.
#' @param output Name of the PDF file, e.g. `"results/locus.pdf"`; the table
#'   files are named after it (`results/locus_SNPs.csv`, ...).
#' @param target Instead of `chr`/`pos`: a SNP name from the GWAS file,
#'   `"chr:pos"`, or `"top"` for the most significant SNP (on `chr`, if given).
#' @param panels What to draw: any of `"manhattan"`, `"regional"`, `"genes"`,
#'   `"ld"` (LD heatmap), or `"none"` for no figure (tables only).
#' @param tables Which tables to write: any of `"genes"`, `"snps"`, `"ld"`,
#'   `"ld_pairs"`, or `"none"`. Default: `c("genes", "snps", "ld")`.
#' @param table_format `"csv"` (one CSV file per table, plus a README text
#'   file) or `"xlsx"` (one Excel workbook `<output>_tables.xlsx` with one
#'   sheet per table and a README sheet; needs the `writexl` package).
#' @param region `"auto"`, `"blocks"`, `"ld_blocks"` or `"window"` (see
#'   "Which region is shown?").
#' @param start,end Your own region (bp); overrides `region`.
#' @param flank Half-width (bp) of the window used without genotypes, or with
#'   `region = "window"`.
#' @param ld_r2 \eqn{r^2} threshold of the LD interval used when the target
#'   is in no LD block (see [define_region()]).
#' @param min_width Minimum region width (bp); 0 (default) shows the block
#'   as it is.
#' @param margin Fraction of the region width added on each side.
#' @param extend_to_genes Widen the region so that genes crossing its
#'   borders are drawn completely.
#' @param chrom_sizes Chromosome lengths for the Manhattan plot. Default: the
#'   species chosen with [set_species()]; otherwise VCF `##contig` lengths,
#'   otherwise estimated from the GWAS positions. Also accepts anything
#'   [chrom_sizes()] accepts (`"rice"`, a `.fai` file, ...).
#' @param threshold Significance line: `"bonferroni"` (0.05 / number of SNPs,
#'   default), a P value such as `5e-8`, a -log10 value such as `7`, or
#'   `NULL` for none.
#' @param suggestive Optional second (suggestive) line, same format.
#' @param min_maf,max_missing Genotype filters used before LD is calculated.
#' @param ld_method,ld_stat LD definition, see [calc_ld()]. `ld_stat =
#'   "dprime"` draws and writes D' instead of \eqn{r^2}.
#' @param max_ld_snps Maximum number of SNPs drawn in the LD heatmap (more
#'   are thinned evenly; the LD table still has all SNPs).
#' @param max_table_snps Maximum number of SNPs in the LD table.
#' @param label_ld_snps SNP names under the heatmap: `"lead"`, `"none"` or
#'   `"all"`.
#' @param block_search Half-width (bp) of the genotype window searched for
#'   LD blocks.
#' @param max_block_kb Maximum SNP distance within a block (kb).
#' @param max_block_snps Maximum number of SNPs used for the block search
#'   (denser data are thinned evenly; keeps whole-genome sequencing data fast).
#' @param transcripts `"canonical"` (one model per gene) or `"all"`
#'   transcripts in the gene panel.
#' @param gene_features Gene list file: `"gene"` (gene lines only) or `"all"`
#'   (gene, mRNA, exon, CDS and UTR lines).
#' @param upstream_bp,downstream_bp Size of the promoter (upstream) and
#'   downstream windows used to describe where each SNP lies (default 3 kb
#'   and 1 kb).
#' @param width,height Figure size in mm (height from the panels when `NULL`).
#' @param base_size Font size (pt).
#' @param ld_file Optional LD matrix made by PLINK (`plink --r2 square` with
#'   `ld_bim`, or the `plink --r2` table) used instead of calculating LD.
#' @param ld_bim The `.bim` file belonging to a square `ld_file`.
#' @param chr_map Optional chromosome name translation (see [normalize_chr()]).
#' @param gene_file,ld_csv Older ways to choose the tables (still work):
#'   `gene_file = FALSE` drops the gene list; `ld_csv = "matrix"`, `"pairs"`,
#'   `"both"` or `"none"` chooses the LD tables.
#' @param verbose Print progress messages.
#' @param ... Column names of the GWAS file when they are not recognised
#'   automatically (`snp_col`, `chr_col`, `pos_col`, `p_col`; passed to
#'   [read_gwas()]). Old argument names of version 0.1 (`gwas_file`,
#'   `gff_file`, `vcf_file`, `block_file_path`, `target_snp_chr`,
#'   `target_snp_pos`, `output_name`) are also still accepted.
#' @return Invisibly, a `quard_result` list: the figure (`plot`), `region`,
#'   `lead` SNP, GWAS rows in the region (`gwas`), `genes`, `genotypes`, LD
#'   matrices (`ld` = heatmap, `ld_all` = table), `snp_table`, `tables`,
#'   `snp_annotation`, `r2_lead`, `blocks`, `threshold`, `panels` and the
#'   `files` written.
#' @examples
#' \donttest{
#' ex <- simulate_example_data(tempfile("qp"))
#' set_species("rice")
#' res <- quard_plot(ex$gwas, ex$gff, ex$vcf, target = "top",
#'                   output = file.path(tempdir(), "locus.pdf"))
#' res
#' # no LD heatmap in the figure, but the LD table as an Excel sheet:
#' if (requireNamespace("writexl", quietly = TRUE)) {
#'   quard_plot(ex$gwas, ex$gff, ex$vcf, target = "top",
#'              panels = c("manhattan", "regional", "genes"), table_format = "xlsx",
#'              output = file.path(tempdir(), "locus2.pdf"))
#' }
#' }
#' @export
quard_plot <- function(gwas, gff = NULL, vcf = NULL, blocks = NULL, chr = NULL, pos = NULL,
                       output = "quard_plot.pdf", target = NULL,
                       panels = c("manhattan", "regional", "genes", "ld"),
                       tables = c("genes", "snps", "ld"), table_format = c("csv", "xlsx"),
                       region = c("auto", "blocks", "ld_blocks", "window"),
                       start = NULL, end = NULL, flank = 50000, ld_r2 = 0.6, min_width = 0,
                       margin = 0.02, extend_to_genes = FALSE, chrom_sizes = NULL,
                       threshold = "bonferroni", suggestive = NULL,
                       min_maf = 0.05, max_missing = 0.2,
                       ld_method = c("genotype", "haplotype"), ld_stat = c("r2", "dprime"),
                       max_ld_snps = 300, max_table_snps = 5000, label_ld_snps = c("lead", "none", "all"),
                       block_search = 250000, max_block_kb = 200, max_block_snps = 1500,
                       transcripts = c("canonical", "all"), gene_features = c("gene", "all"),
                       upstream_bp = 3000, downstream_bp = 1000,
                       width = 180, height = NULL, base_size = 8,
                       ld_file = NULL, ld_bim = NULL, chr_map = NULL,
                       gene_file = NULL, ld_csv = NULL, verbose = TRUE, ...) {
  dots <- list(...)
  old <- c(gwas_file = "gwas", gff_file = "gff", vcf_file = "vcf", block_file_path = "blocks",
           target_snp_chr = "chr", target_snp_pos = "pos", output_name = "output")
  gw_cols <- c("snp_col", "chr_col", "pos_col", "p_col")
  gw_args <- dots[intersect(names(dots), gw_cols)]
  dots <- dots[setdiff(names(dots), gw_cols)]
  bad <- setdiff(names(dots), names(old))
  if (length(bad) || (length(dots) && is.null(names(dots)))) {
    .stop("Unknown argument(s): ", paste(bad, collapse = ", "), ".")
  }
  if (missing(gwas)) {
    if (is.null(dots$gwas_file)) .stop("`gwas` (the GWAS results file) is required.")
    gwas <- dots$gwas_file
  }
  for (nm in setdiff(names(dots), "gwas_file")) assign(old[[nm]], dots[[nm]])

  # 0. options -------------------------------------------------------------------
  panels <- .parse_choice(panels, c("manhattan", "regional", "genes", "ld"), "panel")
  tables <- .parse_choice(tables, c("genes", "snps", "ld", "ld_pairs"), "table")
  if (!is.null(gene_file)) tables <- if (isTRUE(gene_file)) union(tables, "genes") else setdiff(tables, "genes")
  if (!is.null(ld_csv)) {
    ld_csv <- match.arg(ld_csv, c("matrix", "pairs", "both", "none"))
    tables <- union(setdiff(tables, c("ld", "ld_pairs")),
                    switch(ld_csv, matrix = "ld", pairs = "ld_pairs", both = c("ld", "ld_pairs"),
                           none = character()))
  }
  tables <- intersect(c("genes", "snps", "ld", "ld_pairs"), tables)
  if (!length(panels) && !length(tables)) .stop("Nothing to do: `panels` and `tables` are both \"none\".")
  table_format <- match.arg(table_format)
  region <- match.arg(region)
  ld_method <- match.arg(ld_method)
  ld_stat <- match.arg(ld_stat)
  label_ld_snps <- match.arg(label_ld_snps)
  transcripts <- match.arg(transcripts)
  gene_features <- match.arg(gene_features)
  output <- .pdf_name(output)
  base <- sub("\\.pdf$", "", output, ignore.case = TRUE)
  if (table_format == "xlsx" && length(tables) && is.na(.xlsx_engine())) {
    .warn(.xlsx_help, " CSV files are written instead.")
    table_format <- "csv"
  }
  want_ld_table <- any(c("ld", "ld_pairs") %in% tables)

  # 1. GWAS and target -----------------------------------------------------
  gw <- if (inherits(gwas, "qp_gwas")) gwas else
    do.call(read_gwas, c(list(gwas, chr_map = chr_map, verbose = verbose), gw_args))
  tgt <- .resolve_target(gw, target, chr, pos, verbose = verbose)

  # 2. chromosome lengths ----------------------------------------------------
  cs <- chrom_sizes %||% .species_sizes()
  if (is.null(cs) && is.character(vcf) && length(vcf) == 1L && grepl("\\.vcf", tolower(vcf))) {
    cs_vcf <- tryCatch(chrom_sizes(vcf, chr_map = chr_map), error = function(e) NULL)
    if (!is.null(cs_vcf) && all(unique(gw$chr) %in% cs_vcf$chr)) cs <- cs_vcf
  }
  sizes <- if (is.data.frame(cs) && !is.null(attr(cs, "source"))) cs else chrom_sizes(cs, gwas = gw, chr_map = chr_map)
  .msg(verbose, "Chromosome lengths: ", attr(sizes, "source"), ".")
  chr_len <- sizes$length[match(tgt$chr, sizes$chr)]

  # 3. genes -------------------------------------------------------------------
  genes <- NULL
  if (!is.null(gff)) {
    genes <- if (inherits(gff, "qp_genes")) gff else read_gff(gff, chr = tgt$chr, chr_map = chr_map,
                                                               verbose = verbose)
  }
  if ("genes" %in% panels && is.null(genes)) {
    .msg(verbose, "No annotation (`gff`) given: the gene panel is skipped.")
    panels <- setdiff(panels, "genes")
  }

  # 4. region --------------------------------------------------------------------
  mode <- region
  if (!is.null(start) && !is.null(end)) {
    mode <- "user"
  } else if (mode == "auto") {
    mode <- if (!is.null(blocks)) "blocks" else if (!is.null(vcf)) "ld_blocks" else "window"
  }
  blk <- NULL
  geno_win <- NULL
  gq <- NULL
  win <- c(-Inf, Inf)
  if (mode == "blocks" && is.null(blocks)) .stop("region = 'blocks' needs an LD-block file (`blocks`).")
  if (mode == "ld_blocks" && is.null(vcf)) .stop("region = 'ld_blocks' needs genotypes (`vcf`).")
  if (mode %in% c("blocks", "ld_blocks") && !is.null(vcf)) {
    # genotypes around the target: used to find the LD block (or the LD
    # interval of the target when it lies in no block)
    win <- c(max(1, tgt$pos - block_search), tgt$pos + block_search)
    geno_win <- read_genotypes(vcf, chr = tgt$chr, start = win[1], end = win[2], chr_map = chr_map,
                               verbose = verbose)
    gq <- genotype_qc(geno_win, min_maf = min_maf, max_missing = max_missing, verbose = FALSE)
    if (nrow(gq$info) > max_block_snps) {
      nv <- nrow(gq$info)
      near <- which.min(abs(gq$info$pos - tgt$pos))
      sel <- sort(unique(c(round(seq(1, nv, length.out = max_block_snps)), near)))
      gq <- .geno_subset(gq, sel)
      .msg(verbose, "LD-block search uses ", .fmt(length(sel)), " evenly spaced of ", .fmt(nv),
           " variants (max_block_snps = ", .fmt(max_block_snps), ").")
    }
    if (nrow(gq$info) < 2L) gq <- NULL
  }
  if (mode == "blocks") {
    blk <- read_blocks(blocks, chr_map = chr_map)
  } else if (mode == "ld_blocks" && !is.null(gq)) {
    blk <- tryCatch(find_ld_blocks(gq, max_kb = max_block_kb, verbose = verbose),
                    error = function(e) {
                      .warn("LD-block detection failed (", conditionMessage(e), ").")
                      NULL
                    })
    if (is.null(blk)) blk <- .empty_blocks()
  }
  reg <- define_region(tgt$chr, tgt$pos, blocks = blk, geno = gq, ld_r2 = ld_r2, flank = flank,
                       start = start, end = end, min_width = min_width, margin = margin,
                       genes = if (extend_to_genes) genes else NULL, chrom_length = chr_len,
                       verbose = verbose)

  in_reg <- gw$chr == reg$chr & gw$pos >= reg$start & gw$pos <= reg$end
  if (!any(in_reg) && "regional" %in% panels) {
    .stop("No GWAS SNPs in chr ", reg$chr, ":", .fmt(reg$start), "-", .fmt(reg$end),
          ". Check the target position and that all files use the same genome assembly.")
  }
  gwr <- gw[in_reg, , drop = FALSE]
  lead_snp <- if (!is.na(tgt$snp) && tgt$snp %in% gwr$snp) tgt$snp else if (nrow(gwr)) gwr$snp[which.min(gwr$p)] else NA_character_

  # 5. genotypes and LD ----------------------------------------------------------
  geno <- NULL
  ldmat <- NULL
  ld_full <- NULL
  r2lead <- NULL
  lead_ld_id <- NULL
  need_geno <- any(c("regional", "ld") %in% panels) || want_ld_table || "snps" %in% tables
  if (!is.null(vcf) && need_geno) {
    geno <- if (!is.null(geno_win) && reg$start >= win[1] && reg$end <= win[2]) {
      .geno_subset(geno_win, which(geno_win$info$pos >= reg$start & geno_win$info$pos <= reg$end))
    } else {
      read_genotypes(vcf, chr = reg$chr, start = reg$start, end = reg$end, chr_map = chr_map, verbose = verbose)
    }
    geno <- genotype_qc(geno, min_maf = min_maf, max_missing = max_missing, verbose = verbose)
    if (nrow(geno$info) < 2L) {
      .warn("Fewer than two variants passed genotype QC in the region: no LD colouring, LD heatmap or LD table.")
      geno <- NULL
    }
  }
  if (!is.null(geno)) {
    geno$info$gwas_snp <- gwr$snp[match(geno$info$pos, gwr$pos)]
    li <- if (!is.na(lead_snp)) match(lead_snp, geno$info$gwas_snp) else NA_integer_
    if (is.na(li)) {
      cand <- which(!is.na(geno$info$gwas_snp))
      if (length(cand)) {
        li <- cand[which.min(gwr$p[match(geno$info$gwas_snp[cand], gwr$snp)])]
      } else {
        li <- which.min(abs(geno$info$pos - tgt$pos))
      }
      # The lead SNP stays the labelled SNP of the figure; only the LD colours,
      # the LD heatmap mark and the r2 column use the genotyped proxy.
      .msg(verbose, "Lead SNP ", lead_snp, " is not among the genotyped variants that passed QC ",
           "(not in the genotype file, or removed by the MAF/missing-call filters); r2 colours are ",
           "shown relative to ", ifelse(is.na(geno$info$gwas_snp[li]), geno$info$id[li], geno$info$gwas_snp[li]),
           " instead.")
    }
    r2v <- ld_with_lead(geno, li)
    has <- !is.na(geno$info$gwas_snp)
    r2lead <- stats::setNames(unname(r2v[has]), geno$info$gwas_snp[has])
    if ("ld" %in% panels || want_ld_table) {
      g_ld <- geno
      g_ld$info$id <- ifelse(is.na(geno$info$gwas_snp), geno$info$id, geno$info$gwas_snp)
      n_v <- nrow(g_ld$info)
      lead_ld_id <- g_ld$info$id[li]
      if ("ld" %in% panels) {
        if (n_v > max_ld_snps) {
          sel <- sort(unique(c(round(seq(1, n_v, length.out = max_ld_snps - 1L)), li)))
          .msg(verbose, "LD heatmap: ", n_v, " variants thinned to ", length(sel),
               " evenly spaced variants (max_ld_snps = ", max_ld_snps, ").")
          ldmat <- calc_ld(.geno_subset(g_ld, sel), stat = ld_stat, method = ld_method)
        } else {
          ldmat <- calc_ld(g_ld, stat = ld_stat, method = ld_method)
        }
      }
      if (want_ld_table) {
        if (!is.null(ldmat) && n_v <= max_ld_snps) {
          ld_full <- ldmat
        } else {
          g_full <- g_ld
          if (n_v > max_table_snps) {
            sel5 <- sort(unique(c(round(seq(1, n_v, length.out = max_table_snps)), li)))
            g_full <- .geno_subset(g_full, sel5)
            .msg(verbose, "LD table limited to ", .fmt(length(sel5)), " evenly spaced of ", .fmt(n_v),
                 " variants (max_table_snps = ", .fmt(max_table_snps), ").")
          }
          ld_full <- calc_ld(g_full, stat = ld_stat, method = ld_method)
        }
      }
    }
  }

  # 5b. LD from a PLINK file ------------------------------------------------------
  if (!is.null(ld_file)) {
    lm <- read_plink_ld(ld_file, bim = ld_bim)
    lp <- attr(lm, "pos")
    if (is.null(lp)) {
      .stop("`ld_file` needs variant positions: use the PLINK --r2 table output, or give `ld_bim`.")
    }
    inr <- which(lp >= reg$start & lp <= reg$end)
    if (length(inr) < 2L) {
      .warn("Fewer than two variants of `ld_file` lie in the region; it is not used.")
    } else {
      pos_l <- lp[inr]
      ids_l <- rownames(lm)[inr]
      gid <- gwr$snp[match(pos_l, gwr$pos)]
      ids_l <- ifelse(is.na(gid), ids_l, gid)
      L <- lm[inr, inr, drop = FALSE]
      dimnames(L) <- list(ids_l, ids_l)
      lead_ld_id <- if (!is.na(lead_snp) && lead_snp %in% ids_l) lead_snp else ids_l[which.min(abs(pos_l - tgt$pos))]
      if (is.null(r2lead)) {
        r2lead <- stats::setNames(unname(L[lead_ld_id, ]), ids_l)
        r2lead <- r2lead[!is.na(gid)]
      }
      li2 <- match(lead_ld_id, ids_l)
      ld_full <- L
      attr(ld_full, "pos") <- pos_l
      attr(ld_full, "chr") <- reg$chr
      attr(ld_full, "stat") <- "r2"
      attr(ld_full, "method") <- "PLINK file"
      if (length(inr) > max_ld_snps) {
        sel <- sort(unique(c(round(seq(1, length(inr), length.out = max_ld_snps - 1L)), li2)))
        L <- L[sel, sel, drop = FALSE]
        pos_l <- pos_l[sel]
      }
      attr(L, "pos") <- pos_l
      attr(L, "chr") <- reg$chr
      attr(L, "stat") <- "r2"
      attr(L, "method") <- "PLINK file"
      ldmat <- L
      .msg(verbose, "LD heatmap: ", nrow(L), " variants taken from ", basename(ld_file), ".")
    }
  }

  # 6. figure -------------------------------------------------------------------
  th <- gwas_threshold(gw, threshold)
  if (!is.null(th)) .msg(verbose, "Significance line: ", th$label, " (-log10 P = ", round(th$logp, 2), ").")
  plist <- list()
  hts <- numeric()
  if ("manhattan" %in% panels) {
    pm <- plot_manhattan(gw, chrom_sizes = sizes, threshold = th, suggestive = suggestive,
                         highlight = reg, target = lead_snp, base_size = base_size)
    plist$manhattan <- pm
    hts <- c(hts, 42)
    if (any(c("regional", "genes", "ld") %in% panels)) {
      plist$connector <- plot_zoom_connector(pm, reg)
      hts <- c(hts, 7)
    }
  }
  if ("regional" %in% panels) {
    plist$regional <- plot_regional(gw, reg, lead = lead_snp, r2 = r2lead, threshold = th,
                                    base_size = base_size)
    hts <- c(hts, 48)
  }
  if ("genes" %in% panels) {
    pg <- plot_genes(genes, reg, transcripts = transcripts, base_size = base_size)
    plist$genes <- pg
    hts <- c(hts, min(60, 9 + 6.5 * attr(pg, "quard_rows")))
  }
  if ("ld" %in% panels) {
    if (is.null(ldmat)) {
      .msg(verbose, "No genotypes for the region: the LD heatmap is skipped.")
    } else {
      blk_reg <- if (!is.null(blk) && nrow(blk)) blk[blk$chr == reg$chr & blk$end >= reg$start & blk$start <= reg$end, , drop = FALSE] else NULL
      plist$ld <- plot_ld(ldmat, reg, lead = lead_ld_id, max_snps = max_ld_snps,
                          label_snps = label_ld_snps, blocks = blk_reg, base_size = base_size)
      hts <- c(hts, 72)
    }
  }
  files <- character()
  desc <- character()
  fig <- NULL
  if (!length(panels)) {
    .msg(verbose, "No figure (panels = \"none\"): only the tables are written.")
  } else if (!length(plist)) {
    if (!length(tables)) .stop("Nothing to plot: check `panels` and the inputs.")
    .msg(verbose, "None of the chosen panels can be drawn with these inputs: only the tables are written.")
  } else {
    fig <- combine_panels(plist, heights = hts)
    if (is.null(height)) height <- min(sum(hts) + 8, 260)
    save_pdf(fig, output, width = width, height = height)
    files["pdf"] <- output
    shown <- intersect(c("manhattan", "regional", "genes", "ld"), names(plist))
    desc["pdf"] <- sprintf("the figure (%s)", paste(c(manhattan = "Manhattan plot", regional = "regional plot",
                                                       genes = "genes", ld = "LD heatmap")[shown], collapse = ", "))
  }

  # 7. tables ------------------------------------------------------------------------
  genes_reg <- if (!is.null(genes)) genes_in_region(genes, reg) else NULL
  lead_pos <- if (!is.na(lead_snp) && lead_snp %in% gwr$snp) gwr$pos[match(lead_snp, gwr$snp)] else tgt$pos
  tabs <- list()
  snp_tab <- NULL
  gene_tab <- NULL
  if ("genes" %in% tables) {
    if (is.null(genes)) {
      .msg(verbose, "No GFF file (`gff`) given: no gene list.")
    } else {
      gene_path <- paste0(base, "_genes.txt")
      write_region_genes(genes, gene_path, region = reg, features = gene_features)
      files["genes"] <- gene_path
      desc["genes"] <- if (nrow(genes_reg$genes)) {
        sprintf("%d gene(s) in the region, as lines of the GFF file", nrow(genes_reg$genes))
      } else {
        "no genes in the region (the file only has a header)"
      }
      gene_tab <- .gene_table(genes_reg, lead_pos)
    }
  }
  if ("snps" %in% tables) {
    if (!nrow(gwr)) {
      .msg(verbose, "No GWAS SNPs in the region: no SNP table.")
    } else {
      snp_tab <- .snp_table(gwr, reg, lead_snp, r2lead, th, genes, upstream_bp, downstream_bp)
      tabs$SNPs <- snp_tab
    }
  }
  ld_tag <- if (!is.null(ld_full) && identical(attr(ld_full, "stat"), "dprime")) "Dprime" else "r2"
  if (want_ld_table && is.null(ld_full)) {
    .msg(verbose, if (is.null(vcf) && is.null(ld_file)) "No genotypes (`vcf`) given: the LD table needs genotypes and is skipped." else
      "Not enough genotyped SNPs in the region: no LD table.")
  }
  if ("ld" %in% tables && !is.null(ld_full)) {
    tabs[[paste0("LD_", ld_tag)]] <- .ld_table(ld_full, "matrix", genes = genes, upstream_bp = upstream_bp,
                                               downstream_bp = downstream_bp)
  }
  if ("ld_pairs" %in% tables && !is.null(ld_full)) {
    tabs[[paste0("LD_", ld_tag, "_pairs")]] <- .ld_table(ld_full, "pairs", genes = genes,
                                                         upstream_bp = upstream_bp, downstream_bp = downstream_bp)
  }
  stat_lab <- if (ld_tag == "Dprime") "D'" else "r2"
  tab_desc <- function(nm) {
    if (nm == "SNPs") {
      parts <- c("P value", if ("R2_WITH_LEAD" %in% names(snp_tab)) "r2 with the lead SNP",
                 if ("GENE" %in% names(snp_tab)) "closest gene and where each SNP lies")
      sprintf("%d GWAS SNPs in the region: %s", nrow(snp_tab),
              paste(parts, collapse = ", "))
    } else if (nm == "Genes") {
      sprintf("%d gene(s) in the region: ID, position, strand, distance to the lead SNP, description", nrow(gene_tab))
    } else if (grepl("_pairs$", nm)) {
      sprintf("the same %s values as a list of SNP pairs, with their distance", stat_lab)
    } else {
      sprintf("%s between every pair of the %d genotyped SNPs (matrix), with the gene location of each SNP",
              stat_lab, nrow(ld_full))
    }
  }
  key_of <- function(nm) if (nm == "SNPs") "snps" else if (grepl("_pairs$", nm)) "ld_pairs" else "ld_matrix"
  info <- list(reg = reg, lead = lead_snp, ld_ref = lead_ld_id, lead_p = if (!is.na(lead_snp)) gwr$p[match(lead_snp, gwr$snp)] else NA,
               th = th, up = upstream_bp, down = downstream_bp, ld_tag = ld_tag,
               settings = sprintf(paste("region = %s; min_maf = %s; max_missing = %s; ld_method = %s;",
                                        "ld_stat = %s; ld_r2 = %s; max_ld_snps = %s"),
                                  mode, min_maf, max_missing, ld_method, ld_stat, ld_r2, max_ld_snps))
  if (length(tabs) && table_format == "csv") {
    tab_files <- character()
    for (nm in names(tabs)) {
      f <- paste0(base, "_", nm, ".csv")
      write_table(tabs[[nm]], f)
      files[key_of(nm)] <- f
      desc[key_of(nm)] <- tab_desc(nm)
      tab_files[nm] <- f
    }
    f <- paste0(base, "_README.txt")
    files["readme"] <- f
    desc["readme"] <- "what every file and column means"
    writeLines(.readme_text(info, tabs, tab_files, files, desc), f)
  } else if ((length(tabs) || !is.null(gene_tab)) && table_format == "xlsx") {
    sheets <- c(tabs, list(Genes = gene_tab))
    sheets <- sheets[!vapply(sheets, is.null, NA)]
    for (nm in names(sheets)[vapply(sheets, nrow, 0L) > 1048575L]) {
      f <- paste0(base, "_", nm, ".csv")
      write_table(sheets[[nm]], f)
      files[key_of(nm)] <- f
      desc[key_of(nm)] <- paste0(tab_desc(nm), " (too many rows for Excel, so written as CSV)")
      sheets[[nm]] <- NULL
    }
    xl <- paste0(base, "_tables.xlsx")
    files["xlsx"] <- xl
    desc["xlsx"] <- sprintf("Excel workbook with the sheets %s", paste(c("README", names(sheets)), collapse = ", "))
    sheet_desc <- vapply(names(sheets), tab_desc, "")
    readme <- .readme_sheet(info, sheets, sheet_desc, files, desc)
    write_table(c(list(README = readme), sheets), xl)
  }

  snp_ann <- NULL
  if (!is.null(genes) && !is.null(ld_full)) {
    snp_ann <- annotate_variants(ld_full, genes, upstream_bp = upstream_bp, downstream_bp = downstream_bp)
  }
  if (verbose && length(files)) {
    message("\nFiles written:\n", paste0("  ", format(files), "  ", desc[names(files)], collapse = "\n"))
  }
  res <- structure(list(plot = fig, region = reg, target = tgt, lead = lead_snp, gwas = gwr,
                        genes = genes_reg, genotypes = geno, ld = ldmat, ld_all = ld_full,
                        snp_table = snp_tab, tables = tabs, snp_annotation = snp_ann,
                        r2_lead = r2lead, blocks = blk, threshold = th, panels = plist,
                        files = files, file_descriptions = desc[names(files)]),
                   class = "quard_result")
  invisible(res)
}

# Options given as a vector or a comma-separated string; "none" = nothing.
.parse_choice <- function(x, choices, what) {
  if (is.null(x) || !length(x)) return(character())
  x <- tolower(trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE))))
  x <- x[nzchar(x)]
  if (!length(x) || all(x == "none")) return(character())
  hit <- choices[pmatch(x, choices, duplicates.ok = TRUE)]
  bad <- x[is.na(hit)]
  if (length(bad)) {
    .stop("Unknown ", what, "(s): ", paste0("\"", bad, "\"", collapse = ", "), ". Choose from ",
          paste0("\"", choices, "\"", collapse = ", "), ", or \"none\".")
  }
  intersect(choices, hit)
}

# The figure is always a PDF: add .pdf when there is no extension.
.pdf_name <- function(output) {
  if (!is.character(output) || length(output) != 1L || is.na(output) || !nzchar(output)) {
    .stop("`output` must be a file name, e.g. \"my_locus.pdf\".")
  }
  ext <- .file_ext(output)
  if (ext %in% c("png", "jpg", "jpeg", "tif", "tiff", "svg", "eps", "ps", "bmp", "gif", "html", "htm",
                 "csv", "xlsx", "xls", "txt")) {
    .stop("Only PDF output is supported: `output` must be a file name ending in .pdf ",
          "(e.g. \"my_locus.pdf\"); the table files are named after it.")
  }
  if (ext != "pdf") output <- paste0(output, ".pdf")
  output
}

# GWAS SNPs of the region with LD to the lead SNP and their gene location.
.snp_table <- function(gwr, reg, lead, r2lead, th, genes, upstream_bp, downstream_bp) {
  d <- gwr[order(gwr$pos), , drop = FALSE]
  out <- data.frame(SNP = d$snp, CHR = d$chr, POS = d$pos, P = d$p,
                    MINUS_LOG10_P = round(-log10(d$p), 3), stringsAsFactors = FALSE)
  if (!is.null(th)) out$SIGNIFICANT <- ifelse(d$p <= th$p, "yes", "no")
  out$IS_LEAD <- ifelse(!is.na(lead) & d$snp == lead, "yes", "no")
  if (!is.null(r2lead)) out$R2_WITH_LEAD <- round(unname(r2lead[d$snp]), 4)
  if (!is.null(reg$source) && grepl("LD", reg$source)) {
    core <- reg$core %||% c(reg$start, reg$end)
    out$IN_LD_BLOCK <- ifelse(d$pos >= core[1] & d$pos <= core[2], "yes", "no")
  }
  if (!is.null(genes)) {
    a <- annotate_variants(data.frame(snp = d$snp, chr = d$chr, pos = d$pos), genes,
                           upstream_bp = upstream_bp, downstream_bp = downstream_bp)
    out$GENE <- a$gene_id
    out$LOCATION <- a$location
    out$FEATURE <- a$feature
    out$DISTANCE_BP <- a$distance_bp
    out$GENE_STRAND <- a$strand
    out$GENE_DESCRIPTION <- a$description
    out$OTHER_GENES <- a$other_genes
  }
  rownames(out) <- NULL
  out
}

# Readable table of the genes in the region.
.gene_table <- function(genes_reg, lead_pos) {
  G <- genes_reg$genes
  if (is.null(G) || !nrow(G)) return(NULL)
  data.frame(GENE_ID = G$gene_id, NAME = G$name, CHR = G$chr, START = G$start, END = G$end,
             STRAND = G$strand, LENGTH_BP = G$end - G$start + 1,
             DISTANCE_TO_LEAD_BP = pmax(G$start - lead_pos, lead_pos - G$end, 0),
             DESCRIPTION = G$note, stringsAsFactors = FALSE)
}

# Plain-language meaning of every table column.
.column_help <- function(up, down) {
  c(SNP = "SNP (marker) name, as in your GWAS file",
    CHR = "Chromosome",
    POS = "Position on the chromosome (bp)",
    P = "GWAS P value",
    MINUS_LOG10_P = "-log10(P): the height of the SNP in the plots (higher = stronger association)",
    SIGNIFICANT = "yes = the P value passes the significance line",
    IS_LEAD = "yes = the lead SNP (all r2 values are calculated with this SNP)",
    R2_WITH_LEAD = paste("LD (r2) between this SNP and the lead SNP: 1 = always inherited together,",
                         "0 = independent. Empty = the SNP is not in the genotype file or was removed by the",
                         "genotype filters (minor allele frequency, missing data)"),
    IN_LD_BLOCK = "yes = inside the LD block of the target SNP (the plotted region without its margins)",
    GENE = "Closest gene (gene ID from the GFF file)",
    LOCATION = sprintf(paste("Where the SNP lies relative to that gene: exon; intron; upstream = promoter",
                             "region, up to %s bp before the start of the gene; downstream = up to %s bp after",
                             "the end of the gene; intergenic = further away from any gene"),
                       .fmt(up), .fmt(down)),
    FEATURE = "For exon SNPs, the part of the exon: CDS (protein-coding), 5'UTR, 3'UTR, UTR or non-coding exon",
    DISTANCE_BP = paste("Distance to the gene (bp): 0 inside the gene; upstream SNPs: to the gene start;",
                        "downstream SNPs: to the gene end; intergenic SNPs: to the nearest end of the gene"),
    GENE_STRAND = "Strand of the gene (+ or -): the direction in which the gene is read",
    GENE_DESCRIPTION = "Gene description from the GFF file",
    OTHER_GENES = "Other genes that also contain the SNP or have it in their promoter / downstream window",
    SNP_A = "First SNP of the pair", POS_A = "Position of the first SNP (bp)",
    SNP_B = "Second SNP of the pair", POS_B = "Position of the second SNP (bp)",
    DIST_BP = "Distance between the two SNPs (bp)",
    R2 = "LD r2 between the two SNPs (0 = independent, 1 = always inherited together)",
    Dprime = "LD D' between the two SNPs (0 to 1)",
    GENE_A = "Closest gene of the first SNP", LOCATION_A = "Location of the first SNP relative to that gene",
    GENE_B = "Closest gene of the second SNP", LOCATION_B = "Location of the second SNP relative to that gene",
    GENE_ID = "Gene ID from the GFF file", NAME = "Gene name", START = "Gene start (bp)", END = "Gene end (bp)",
    STRAND = "Strand (+ or -)", LENGTH_BP = "Gene length (bp)",
    DISTANCE_TO_LEAD_BP = "Distance from the lead SNP to the gene (0 = the lead SNP is inside the gene)",
    DESCRIPTION = "Gene description from the GFF file")
}

# General information rows of the README.
.readme_head <- function(info) {
  reg <- info$reg
  core <- reg$core %||% c(reg$start, reg$end)
  sp <- get_species()
  rows <- list(
    c("Made with", sprintf("quard_plot %s (R %s) on %s", .qp_version, getRversion(), format(Sys.Date()))),
    if (!is.null(info$settings)) c("Settings", info$settings),
    if (!is.null(sp)) c("Species", sprintf("%s, reference genome %s", sp$organism, sp$genome)),
    c("Region", sprintf("chromosome %s: %s - %s bp (%s); plotted with margins: %s - %s bp", reg$chr,
                        .fmt(round(core[1])), .fmt(round(core[2])), reg$source %||% "user-defined",
                        .fmt(round(reg$start)), .fmt(round(reg$end)))),
    c("Lead SNP", if (is.na(info$lead)) "none" else sprintf("%s (P = %s)", info$lead, format(info$lead_p, digits = 3))),
    if (!is.null(info$ld_ref) && !is.na(info$lead) && !identical(info$ld_ref, info$lead))
      c("LD reference", sprintf("%s (the lead SNP has no genotypes after QC, so r2 values refer to this SNP)", info$ld_ref)),
    if (!is.null(info$th)) c("Significance line", sprintf("%s (-log10 P = %s)", info$th$label, round(info$th$logp, 2))),
    c("Gene windows", sprintf("promoter (upstream) = %s bp before the gene start; downstream = %s bp after the gene end",
                              .fmt(info$up), .fmt(info$down))))
  do.call(rbind, rows[!vapply(rows, is.null, NA)])
}

# Column explanation rows of one table.
.readme_columns <- function(nm, d, info) {
  help <- .column_help(info$up, info$down)
  cols <- names(d)
  if (grepl("^LD_", nm) && !grepl("_pairs$", nm)) {
    cols <- intersect(cols, names(help))
    extra <- sprintf(paste("then one column per SNP: each cell is the %s between the SNP of the row and the SNP of",
                           "the column (1 on the diagonal)"), if (info$ld_tag == "Dprime") "D'" else "r2")
    return(rbind(cbind(cols, unname(help[cols])), c("(SNP columns)", extra)))
  }
  cbind(cols, vapply(cols, function(cc) help[[cc]] %||% "", ""))
}

.readme_text <- function(info, tabs, tab_files, files, desc) {
  h <- .readme_head(info)
  out <- c("quard_plot - what is in each file", strrep("=", 34), "",
           sprintf("%-18s %s", paste0(h[, 1], ":"), h[, 2]), "", "FILES", strrep("-", 5),
           sprintf("  %s  %s", format(basename(files)), desc[names(files)]), "")
  for (nm in names(tabs)) {
    cc <- .readme_columns(nm, tabs[[nm]], info)
    out <- c(out, sprintf("COLUMNS OF %s", basename(tab_files[[nm]])), strrep("-", 11 + nchar(basename(tab_files[[nm]]))),
             sprintf("  %-20s %s", cc[, 1], cc[, 2]), "")
  }
  out
}

.readme_sheet <- function(info, sheets, sheet_desc, files, desc) {
  h <- .readme_head(info)
  oth <- files[!names(files) %in% "xlsx"]
  rows <- rbind(h, c("", ""), c("SHEETS", ""), c("README", "this explanation"),
                cbind(names(sheets), unname(sheet_desc)))
  if (length(oth)) rows <- rbind(rows, c("", ""), c("OTHER FILES", ""), cbind(basename(oth), unname(desc[names(oth)])))
  for (nm in names(sheets)) {
    rows <- rbind(rows, c("", ""), c(sprintf("COLUMNS OF SHEET %s", nm), ""), .readme_columns(nm, sheets[[nm]], info))
  }
  data.frame(ITEM = rows[, 1], EXPLANATION = rows[, 2], stringsAsFactors = FALSE)
}

#' @export
print.quard_result <- function(x, ...) {
  cat("<quard_result>\n")
  print(x$region)
  cat("Lead SNP:            ", x$lead, "\n")
  cat("GWAS SNPs in region: ", nrow(x$gwas), "\n")
  if (!is.null(x$genes)) cat("Genes in region:     ", nrow(x$genes$genes), "\n")
  if (!is.null(x$ld)) cat("SNPs in LD heatmap:  ", nrow(x$ld), "\n")
  if (!is.null(x$ld_all)) cat("SNPs in LD table:    ", nrow(x$ld_all), "\n")
  if (length(x$files)) {
    cat("Files:\n")
    d <- x$file_descriptions %||% rep("", length(x$files))
    cat(paste0("  ", format(x$files), "  ", d, collapse = "\n"), "\n")
  }
  invisible(x)
}

#' @export
plot.quard_result <- function(x, ...) {
  if (is.null(x$plot)) {
    message("No figure was made (panels = \"none\").")
  } else {
    print(x$plot)
  }
  invisible(x)
}
