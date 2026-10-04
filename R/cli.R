# Command-line interface ---------------------------------------------------------

#' Command-line interface
#'
#' Runs [quard_plot()] from a terminal. Every argument of [quard_plot()] can
#' be given as `--name value` (`-` or `_` both work), plus `--species`
#' (see [set_species()]). Example:
#'
#' ```
#' Rscript quard_plot_v3.R --species rice --gwas gwas.csv --gff genes.gff3 \
#'   --vcf geno.vcf.gz --target top --out results/locus.pdf \
#'   --panels manhattan,regional,genes --tables genes,snps,ld --table-format xlsx
#' ```
#'
#' `--out` is short for `--output`; `--panels` and `--tables` take
#' comma-separated lists (or `none`); `--quiet` turns messages off;
#' `--help` prints the options.
#'
#' @param args Character vector of arguments (default: the command line).
#' @return Invisibly, the `quard_result` from [quard_plot()].
#' @export
quard_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste(
    "Usage: Rscript -e \"quardplot::quard_cli()\" --gwas FILE [options]",
    "   or, with the single-file script: Rscript quard_plot_v3.R --gwas FILE [options]",
    "",
    "Input files",
    "  --gwas FILE            GWAS results (CSV/TSV; columns are recognised automatically,",
    "                         or name them with --snp-col, --chr-col, --pos-col, --p-col)",
    "  --gff FILE             gene annotation (GFF3/GTF)",
    "  --vcf FILE             genotypes (VCF or HapMap), needed for LD",
    "  --blocks FILE          optional LD-block file (PLINK .blocks.det or CHR START END)",
    "Species (chromosome lengths for the Manhattan plot)",
    "  --species NAME         rice, arabidopsis, tomato, human, mouse (or --chrom-sizes FILE)",
    "Which SNP",
    "  --chr CHR --pos BP     chromosome and position, or",
    "  --target ID            a SNP name, chr:pos, or top (= most significant SNP)",
    "Output",
    "  --out FILE.pdf         figure name; tables are named after it",
    "  --panels LIST          manhattan,regional,genes,ld  or  none (tables only)",
    "  --tables LIST          genes,snps,ld,ld_pairs  or  none   (default genes,snps,ld)",
    "  --table-format FMT     csv (default) or xlsx (one Excel file; needs the writexl package)",
    "Options",
    "  --threshold VALUE      bonferroni (default), a P value (5e-8) or -log10 value (7)",
    "  --region MODE          auto, blocks, ld_blocks or window; --flank BP; --start BP --end BP",
    "  --min-maf 0.05 --max-missing 0.2 --ld-r2 0.6 --max-ld-snps 300",
    "  --upstream-bp 3000 --downstream-bp 1000   promoter / downstream windows",
    "  --gene-features gene|all --transcripts canonical|all --quiet",
    "Any other argument of quard_plot() works the same way, e.g. --ld-method haplotype",
    "or --min-width 50000 (see ?quard_plot).",
    sep = "\n")
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  opts <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) .stop("Unexpected argument '", a, "'. Options start with --.\n", usage)
    a <- sub("^--", "", a)
    if (grepl("=", a, fixed = TRUE)) {
      key <- sub("=.*$", "", a)
      val <- sub("^[^=]*=", "", a)
      i <- i + 1L
    } else if (i == length(args) || startsWith(args[i + 1L], "--")) {
      key <- a
      val <- "TRUE"
      i <- i + 1L
    } else {
      key <- a
      val <- args[i + 1L]
      i <- i + 2L
    }
    opts[[gsub("-", "_", key)]] <- val
  }
  if (!is.null(opts$out)) {
    opts$output <- opts$out
    opts$out <- NULL
  }
  if (!is.null(opts$quiet)) {
    opts$verbose <- !as.logical(opts$quiet)
    opts$quiet <- NULL
  }
  if (!is.null(opts$species)) {
    set_species(opts$species, verbose = !identical(opts$verbose, FALSE))
    opts$species <- NULL
  }
  num <- c("pos", "start", "end", "flank", "min_width", "min_maf", "max_missing", "max_ld_snps",
           "max_table_snps", "block_search", "max_block_kb", "max_block_snps", "width", "height",
           "base_size", "ld_r2", "margin", "upstream_bp", "downstream_bp")
  for (k in intersect(names(opts), num)) opts[[k]] <- as.numeric(opts[[k]])
  lgl <- c("gene_file", "extend_to_genes", "verbose")
  for (k in intersect(names(opts), lgl)) opts[[k]] <- as.logical(opts[[k]])
  for (k in intersect(names(opts), c("panels", "tables"))) {
    opts[[k]] <- trimws(strsplit(opts[[k]], ",", fixed = TRUE)[[1]])
  }
  for (k in intersect(names(opts), c("threshold", "suggestive"))) {
    v <- suppressWarnings(as.numeric(opts[[k]]))
    if (!is.na(v)) opts[[k]] <- v
  }
  invisible(do.call(quard_plot, opts))
}
