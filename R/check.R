# Input checks ---------------------------------------------------------------------

# Short list of chromosome names: "1, 2, 3, ..., 12".
.chr_list <- function(x) {
  u <- sort_chr(x)
  if (length(u) <= 8L) return(paste(u, collapse = ", "))
  paste(c(u[1:4], "...", utils::tail(u, 2L)), collapse = ", ")
}

# Normalised chromosome names in the first data lines of a VCF or HapMap file.
.first_chroms <- function(file, type, n = 5000L) {
  con <- file(file, "r")
  on.exit(close(con))
  x <- readLines(con, n = n + 2000L, warn = FALSE)
  if (type == "VCF") {
    x <- x[!startsWith(x, "#")]
    ch <- .field(x, 1L, "\t")
  } else {
    lay <- .hapmap_layout(file)
    x <- x[-1L]
    ch <- gsub("\"", "", .field(x, lay$r + 2L, lay$sep))
  }
  unique(normalize_chr(ch[nzchar(ch)]))
}

#' Check your input files
#'
#' Reads each file the way [quard_plot()] does and tells you, in plain
#' words, whether it has the expected structure and whether the files fit
#' together: chromosome names that can be matched, positions inside the
#' chromosomes of the species chosen with [set_species()] (a sign that all
#' files use the same reference genome version), and GWAS SNPs that are
#' present in the genotype files. Run it once for every new data set or
#' species, before making plots. The expected structure of every file is
#' described in the step-by-step guide (`quard_plot_examples.R`), and the
#' practice files written by [simulate_example_data()] show it as well.
#'
#' @param gwas GWAS results file, or a list with the elements `gwas`, `gff`,
#'   `vcf` and `hapmap` (for example `my` in the guide).
#' @param gff,vcf,hapmap Gene annotation (GFF3/GTF), VCF and HapMap files;
#'   each is optional.
#' @param window Half-width (bp) of the region around the most significant
#'   GWAS SNP that is read from the genotype files for the checks.
#' @param verbose Print the report.
#' @return Invisibly, a data frame with `file`, `status` (`"OK"`, `"NOTE"` or
#'   `"PROBLEM"`) and `message`.
#' @examples
#' \donttest{
#' ex <- simulate_example_data(tempfile("qp"))
#' set_species("rice")
#' check_inputs(list(gwas = ex$gwas, gff = ex$gff, vcf = ex$vcf, hapmap = ex$complete_hapmap))
#' }
#' @export
check_inputs <- function(gwas = NULL, gff = NULL, vcf = NULL, hapmap = NULL, window = 50000, verbose = TRUE) {
  if (is.list(gwas) && !is.data.frame(gwas)) {
    l <- gwas
    gwas <- l$gwas
    gff <- gff %||% l$gff
    vcf <- vcf %||% l$vcf
    hapmap <- hapmap %||% l$hapmap
  }
  rows <- list()
  add <- function(file, status, ...) rows[[length(rows) + 1L]] <<- c(file, status, paste0(...))
  run <- function(expr, file) {
    w <- character()
    val <- withCallingHandlers(
      tryCatch(expr, error = function(e) {
        add(file, "PROBLEM", conditionMessage(e))
        NULL
      }),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      },
      message = function(cnd) invokeRestart("muffleMessage"))
    for (x in unique(w)) add(file, "NOTE", x)
    val
  }
  sp <- get_species()
  sizes <- .species_sizes()
  genome <- if (!is.null(sp)) sp$genome else NULL
  if (is.null(sp)) {
    add("Species", "NOTE", "No species chosen: run set_species() first (e.g. set_species(\"rice\")) so that ",
        "positions can be checked against the chromosome lengths.")
  } else if (is.null(sizes)) {
    add("Species", "NOTE", sp$organism, ": no chromosome lengths given, so they are estimated from the data ",
        "and positions cannot be checked against the reference genome.")
  } else {
    add("Species", "OK", sp$organism, ", reference genome ", genome, " (", nrow(sizes), " chromosomes: ",
        .chr_list(sizes$chr), ").")
  }
  beyond <- function(chr, pos, what) {
    if (is.null(sizes)) return(invisible(NULL))
    len <- sizes$length[match(chr, sizes$chr)]
    over <- !is.na(len) & pos > len
    if (any(over)) {
      add(what, "PROBLEM", .fmt(sum(over)), " position(s) lie beyond the end of chromosome ",
          paste(utils::head(unique(chr[over]), 5L), collapse = ", "), " of ", genome,
          ": this file probably uses another reference genome version.")
    }
    off <- setdiff(unique(chr), sizes$chr)
    real <- off[grepl("^[0-9]+$", off) & off != "0"]
    if (length(real)) {
      add(what, "PROBLEM", "Chromosome(s) ", paste(utils::head(sort_chr(real), 8L), collapse = ", "),
          if (length(real) > 8L) ", ..." else "", " do not exist in ", genome,
          ": is the right species chosen in set_species()?")
    }
    other <- setdiff(off, real)
    if (length(other)) {
      add(what, "NOTE", "Sequence name(s) that are not chromosomes of ", genome, ": ",
          paste(utils::head(other, 8L), collapse = ", "), if (length(other) > 8L) ", ..." else "",
          ". They are left out of the genome-wide plots (normal for unanchored contigs such as Un, Sy or 0).")
    }
    invisible(NULL)
  }
  chr_sets <- list()

  # GWAS results ----------------------------------------------------------------
  gw <- NULL
  if (is.null(gwas)) {
    add("GWAS", "PROBLEM", "No GWAS file given; it is required for quard_plot().")
  } else {
    gw <- run(read_gwas(gwas, verbose = FALSE), "GWAS")
    if (!is.null(gw)) {
      cl <- attr(gw, "columns")
      add("GWAS", "OK", .fmt(nrow(gw)), " SNPs on ", length(unique(gw$chr)), " chromosome(s) (",
          .chr_list(gw$chr), ")",
          if (!is.null(cl)) sprintf("; columns used: SNP = '%s', chromosome = '%s', position = '%s', P value = '%s'",
                                    cl[["snp"]], cl[["chr"]], cl[["pos"]], cl[["p"]]) else "", ".")
      k <- which.min(gw$p)
      add("GWAS", "OK", "Most significant SNP: ", gw$snp[k], " (chromosome ", gw$chr[k], ", ", .fmt(gw$pos[k]),
          " bp, P = ", format(gw$p[k], digits = 3), ").")
      if (stats::median(gw$p) < 0.01) {
        add("GWAS", "NOTE", "Most P values are very small (median ", format(stats::median(gw$p), digits = 2),
            "): check that the P value column holds P values, not -log10(P).")
      }
      beyond(gw$chr, gw$pos, "GWAS")
      chr_sets$GWAS <- unique(gw$chr)
    }
  }

  # gene annotation -----------------------------------------------------------------
  if (!is.null(gff)) {
    genes <- run(read_gff(gff, verbose = FALSE), "GFF")
    if (!is.null(genes)) {
      G <- genes$genes
      cls <- table(factor(genes$features$class, levels = c("exon", "CDS", "UTR5", "UTR3", "UTR")))
      add("GFF", "OK", .fmt(nrow(G)), " genes (", .fmt(nrow(genes$transcripts)), " transcripts) on ",
          length(unique(G$chr)), " chromosome(s) (", .chr_list(G$chr), "); ", toupper(genes$format), " format.")
      if (all(genes$transcripts$type == "gene_only")) {
        add("GFF", "NOTE", "Only gene lines (no mRNA / exon / CDS records): genes are drawn as single boxes, ",
            "and SNPs inside genes are called 'genic' instead of exon or intron.")
      } else {
        add("GFF", "OK", "Gene structure: ", .fmt(cls[["exon"]]), " exon, ", .fmt(cls[["CDS"]]), " CDS and ",
            .fmt(cls[["UTR5"]] + cls[["UTR3"]] + cls[["UTR"]]), " UTR records.")
        if (cls[["CDS"]] == 0) add("GFF", "NOTE", "No CDS records: exon SNPs cannot be split into CDS and UTR.")
      }
      if (all(is.na(G$note) | !nzchar(G$note))) {
        add("GFF", "NOTE", "No gene descriptions (Note= or description= attributes): the description ",
            "columns of the tables will be empty.")
      }
      beyond(G$chr, G$end, "GFF")
      chr_sets$GFF <- unique(G$chr)
    }
  }

  # genotype files -------------------------------------------------------------------
  geno_check <- function(file, label) {
    if (label == "VCF") {
      h <- run(.vcf_header(file), label)
      if (is.null(h)) return(invisible(NULL))
      if (!length(h) || !startsWith(h[length(h)], "#CHROM")) {
        add(label, "PROBLEM", "No #CHROM header line: this does not look like a VCF file.")
        return(invisible(NULL))
      }
      ns <- length(strsplit(h[length(h)], "\t", fixed = TRUE)[[1]]) - 9L
      add(label, if (ns > 0L) "OK" else "PROBLEM", .fmt(ns), " individuals (sample columns) in the header.")
    } else {
      lay <- run(.hapmap_layout(file), label)
      if (is.null(lay)) return(invisible(NULL))
      add(label, "OK", .fmt(length(lay$samples)), " individuals; ",
          c(`\t` = "tab", `,` = "comma", `;` = "semicolon")[[lay$sep]], "-separated HapMap.")
    }
    ch <- run(.first_chroms(file, label), label)
    if (!is.null(ch)) chr_sets[[label]] <<- ch
    if (is.null(gw)) return(invisible(NULL))
    top <- gw[which.min(gw$p), ]
    g <- run(read_genotypes(file, chr = top$chr, start = max(1, top$pos - window), end = top$pos + window,
                            verbose = FALSE), label)
    if (is.null(g)) return(invisible(NULL))
    nv <- nrow(g$info)
    if (!nv) {
      add(label, "PROBLEM", "No variants within ", .fmt(window), " bp of the most significant GWAS SNP ",
          "(chromosome ", top$chr, ", ", .fmt(top$pos), " bp). Check that this file has the same chromosome ",
          "names and genome version as the GWAS file.")
      return(invisible(NULL))
    }
    gpos <- gw$pos[gw$chr == top$chr & abs(gw$pos - top$pos) <= window]
    hit <- mean(gpos %in% g$info$pos)
    add(label, if (hit >= 0.5) "OK" else if (hit > 0) "NOTE" else "PROBLEM", .fmt(nv), " variants within ",
        .fmt(window), " bp of the most significant GWAS SNP; ", round(100 * hit), "% of the GWAS SNPs there ",
        "have a genotype in this file",
        if (hit == 0) ": the positions do not match (another genome version, or other SNPs than in the GWAS?)." else ".")
    miss <- mean(is.na(g$dosage))
    add(label, if (miss > 0.5) "NOTE" else "OK", format(round(100 * miss, 1), nsmall = 1),
        "% missing genotype calls in that region", if (miss > 0.5) " (very high)." else ".")
    if (!(top$pos %in% g$info$pos)) {
      add(label, "NOTE", "The most significant SNP itself has no genotype in this file: LD will be shown ",
          "relative to the closest genotyped SNP.")
    }
    invisible(NULL)
  }
  if (!is.null(vcf)) geno_check(vcf, "VCF")
  if (!is.null(hapmap)) geno_check(hapmap, "HapMap")

  # chromosome names across files ----------------------------------------------------
  if (!is.null(chr_sets$GWAS)) {
    for (nm in setdiff(names(chr_sets), "GWAS")) {
      if (!length(intersect(chr_sets[[nm]], chr_sets$GWAS))) {
        add(nm, "PROBLEM", "Its chromosome names (", .chr_list(chr_sets[[nm]]), ") do not match those of the ",
            "GWAS file (", .chr_list(chr_sets$GWAS), "). Rename them, or translate them with chr_map ",
            "(see normalize_chr()).")
      }
    }
  }

  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("file", "status", "message")
  if (verbose) {
    cat("Checking your input files\n")
    for (i in seq_len(nrow(out))) {
      txt <- strwrap(out$message[i], width = 100L)
      cat(sprintf("%-8s %-8s %s\n", out$file[i], out$status[i], txt[1]))
      if (length(txt) > 1L) cat(paste0(strrep(" ", 18L), txt[-1L], "\n"), sep = "")
    }
    np <- sum(out$status == "PROBLEM")
    cat(if (np) sprintf("\nResult: %d PROBLEM(S) - please fix these first (see the messages above).\n", np) else
      "\nResult: no problems found - the files have the expected structure.\n")
  }
  invisible(out)
}
