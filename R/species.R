# Species ------------------------------------------------------------------------

# Built-in species: organism, reference genome and the chromosome-length preset.
.qp_species <- list(
  rice = list(organism = "Rice (Oryza sativa, Nipponbare)", genome = "IRGSP-1.0 / MSU7",
              preset = "rice", aliases = c("oryza_sativa", "oryza sativa", "irgsp", "irgsp1", "irgsp-1.0",
                                           "msu7", "msu", "japonica")),
  arabidopsis = list(organism = "Arabidopsis thaliana (Col-0)", genome = "TAIR10",
                     preset = "arabidopsis", aliases = c("arabidopsis_thaliana", "arabidopsis thaliana",
                                                         "athaliana", "a. thaliana", "tair10", "tair")),
  tomato = list(organism = "Tomato (Solanum lycopersicum, Heinz 1706)", genome = "SL4.0 (ITAG4.0 / ITAG4.1)",
                preset = "tomato", aliases = c("solanum_lycopersicum", "solanum lycopersicum", "sl4.0", "sl4",
                                               "itag4", "itag4.0", "itag4.1")),
  human = list(organism = "Human (Homo sapiens)", genome = "GRCh38 / hg38",
               preset = "human_grch38", aliases = c("homo_sapiens", "homo sapiens", "grch38", "hg38",
                                                    "human_grch38")),
  human_grch37 = list(organism = "Human (Homo sapiens)", genome = "GRCh37 / hg19",
                      preset = "human_grch37", aliases = c("grch37", "hg19")),
  mouse = list(organism = "Mouse (Mus musculus)", genome = "GRCm39 / mm39",
               preset = "mouse_grcm39", aliases = c("mus_musculus", "mus musculus", "grcm39", "mm39",
                                                    "mouse_grcm39")),
  mouse_grcm38 = list(organism = "Mouse (Mus musculus)", genome = "GRCm38 / mm10",
                      preset = "mouse_grcm38", aliases = c("grcm38", "mm10"))
)

# Name of a built-in species for any accepted spelling, or NA.
.species_key <- function(x) {
  k <- tolower(trimws(as.character(x)))
  if (k %in% names(.qp_species)) return(k)
  for (nm in names(.qp_species)) if (k %in% .qp_species[[nm]]$aliases) return(nm)
  NA_character_
}

#' Built-in species
#'
#' Lists the species whose chromosome lengths are built in. For any other
#' species, use [set_species()] with your own chromosome lengths (for
#' example the `.fai` index of the reference genome), or let the lengths be
#' estimated from your data.
#'
#' @return A data frame: `species` (the name to use in [set_species()]),
#'   `organism`, `reference_genome`, `chromosomes`, `size_Mb`.
#' @examples
#' list_species()
#' @export
list_species <- function() {
  rows <- lapply(names(.qp_species), function(nm) {
    s <- .qp_species[[nm]]
    v <- .qp_presets[[s$preset]]
    data.frame(species = nm, organism = s$organism, reference_genome = s$genome,
               chromosomes = length(v), size_Mb = round(sum(v) / 1e6, 1), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Choose the species (do this first)
#'
#' Tells all functions which species, and so which chromosome lengths, you
#' are working with. After `set_species("rice")` the Manhattan plot,
#' karyotype plot and HapMap-to-VCF conversion use the rice chromosome
#' lengths without further arguments. A `chrom_sizes` argument given to a
#' function still wins over the species setting.
#'
#' Built-in species: see [list_species()] (rice, arabidopsis, tomato,
#' human, mouse). For any other species give the chromosome lengths with
#' `chrom_sizes`: the reference genome FASTA or its `.fai` index (from
#' `samtools faidx`), a GFF3 file with `##sequence-region` lines, a
#' two-column text file (chromosome, length) or a named vector. Without
#' lengths, they are estimated from your data (the largest SNP position on
#' each chromosome), which is fine for most plots.
#'
#' Your GWAS results, genotype files (VCF/HapMap) and GFF must all use the
#' same reference genome version as the species setting.
#'
#' @param species A built-in species (`"rice"`, `"arabidopsis"`, `"tomato"`,
#'   `"human"`, `"mouse"`, or a genome version such as `"hg19"`, `"mm10"`),
#'   any other name together with `chrom_sizes`, or `NULL` to clear the
#'   setting.
#' @param chrom_sizes Chromosome lengths for a species that is not built in
#'   (anything accepted by [chrom_sizes()]), or to replace the built-in ones.
#' @param verbose Print what was set.
#' @return Invisibly, the chromosome lengths (a data frame), or `NULL`.
#' @examples
#' set_species("rice")
#' get_species()
#' set_species("my_plant", chrom_sizes = c(Chr1 = 52e6, Chr2 = 47e6, Chr3 = 39e6))
#' set_species(NULL)
#' @export
set_species <- function(species, chrom_sizes = NULL, verbose = TRUE) {
  if (missing(species) || is.null(species) || tolower(species) %in% c("none", "")) {
    options(quardplot.species = NULL)
    .msg(verbose, "Species setting cleared: chromosome lengths now come from each function's ",
         "`chrom_sizes` argument or from your data.")
    return(invisible(NULL))
  }
  if (!is.character(species) || length(species) != 1L) .stop("`species` must be one name, e.g. \"rice\".")
  key <- .species_key(species)
  if (!is.na(key)) {
    info <- .qp_species[[key]]
    sizes <- if (is.null(chrom_sizes)) chrom_sizes(info$preset) else chrom_sizes(chrom_sizes)
    src_len <- if (is.null(chrom_sizes)) "built-in" else attr(sizes, "source")
    name <- key
  } else {
    info <- list(organism = species, genome = "your reference genome")
    sizes <- if (is.null(chrom_sizes)) NULL else chrom_sizes(chrom_sizes)
    src_len <- if (is.null(sizes)) NULL else attr(sizes, "source")
    name <- species
  }
  if (!is.null(sizes)) {
    attr(sizes, "source") <- sprintf("species '%s' (%s)", name, info$genome)
  }
  options(quardplot.species = list(name = name, organism = info$organism, genome = info$genome,
                                   sizes = sizes))
  if (isTRUE(verbose)) {
    if (!is.null(sizes)) {
      message(sprintf("Species: %s - reference genome %s (%d chromosomes, %s Mb; lengths: %s).",
                      info$organism, info$genome, nrow(sizes), format(round(sum(sizes$length) / 1e6, 1), nsmall = 1),
                      src_len))
      message("All plots now use these chromosome lengths. Your GWAS, genotype (VCF/HapMap) and GFF ",
              "files must use the same reference genome.")
    } else {
      message("Species: ", species, ". It is not built in (see list_species()), so chromosome lengths ",
              "will be estimated from your data (largest SNP position on each chromosome).")
      message("For exact lengths: set_species(\"", species, "\", chrom_sizes = \"genome.fa.fai\") ",
              "(or the genome FASTA, a GFF3 with ##sequence-region lines, or a chromosome/length table).")
    }
  }
  invisible(sizes)
}

#' The species that is set
#'
#' @return A list (`name`, `organism`, `genome`, `sizes`) or `NULL` if no
#'   species was set with [set_species()].
#' @examples
#' get_species()
#' @export
get_species <- function() {
  getOption("quardplot.species")
}

# Chromosome lengths of the species set with set_species(), or NULL.
.species_sizes <- function() {
  s <- getOption("quardplot.species")
  if (is.null(s)) NULL else s$sizes
}
