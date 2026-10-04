#' quardplot: GWAS locus figures and tables
#'
#' `quardplot` turns GWAS results, a gene annotation (GFF3/GTF) and genotypes
#' (VCF or HapMap) into publication-ready locus figures and tables.
#'
#' Start here:
#' * [quard_guide()] copies the step-by-step guide (with plain-language
#'   explanations of every function) into your folder;
#' * [set_species()] chooses the species (do this first);
#' * [check_inputs()] checks that your files have the expected structure;
#' * [quard_plot()] makes the figure (PDF) and the tables (CSV or Excel).
#'
#' The figure has up to four panels on one aligned genomic axis: a
#' genome-wide Manhattan plot ([plot_manhattan()]), a regional association
#' plot coloured by linkage disequilibrium with the lead SNP
#' ([plot_regional()]), the gene models ([plot_genes()]) and an LD heatmap
#' ([plot_ld()]). LD (\eqn{r^2}, D') and LD blocks are calculated internally
#' ([calc_ld()], [find_ld_blocks()]) with the same definitions as PLINK 1.9.
#' Each SNP is placed relative to the genes: exon (CDS / UTR), intron,
#' promoter (upstream), downstream or intergenic ([annotate_variants()]).
#'
#' Other tools: gene lists for any region ([get_genes()]), single-gene
#' structure plots ([plot_gene_structure()]), HapMap to VCF conversion
#' ([hapmap_to_vcf()]), chromosome maps of SNP and indel density
#' ([plot_karyotype()]) and phylogenetic trees of the individuals
#' ([hapmap_tree()]). Tables are saved with [write_table()] as CSV, text or
#' Excel; figures with [save_pdf()].
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom ggplot2 .data
## usethis namespace: end
NULL
