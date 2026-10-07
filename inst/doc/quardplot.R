## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 8, dev = "png",
                      dpi = 96, fig.retina = 1)
if (Sys.info()[["sysname"]] == "Linux" && isTRUE(capabilities("cairo"))) options(bitmapType = "cairo")

## -----------------------------------------------------------------------------
library(quardplot)
list_species()[, c("species", "reference_genome", "chromosomes")]
set_species("rice")

## -----------------------------------------------------------------------------
ex <- simulate_example_data(file.path(tempdir(), "qp_vignette"))
my <- list(gwas = ex$gwas, gff = ex$gff, vcf = ex$vcf, hapmap = ex$complete_hapmap)
check_inputs(my)

## ----fig.height = 9-----------------------------------------------------------
out <- file.path(tempdir(), "qp_vignette_out")
res <- quard_plot(my$gwas, my$gff, my$vcf, chr = ex$chr, pos = ex$pos,
                  output = file.path(out, "locus.pdf"), verbose = FALSE)
res
res$plot

## -----------------------------------------------------------------------------
head(res$snp_table[, c("SNP", "POS", "P", "R2_WITH_LEAD", "GENE", "LOCATION", "FEATURE")])

## -----------------------------------------------------------------------------
r2 <- quard_plot(my$gwas, my$gff, my$vcf, chr = ex$chr, pos = ex$pos,
                 panels = c("regional", "genes"), tables = c("snps", "ld"),
                 output = file.path(out, "locus_no_heatmap.pdf"), verbose = FALSE)
names(r2$files)

## -----------------------------------------------------------------------------
g <- genotype_qc(read_vcf(my$vcf, chr = 9, verbose = FALSE), verbose = FALSE)
ld <- calc_ld(g)                       # r2 as in PLINK 1.9 --r2
blocks <- find_ld_blocks(g, verbose = FALSE)   # Gabriel method, as in PLINK --blocks
head(blocks[, c("chr", "start", "end", "kb", "nsnps")])

## -----------------------------------------------------------------------------
a <- annotate_variants(g, my$gff)
table(a$location)

## ----fig.height = 4-----------------------------------------------------------
genes <- get_genes(my$gff, chr = 9, pos = ex$pos, flank = 20000, verbose = FALSE)
genes[, c("gene_id", "start", "end", "strand", "distance_bp")]
plot_gene_structure(my$gff, genes$gene_id[1], flank = 2000, variants = my$vcf)

## ----fig.height = 6-----------------------------------------------------------
plot_karyotype(my$hapmap, bin_size = 2e6, verbose = FALSE)

## ----fig.height = 7-----------------------------------------------------------
tr <- hapmap_tree(my$hapmap, file.path(out, "tree.pdf"), k = 4, verbose = FALSE)
tr$plot

