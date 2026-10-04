# Run quardplot from a terminal, for example:
#   Rscript quard_plot.R --species rice --gwas gwas.csv --gff genes.gff3 \
#           --vcf geno.vcf --target top --out results/locus.pdf
#   Rscript quard_plot.R --help
# Find this file with: system.file("scripts", "quard_plot.R", package = "quardplot")
quardplot::quard_cli()
