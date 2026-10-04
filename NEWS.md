# quardplot 2.0.0

First release as an R package. Version 1.0.0 was the single script
`quard_plot.r` evaluated in the first submission of the paper; the main
function keeps its name, `quard_plot()`, and still accepts the argument names
of that script.

## Compared with the 1.0.0 script

* **Installable R package** with a help page for every function, a vignette,
  a step-by-step guide (`quard_guide()`), simulated example data
  (`simulate_example_data()`), automated tests and a command-line interface
  (`quard_cli()`).
* **LD is calculated within R as in PLINK 1.9**: genotype-based r2 as
  `--r2`, and haplotype-based r2 and D' (EM algorithm) as `--r2 dprime`.
  Missing calls stay missing, and variants are filtered by minor allele
  frequency (0.05) and missing rate (0.2) before LD is calculated.
* **No PLINK block file is needed**: LD blocks are found from the genotypes
  with the method of Gabriel et al. (as PLINK 1.9 `--blocks`); a PLINK
  `.blocks.det` file can still be given, as can a window or exact
  coordinates.
* **Any species without editing the code**: `set_species()` gives the
  chromosome lengths of rice, Arabidopsis, tomato, human and mouse genomes,
  or they are read from a FASTA index, the genome FASTA, GFF3 or VCF headers.
  Chromosome names such as `Chr1`, `chr01` and `1` are matched automatically.
* **Choose panels and tables**: `panels` selects the figure panels, `tables`
  the gene list, the SNP table (P value, r2 with the lead SNP, nearest gene
  and location in exon, intron, promoter, ...) and the LD table, written as
  CSV or one Excel file (`table_format = "xlsx"`).
* **Journal-ready figure**: 180 mm wide with 8-point text, vector PDF kept
  small by thinning overlapping Manhattan points, colour-blind-safe colours
  that keep their order in grey, and SNPs without LD data drawn as crosses.
* **Significance line**: Bonferroni by default, or any P or -log10 P value,
  plus an optional suggestive line.
* **Input files**: GWAS results of GAPIT, TASSEL, PLINK, GEMMA and rMVP are
  read without changes; GFF3 or GTF; VCF (plain, bgzip-compressed, optionally
  tabix-indexed) or HapMap. `check_inputs()` checks the files before
  plotting, and error messages say what to fix (including Git LFS
  placeholder files).
* **Extra tools**: gene lists for any regions (`get_genes()`), single-gene
  plots with SNPs (`plot_gene_structure()`), HapMap to VCF conversion
  (`hapmap_to_vcf()`), chromosome maps of SNP and indel density
  (`plot_karyotype()`) and trees of the samples with their names
  (`hapmap_tree()`).
* **Reproducibility**: the benchmark and the PLINK comparison of the paper
  can be repeated with the scripts in `inst/benchmark`, and every run of
  `quard_plot()` writes a README with the package and R versions and the
  settings used.
