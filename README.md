# quard_plot

A biologist-friendly GWAS-LD-Gene visualization tool.

`quard_plot` is an R script designed to simplify genomic data exploration by fusing genome-wide association study (GWAS) results, GFF3 gene annotations, and Variant Call Format (VCF) data into a highly interpretable, publication-quality 4-panel vector graphic. It allows biologists and agricultural researchers to rapidly pinpoint functional candidate genes within specific linkage disequilibrium (LD) blocks.

---

## Key Visual Panels

The tool automates base R graphic placement to stack four critical layers into a single `.pdf` output:
1. **Global Manhattan Plot:** Provides full genome-wide context while highlighting your target locus region using light-grey alignment paths.
2. **Regional Zoom-In Plot:** Focuses tightly on the target chromosome window to expose local association signals (`-log10(p)`).
3. **Gene Track:** Parses regional GFF3 files to draw clear gene structures, exon blocks, and coding directions ($+$ vs $-$ strands) for intuitive candidate screening.
4. **LD Heatmap Block:** Extracts genotypic variants from a VCF matrix, calculates pairwise correlation coefficients ($r^2$), and aligns the resulting diamond-mesh matrix precisely beneath the local signals.

---

## Required Setup

This script uses optimized core data tables. Before running the function, ensure you have installed the required dependencies in your R environment:

```R
install.packages(c("data.table", "utils", "stats", "grDevices", "graphics"))


Input File Formats
To ensure flawless parsing, prepare your data matrices to match the following structural formats:

GWAS File (.csv / Comma-separated): Contains structural mapping coordinates. The tool automatically detects single-string arrays (e.g., Chr1,1250432,1e-5) as well as clear multi-column headers.

GFF3 File (.gff / .gff3): Standard tab-delimited annotation file. Ensure the sequence IDs in column 1 use the standard naming scheme matching your target call (e.g., Chr1, Chr2).

VCF File (.vcf): Standard variant call matrix containing sample genotypes (e.g., 0/0, 0|1, 1|1). The script converts these positions into numerical alleles automatically.

Block File (.txt / .csv): A custom structured matrix identifying localized boundary ranges. It must include three specific header columns: CHR, START, and END.

Quick Start Example
Since quard_plot is contained within a single .r file, you don't need to configure a complex build system. Simply download or clone quard_plot.r, move it to your working directory, and run:


# 1. Source the visualization tool
source("quard_plot.r")

# 2. Configure paths to your genomic dataset files
gwas_input  <- "data/rice_gwas_results.csv"
gff_input   <- "data/Oryza_sativa.IRGSP-1.0.gff3"
vcf_input   <- "data/target_region.vcf"
block_input <- "data/ld_blocks.txt"

# 3. Generate your 4-panel locus figure (e.g., targeting a locus on Chromosome 1)
quard_plot(
  gwas_file        = gwas_input,
  gff_file         = gff_input,
  vcf_file         = vcf_input,
  block_file_path  = block_input,
  target_snp_chr   = 1,
  target_snp_pos   = 23405000,
  output_name      = "Chr1_Candidate_Locus_Plot.pdf"
)


+-------------------+-------------+---------------------------------------------------------------------------------+
| Variable Name     | Class       | Description                                                                     |
+-------------------+-------------+---------------------------------------------------------------------------------+
| gwas_file         | character   | File path pointing to your background GWAS data.                                |
| gff_file         | character   | File path pointing to your reference structural annotations.                   |
| vcf_file         | character   | File path pointing to your genotypic matrix data.                               |
| block_file_path   | character   | File path pointing to estimated physical linkage boundaries.                    |
| target_snp_chr    | numeric     | Target chromosome index (e.g., 1 to 12 for rice).                               |
| target_snp_pos    | numeric     | Precise base pair coordinate of your target variant.                            |
| output_name       | character   | The destination filename for your rendered high-resolution PDF file.           |
|                   |             | Default is "GWAS_Locus_LD_Plot.pdf".                                            |
+-------------------+-------------+---------------------------------------------------------------------------------+
