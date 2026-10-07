# =============================================================================
#  quard_plot 2.0.0 - STEP-BY-STEP GUIDE
#
#  This file shows every function of quard_plot, using YOUR data files.
#  You do not need to know R well:
#    * Lines that start with # are explanations; R ignores them.
#    * Run the file from top to bottom (RStudio: the "Source" button), or put
#      the cursor on a line and press Ctrl+Enter (Mac: Cmd+Enter) to run the
#      lines one at a time and see what each one does.
#    * Everything is saved in the folder "my_results".
#
#  WORDS USED IN THIS GUIDE
#    SNP          a single-base DNA difference between individuals (a marker).
#                 An indel is a short insertion or deletion.
#    GWAS         genome-wide association study: every SNP gets a P value.
#                 The smaller P, the stronger the link with the trait.
#    -log10(P)    the height of a SNP in the plots: P = 0.001 -> 3, P = 1e-8 -> 8.
#    Lead SNP     the SNP you are looking at, usually the most significant one.
#    LD           linkage disequilibrium: how often two SNPs are inherited
#                 together. r2 = 1: always together; r2 = 0: independent.
#    LD block     a stretch of chromosome that is inherited as one piece
#                 (all SNPs in strong LD). The figure shows the block of your SNP.
#    GFF file     the gene annotation: where every gene, exon, intron and UTR is.
#    VCF, HapMap  genotype files: which variant every individual carries at
#                 every SNP. They are needed to calculate LD.
#    Exon/intron  parts of a gene that are kept in / removed from the mRNA.
#    CDS, UTR     CDS = protein-coding part of the exons; UTR = untranslated
#                 part at the start (5'UTR) or end (3'UTR) of the mRNA.
#    Promoter     the region just before the start of a gene (here: 3,000 bp);
#                 SNPs there are called "upstream".
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 0. Load quardplot
# -----------------------------------------------------------------------------
# Needed once (see README): install.packages(c("data.table", "ggplot2", "patchwork"))
#                           install.packages("quardplot_2.0.0.tar.gz", repos = NULL, type = "source")
# Optional:  install.packages("writexl")  # Excel (.xlsx) files; without it
#                                         # you get CSV files (they also open in Excel)
#            install.packages("ape")      # faster trees for more than 400 samples
# Put this file and your data files in one folder and make it the working
# directory (RStudio: Session > Set Working Directory > To Source File Location).
library(quardplot)                        # loads all quardplot functions
getwd()                                   # the working directory (the folder R is looking in)


# -----------------------------------------------------------------------------
# STEP 1. Choose your species  (always do this first)
# -----------------------------------------------------------------------------
# The species tells the functions the chromosome lengths of your reference
# genome (for the genome-wide plots). Built-in species:
list_species()

set_species("rice")                       # rice, Nipponbare IRGSP-1.0 / MSU7

# Other built-in choices:
#   set_species("arabidopsis")            # TAIR10
#   set_species("tomato")                 # SL4.0 (ITAG4.0 / ITAG4.1)
#   set_species("human")                  # GRCh38;  set_species("hg19") for GRCh37
#   set_species("mouse")                  # GRCm39;  set_species("mm10") for GRCm38
# Any other species: give the chromosome lengths of its reference genome,
#   set_species("pepper", chrom_sizes = "pepper_genome.fa.fai")  # FASTA index (samtools faidx)
#   set_species("pepper", chrom_sizes = "pepper_genome.fa")      # or the genome FASTA itself
#   set_species("pepper")       # or no lengths: they are estimated from your data
# IMPORTANT: your GWAS, VCF/HapMap and GFF files must all use the same
# reference genome version.


# -----------------------------------------------------------------------------
# STEP 2. Your files: what they must look like, and where they are
# -----------------------------------------------------------------------------
# Your files must have the SAME STRUCTURE as the rice files below. Files of
# any species that follow this structure work without any change, so prepare
# them like this first (the practice files of simulate_example_data() have
# exactly this structure and can be used as templates).
#
# 1. GWAS RESULTS (my$gwas): a table with one row per SNP (.csv, .txt or .gz;
#    comma-, tab- or space-separated). Four columns are needed; extra columns
#    are ignored. The column names are recognised automatically:
#      SNP name    SNP, rs, rsid, marker, variant_id, ID or Taxa
#      chromosome  Chr, CHR, chrom, chromosome or seqid
#      position    Pos, BP, position, ps or base_pair_location
#      P value     P.value, P, p_value, pvalue, pval, p_wald or p_lrt
#    Rice example (GAPIT output, Internal_CV.csv):
#      SNP,Chr,Pos,P.value,MAF,nobs,H&B.P.Value,Effect
#      S1_4734,1,4734,0.944,0.495,120,0.9968,0.036
#      S1_42827,1,42827,0.6736,0.229,120,0.9834,0.0203
#    * P values must be normal P values (between 0 and 1), not -log10(P).
#    * Other column names? Tell read_gwas() which ones to use, e.g.
#        gw <- read_gwas("x.csv", snp_col = "Marker", chr_col = "Chrom",
#                        pos_col = "Position", p_col = "Pvalue")
#      and give gw to quard_plot() instead of the file name.
#
# 2. GENE ANNOTATION (my$gff): a GFF3 file (GTF also works), 9 columns
#    separated by TABS: chromosome, source, type, start, end, score, strand,
#    phase, attributes. Every gene has a "gene" line with ID=...; its
#    transcripts ("mRNA") have Parent=<gene ID>; the "exon", "CDS",
#    "five_prime_UTR" and "three_prime_UTR" lines have Parent=<mRNA ID>.
#    Note=... or description=... gives the gene description.
#    Rice example (MSU7, msu_rice.gff3):
#      Chr9  MSU_osa1r7  gene  12403000  12406546  .  +  .  ID=LOC_Os09g20010;Name=LOC_Os09g20010;Note=MYB family transcription factor
#      Chr9  MSU_osa1r7  mRNA  12403000  12406546  .  +  .  ID=LOC_Os09g20010.1;Parent=LOC_Os09g20010
#      Chr9  MSU_osa1r7  exon  12403000  12403278  .  +  .  ID=LOC_Os09g20010.1:exon_1;Parent=LOC_Os09g20010.1
#      Chr9  MSU_osa1r7  CDS   12403150  12403278  .  +  0  Parent=LOC_Os09g20010.1
#    (A GFF with gene lines only also works: genes are then drawn as boxes.)
#
# 3. GENOTYPES, VCF (my$vcf): a standard VCF file (.vcf or .vcf.gz). After the
#    ## header lines comes the line "#CHROM POS ID REF ALT QUAL FILTER INFO
#    FORMAT" followed by one column per individual. Only the GT part of each
#    genotype is used: 0/0, 0/1, 1/1 (or 0|1), ./. = missing; AD, DP, ... are
#    ignored. Rice example (inputfinal.vcf):
#      #CHROM  POS   ID        REF  ALT  QUAL  FILTER  INFO  FORMAT         Acc001                 Acc002
#      Chr1    5416  S_1_5416  G    A    .     PASS    .     GT:AD:DP:GQ:PL 1/1:5,3:8:30:0,30,300  0/0:6,0:6:30:0,30,300
#    Large VCF? Compress it with bgzip and index it with tabix: much faster.
#
# 4. GENOTYPES, HAPMAP (my$hapmap): the TASSEL/GAPIT HapMap format, tab- or
#    comma-separated: 11 columns (rs#, alleles, chrom, pos, strand, assembly#,
#    center, protLSID, assayLSID, panelLSID, QCcode), then one column per
#    individual. Genotypes as one letter (A, C, G, T; R, Y, S, W, K, M =
#    heterozygous; N = missing) or two letters (AA, AG, NN); + and - for
#    insertions/deletions. Rice example (geno_complete.csv):
#      rs#,alleles,chrom,pos,strand,assembly#,center,protLSID,assayLSID,panelLSID,QCcode,Acc001,Acc002
#      S_1_5416,G/A,1,5416,+,NA,NA,NA,NA,NA,NA,A,R
#    You need a VCF OR a HapMap file for LD; giving both is not necessary.
#
# RULES FOR ALL FILES
#   * All files use the same reference genome version as the species of
#     step 1 (rice: IRGSP-1.0 = MSU7). Positions are in bp, counted from 1.
#   * Chromosome names may be written differently in each file: Chr1, chr01
#     and 1 are matched automatically. Accession-style names (e.g.
#     NC_029256.1) must be renamed, or translated with chr_map (see normalize_chr()).
#   * The genotypes must come from the same individuals (population) as the GWAS.
#
# OPTIONAL FILES (only for some functions)
#   * LD blocks from PLINK (quard_plot(blocks = ...)): the .blocks.det file of
#     plink --blocks, or a table with the columns CHR, START, END.
#   * Chromosome lengths of another species (set_species(..., chrom_sizes = ...)):
#     the .fai index of the genome (samtools faidx genome.fa), the genome FASTA,
#     or a table with two columns: chromosome, length.
#   * Groups of individuals for the tree (hapmap_tree(groups = ...)): two
#     columns, the individual name as in the VCF/HapMap and its group:
#        Taxa,Group
#        Acc001,indica
#        Acc002,japonica

my <- list(gwas   = "Internal_CV.csv",      # GWAS results
           gff    = "msu_rice.gff3",        # gene annotation (GFF3 or GTF)
           vcf    = "inputfinal.vcf",       # genotypes, VCF format (.vcf or .vcf.gz)
           hapmap = "geno_complete.csv")    # genotypes, HapMap format (tab- or comma-separated)
# Files in another folder: write the full path with / (not \), for example
#   gwas = "D:/GWAS/Internal_CV.csv"

# No data at hand? Change FALSE to TRUE to practise with a small simulated
# rice data set (written to the folder "practice_data").
USE_PRACTICE_DATA <- FALSE
if (USE_PRACTICE_DATA) {
  ex <- simulate_example_data("practice_data")
  my <- list(gwas = ex$gwas, gff = ex$gff, vcf = ex$vcf, hapmap = ex$complete_hapmap)
}

out <- "my_results"                        # all results are saved in this folder
dir.create(out, showWarnings = FALSE)

not_found <- unlist(my)[!file.exists(unlist(my))]
if (length(not_found)) {
  stop("Not found in the folder ", getwd(), ": ", paste(not_found, collapse = ", "))
}

# Check that your files have the expected structure and fit together
# (same chromosome names, same genome version, GWAS SNPs in the genotype
# files). Fix every line marked PROBLEM before you continue.
check_inputs(my)


# -----------------------------------------------------------------------------
# STEP 3. Look at the GWAS results and choose the SNP to study
# -----------------------------------------------------------------------------
gw <- read_gwas(my$gwas)        # SNP, chromosome, position and P columns are found automatically
head(gw)                        # the first rows: snp, chr, pos, p, logp (= -log10 P)
gwas_threshold(gw)              # default significance line: Bonferroni = 0.05 / number of SNPs

top <- gw[which.min(gw$p), ]    # the most significant SNP
top

CHR <- top$chr                  # chromosome of the SNP to study  (or type it, e.g. CHR <- 9)
POS <- top$pos                  # position of the SNP to study    (or type it, e.g. POS <- 12559349)
if (USE_PRACTICE_DATA) {        # practice data: use the simulated causal SNP
  CHR <- ex$chr
  POS <- ex$pos
}

genes <- read_gff(my$gff, chr = CHR)                              # the genes of that chromosome
near <- get_genes(genes, chr = CHR, pos = POS, flank = 100000)    # genes within 100 kb of the SNP
near                                                              # distance_bp = distance to the SNP
GENE <- if (nrow(near)) near$gene_id[which.min(near$distance_bp)] else NA   # the closest gene (step 10)
# To study another gene, type its ID as in the GFF file (see the gene_id column of near), e.g.
# GENE <- "LOC_Os09g20020"     # check: GENE %in% genes$genes$gene_id should give TRUE
GENE


# -----------------------------------------------------------------------------
# STEP 4. The main figure and tables: quard_plot()
# -----------------------------------------------------------------------------
# For the SNP you chose, quard_plot() makes:
#
#   A PDF FIGURE. Choose the panels with `panels =`:
#     "manhattan"  all SNPs of the genome (the studied region is marked)
#     "regional"   the SNPs of the region, coloured by LD (r2) with the lead SNP
#     "genes"      the genes of the region (boxes = exons, lines = introns)
#     "ld"         the LD heatmap: r2 between every pair of SNPs
#     "none"       no figure, only tables
#
#   TABLES, named after the PDF. Choose them with `tables =`:
#     "genes"      ..._genes.txt     the genes of the region, as lines of the GFF file
#     "snps"       ..._SNPs.csv      every GWAS SNP of the region: P value, r2 with the
#                                    lead SNP, closest gene, and whether the SNP lies in
#                                    an exon, intron, promoter (upstream), downstream
#                                    region or between genes (intergenic)
#     "ld"         ..._LD_r2.csv     r2 between every pair of SNPs (a matrix)
#     "ld_pairs"   ..._LD_r2_pairs.csv  the same r2 values, one row per pair of SNPs
#     "none"       no tables
#   plus ..._README.txt, which explains every column.
#   The default is tables = c("genes", "snps", "ld").
#
#   EXCEL: table_format = "xlsx" puts all tables in ONE Excel file
#   (..._tables.xlsx, one sheet per table + a README sheet).
#
# The tables do not depend on the panels: you can leave the LD heatmap out of
# the figure and still get the LD table (example 4b).
#
# Which region is shown? The LD block that contains your SNP (found from the
# genotypes). If the SNP is in no block: all SNPs with r2 >= 0.6 with it.
# Without a genotype file: 50 kb on each side. Or set start = and end =.

# 4a. Everything automatic: 4 panels + gene list + SNP table + LD table
res <- quard_plot(gwas = my$gwas, gff = my$gff, vcf = my$vcf, chr = CHR, pos = POS,
                  output = file.path(out, "4a_locus.pdf"))
res                                     # summary: region, lead SNP and the files written
head(res$snp_table, 20)                 # the SNP table (also saved as 4a_locus_SNPs.csv)
table(res$snp_table$LOCATION)           # number of SNPs in exons, introns, promoters, ...

# 4b. Figure WITHOUT the LD heatmap, but WITH the LD table
quard_plot(my$gwas, my$gff, my$vcf, chr = CHR, pos = POS,
           panels = c("manhattan", "regional", "genes"),     # no "ld" = no heatmap in the figure
           tables = c("genes", "snps", "ld"),                # the LD table is still written
           output = file.path(out, "4b_no_heatmap.pdf"))

# 4c. All tables in one Excel file
quard_plot(my$gwas, my$gff, my$vcf, chr = CHR, pos = POS,
           tables = c("genes", "snps", "ld", "ld_pairs"),
           table_format = "xlsx",                            # -> 4c_excel_tables.xlsx
           output = file.path(out, "4c_excel.pdf"))
# (Without the writexl package you get CSV files instead, and a message.)

# 4d. Only the tables, no figure
quard_plot(my$gwas, my$gff, my$vcf, chr = CHR, pos = POS,
           panels = "none", tables = c("snps", "ld"), table_format = "xlsx",
           output = file.path(out, "4d_tables_only"))

# 4e. Only some panels (here: no genome-wide Manhattan plot)
quard_plot(my$gwas, my$gff, my$vcf, chr = CHR, pos = POS,
           panels = c("regional", "genes", "ld"),
           output = file.path(out, "4e_three_panels.pdf"))

# 4f. Your own region (here 100 kb on each side) and your own significance line
quard_plot(my$gwas, my$gff, my$vcf, chr = CHR, pos = POS,
           start = POS - 100000, end = POS + 100000,
           threshold = 1e-7,          # "bonferroni" (default), a P value (1e-7) or -log10 P (7)
           suggestive = 1e-5,         # optional second, dotted line
           output = file.path(out, "4f_200kb.pdf"))

# 4g. Let quard_plot pick the most significant SNP, or name a SNP
quard_plot(my$gwas, my$gff, my$vcf, target = "top",
           output = file.path(out, "4g_top_snp.pdf"))
quard_plot(my$gwas, my$gff, my$vcf, target = top$snp, panels = c("regional", "genes"),
           output = file.path(out, "4g_by_name.pdf"))

# 4h. Genotypes from the HapMap file instead of the VCF
quard_plot(my$gwas, my$gff, my$hapmap, chr = CHR, pos = POS,
           output = file.path(out, "4h_hapmap.pdf"))

# 4i. Without a genotype file (no LD: Manhattan, regional and gene panels)
quard_plot(my$gwas, my$gff, chr = CHR, pos = POS,
           output = file.path(out, "4i_no_genotypes.pdf"))

# 4j. More settings
quard_plot(my$gwas, my$gff, my$vcf, chr = CHR, pos = POS,
           ld_stat = "dprime",          # D' instead of r2 in the heatmap and LD table
           transcripts = "all",         # draw all transcripts of each gene, not one per gene
           gene_features = "all",       # gene list with the gene, mRNA, exon, CDS and UTR lines
           extend_to_genes = TRUE,      # widen the region so genes at its borders are complete
           ld_r2 = 0.8,                 # stricter LD interval (used when the SNP is in no block)
           upstream_bp = 2000,          # promoter = 2,000 bp before a gene (default 3,000)
           downstream_bp = 500,         # downstream = 500 bp after a gene (default 1,000)
           min_maf = 0.05, max_missing = 0.2,   # genotype filters before LD (the defaults)
           output = file.path(out, "4j_more_settings.pdf"))

# 4k. Size of the figure (mm) and of the text (pt)
quard_plot(my$gwas, my$gff, my$vcf, chr = CHR, pos = POS,
           width = 170, height = 210, base_size = 9,
           output = file.path(out, "4k_size.pdf"))

# 4l. With an LD-block file made by PLINK (plink --blocks), or LD values from PLINK:
# quard_plot(my$gwas, my$gff, my$vcf, blocks = "plink.blocks.det", chr = CHR, pos = POS,
#            output = file.path(out, "4l_plink_blocks.pdf"))
# quard_plot(my$gwas, my$gff, chr = CHR, pos = POS, ld_file = "plink.ld", ld_bim = "plink.bim",
#            output = file.path(out, "4l_plink_ld.pdf"))


# -----------------------------------------------------------------------------
# STEP 5. Read and check your input files
# -----------------------------------------------------------------------------
chrom_sizes()                             # chromosome lengths of the species from step 1
try(chrom_sizes(my$vcf))                  # lengths from the VCF header; an error here only means
                                          # that the VCF header has no lengths (that is fine)
normalize_chr(c("Chr1", "chr01", "1", "SL4.0ch01"))   # chromosome names are matched automatically

genes                                     # the genes read in step 3
g_vcf <- read_vcf(my$vcf, chr = CHR, start = POS - 100000, end = POS + 100000)        # 200 kb of the VCF
g_hmp <- read_hapmap(my$hapmap, chr = CHR, start = POS - 100000, end = POS + 100000)  # same, from HapMap
g_vcf                                     # number of SNPs and individuals

# Genotype filters used before LD: two alleles, at most 20% missing data,
# minor allele frequency (MAF) at least 5%.
g <- genotype_qc(g_vcf, min_maf = 0.05, max_missing = 0.2)
g$qc                                      # how many SNPs each filter removed


# -----------------------------------------------------------------------------
# STEP 6. LD and LD blocks on their own
# -----------------------------------------------------------------------------
lead_id <- g$info$id[which.min(abs(g$info$pos - POS))]    # the genotyped SNP closest to POS
r2_lead <- ld_with_lead(g, lead_id)                       # r2 of every SNP with that SNP
head(sort(r2_lead, decreasing = TRUE), 10)                # the SNPs in strongest LD with it

ld <- calc_ld(g)                                          # r2 between all pairs (same as PLINK --r2)
ld_dp <- calc_ld(g, stat = "dprime")                      # D' (same as PLINK --r2 dprime)
write_ld(ld, file.path(out, "6_LD_r2_matrix.csv"))                        # the matrix as CSV
write_ld(ld, file.path(out, "6_LD_r2_with_genes.xlsx"), genes = genes)   # + gene location (Excel)
write_ld(ld, file.path(out, "6_LD_pairs_r2_above_0.8.csv"), format = "pairs", min_value = 0.8)

blocks <- find_ld_blocks(g)                               # LD blocks (same method as PLINK --blocks)
blocks[, c("chr", "start", "end", "kb", "nsnps")]
write_blocks(blocks, file.path(out, "6_blocks.blocks.det"))
region <- define_region(CHR, POS, blocks = blocks, geno = g)   # the block (region) of your SNP
region

# To compare with PLINK (for the paper), run PLINK 1.9 on the same region first:
#   plink --vcf inputfinal.vcf --chr 9 --from-bp START --to-bp END --maf 0.05 --geno 0.2 \
#         --make-bed --out locus --allow-extra-chr --double-id
#   plink --bfile locus --r2 square --out locus --allow-extra-chr
# compare_ld(calc_ld(g), read_plink_ld("locus.ld", bim = "locus.bim"))


# -----------------------------------------------------------------------------
# STEP 7. Where do SNPs lie relative to genes?
# -----------------------------------------------------------------------------
# For each SNP: the closest gene, and whether the SNP is in an exon (CDS,
# 5'UTR or 3'UTR), an intron, the promoter (upstream, 3 kb), the downstream
# region (1 kb) or between genes (intergenic), with the distance to the gene.
ann <- annotate_variants(g, genes)                        # all SNPs of the region
head(ann)
table(ann$location)
annotate_variants(c(POS, POS - 2000, POS + 2000), genes, chr = CHR)   # any positions

sig <- gw[gw$p <= gwas_threshold(gw)$p, ]                 # all significant SNPs, all chromosomes
if (nrow(sig)) {
  sig_ann <- annotate_variants(sig, my$gff)               # the P values are kept
  print(head(sig_ann))
  write_table(sig_ann, file.path(out, "7_significant_SNPs_genes.xlsx"))   # or ".csv"
}
best100 <- gw[order(gw$p)[1:min(100, nrow(gw))], ]       # or the 100 best SNPs,
annotate_variants(best100, my$gff, upstream_bp = 2000, downstream_bp = 500)[1:5, ]   # other windows


# -----------------------------------------------------------------------------
# STEP 8. Single panels and your own combinations
# -----------------------------------------------------------------------------
# Every panel is a separate function. Here they use the results of 4a (res).
p_man <- plot_manhattan(gw, highlight = res$region, target = res$lead)   # genome-wide
p_con <- plot_zoom_connector(p_man, res$region)                         # zoom lines
p_reg <- plot_regional(gw, res$region, lead = res$lead, r2 = res$r2_lead)   # the region
p_gen <- plot_genes(genes, res$region)                                  # the genes
p_ld  <- plot_ld(res$ld, res$region, lead = res$lead, blocks = res$blocks)   # LD heatmap

save_pdf(p_man, file.path(out, "8_manhattan.pdf"), width = 180, height = 70)
save_pdf(combine_panels(p_reg, p_gen, heights = c(2, 1)),               # any panels, any order
         file.path(out, "8_regional_and_genes.pdf"), width = 180, height = 110)

# Add your own track, e.g. QTL found before (start, end, label):
qtl <- data.frame(start = POS - 20000, end = POS + 20000, label = "my QTL")
p_qtl <- plot_track(qtl, res$region, title = "QTL")
save_pdf(combine_panels(p_man, p_con, p_reg, p_qtl, p_gen, p_ld, heights = c(40, 7, 45, 10, 25, 70)),
         file.path(out, "8_custom_figure.pdf"), width = 180, height = 210)

# Panels are ggplot2 objects, so titles and fonts can be changed:
p_title <- p_reg + ggplot2::ggtitle("My locus") + theme_quard(base_size = 10)
save_pdf(p_title, file.path(out, "8_regional_with_title.pdf"), width = 120, height = 70)


# -----------------------------------------------------------------------------
# STEP 9. Gene lists for any region
# -----------------------------------------------------------------------------
get_genes(my$gff, chr = CHR, start = POS - 200000, end = POS + 200000,
          output = file.path(out, "9_genes_400kb.xlsx"))                 # table (.xlsx, .csv or .txt)
get_genes(my$gff, chr = CHR, start = POS - 200000, end = POS + 200000,
          output = file.path(out, "9_genes_400kb.gff3"), format = "gff")  # original GFF lines
# Several regions at once (e.g. all your QTL):
my_regions <- data.frame(chr = CHR, start = c(POS - 300000, POS + 100000),
                         end = c(POS - 100000, POS + 300000), name = c("left", "right"))
get_genes(my$gff, regions = my_regions, output = file.path(out, "9_genes_two_regions.csv"))


# -----------------------------------------------------------------------------
# STEP 10. Structure of one gene, with its SNPs
# -----------------------------------------------------------------------------
# GENE was chosen in step 3 (the gene closest to your SNP). For another gene, type
# its ID here, e.g. GENE <- "LOC_Os09g20020"  (a gene ID, not a transcript ID such as ...20020.1)
if (!is.na(GENE)) {
  p1 <- plot_gene_structure(genes, GENE)                         # exons, introns, UTRs
  save_pdf(p1, file.path(out, "10_gene_structure.pdf"), width = 170, height = 70)

  p2 <- plot_gene_structure(genes, GENE, flank = 2000,           # + 2 kb on each side
                            variants = my$vcf, gwas = gw)        # + SNPs and their P values
  save_pdf(p2, file.path(out, "10_gene_with_snps.pdf"), width = 170, height = 100)
  print(attr(p2, "variants"))                                    # each SNP: promoter, CDS, intron, ...

  p3 <- plot_gene_structure(genes, GENE, flank = 2000, variants = my$hapmap,
                            coordinates = "relative")            # positions from the gene start
  save_pdf(p3, file.path(out, "10_gene_hapmap.pdf"), width = 170, height = 100)
}


# -----------------------------------------------------------------------------
# STEP 11. Convert HapMap to VCF
# -----------------------------------------------------------------------------
hapmap_to_vcf(my$hapmap, file.path(out, "11_region_from_hapmap.vcf"),
              chr = CHR, start = POS - 500000, end = POS + 500000, chr_prefix = "Chr")
# The whole file (can take a few minutes). With the reference genome FASTA,
# REF is the true reference base:
# hapmap_to_vcf(my$hapmap, file.path(out, "geno_complete.vcf.gz"), chr_prefix = "Chr",
#               ref_fasta = "IRGSP-1.0_genome.fasta")


# -----------------------------------------------------------------------------
# STEP 12. Chromosome map: SNP (red) and indel (yellow) density
# -----------------------------------------------------------------------------
k1 <- plot_karyotype(my$hapmap, bin_size = 1e6)                  # 1-Mb windows
save_pdf(k1, file.path(out, "12_chromosome_map.pdf"), width = 180, height = 130)
head(attr(k1, "density"))                                        # counts per window
k2 <- plot_karyotype(my$vcf, bin_size = 1e6, layout = "horizontal")   # from the VCF
save_pdf(k2, file.path(out, "12_chromosome_map_vcf.pdf"), width = 180, height = 150)


# -----------------------------------------------------------------------------
# STEP 13. Phylogenetic tree of the individuals (like GAPIT)
# -----------------------------------------------------------------------------
# Reads up to 10,000 evenly spaced SNPs, calculates genetic distances and
# draws a neighbour-joining tree with the name of every individual at its
# leaf. k = 4 colours 4 groups of related individuals.
# The PDF grows with the number of individuals so that every name is
# readable (6 pt); zoom in on the PDF to read them.
tr <- hapmap_tree(my$hapmap, file.path(out, "13_tree_circular.pdf"), k = 4)
hapmap_tree(my$hapmap, file.path(out, "13_tree_rectangular.pdf"), layout = "rectangular", k = 4)
hapmap_tree(my$vcf, file.path(out, "13_tree_from_vcf.pdf"), layout = "unrooted", k = 4)
# More settings for the names:
hapmap_tree(my$hapmap, file.path(out, "13_tree_big_names.pdf"), k = 4,
            label_size = 8,               # font size of the names (pt)
            label_color = "group")        # names in the colour of their group
# hapmap_tree(..., width = 180, height = 200)   # fixed page size: names shrink to fit
# hapmap_tree(..., labels = FALSE)              # no names
# Colour by your own groups (a file with two columns: individual, group):
# hapmap_tree(my$hapmap, file.path(out, "13_tree_groups.pdf"), groups = "subpopulations.csv")
write_newick(tr$tree, file.path(out, "13_tree.nwk"))   # open in FigTree, iTOL or MEGA


# -----------------------------------------------------------------------------
# STEP 14. From a terminal (without opening R)
# -----------------------------------------------------------------------------
#   Rscript -e "quardplot::quard_cli()" --help
#   Rscript -e "quardplot::quard_cli()" --species rice --gwas Internal_CV.csv --gff msu_rice.gff3 \
#           --vcf inputfinal.vcf --target top --out my_results/14_terminal.pdf \
#           --panels manhattan,regional,genes --tables genes,snps,ld --table-format xlsx
# The same from R:
quard_cli(c("--species", "rice", "--gwas", my$gwas, "--gff", my$gff, "--vcf", my$vcf,
            "--target", "top", "--out", file.path(out, "14_terminal.pdf"),
            "--panels", "manhattan,regional,genes", "--tables", "genes,snps,ld"))

cat("\nFinished. All results are in", normalizePath(out), "\n")
