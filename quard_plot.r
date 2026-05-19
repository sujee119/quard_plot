#' @export
quard_plot <- function(gwas_file, gff_file, vcf_file, block_file_path, 
                       target_snp_chr, target_snp_pos, output_name = "GWAS_Locus_LD_Plot.pdf") {
  
  # --- 1. DATA LOADING & ROBUST PARSING (Restored Original Logic) ---
  raw_gwas <- utils::read.table(gwas_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  df_list <- lapply(strsplit(as.character(raw_gwas[[1]]), ","), function(x) {
    if(length(x) >= 4) {
      data.frame(snp   = x[1], 
                 chrom = as.numeric(gsub("[Cc]hr", "", x[2])),
                 pos   = as.numeric(x[3]),
                 p     = as.numeric(x[4]))
    } else {
      NULL
    }
  })
  
  df_cmplot <- do.call(rbind, df_list)
  
  if(is.null(df_cmplot)) {
    df_cmplot <- data.frame(
      snp   = as.character(raw_gwas$SNP),
      chrom = as.numeric(gsub("[Cc]hr", "", as.character(raw_gwas[[2]]))),
      pos   = as.numeric(raw_gwas[[3]]),
      p     = as.numeric(raw_gwas[[4]])
    )
  }
  
  df_cmplot <- stats::na.omit(df_cmplot)
  
  if(is.null(df_cmplot) || nrow(df_cmplot) == 0) stop("No valid GWAS data found.")
  
  # Identify target SNP name for Panel 1 label
  target_snp_id <- df_cmplot$snp[df_cmplot$chrom == target_snp_chr & df_cmplot$pos == target_snp_pos][1]
  
  gff_all <- utils::read.table(gff_file, sep = "\t", quote = "", stringsAsFactors = FALSE, comment.char = "#")
  block_info <- utils::read.table(block_file_path, header = TRUE)
  
  this_block <- block_info[block_info$CHR == target_snp_chr &
                           block_info$START <= target_snp_pos &
                           block_info$END >= target_snp_pos, ]
  
  if(nrow(this_block) > 0) {
    b_start <- this_block$START[1]
    b_end   <- this_block$END[1]
  } else {
    b_start <- target_snp_pos - 10000
    b_end   <- target_snp_pos + 10000
  }
  
  if((b_end - b_start) < 20000) {
    mid <- (b_start + b_end) / 2
    b_start <- mid - 10000
    b_end   <- mid + 10000
  }
  
  gff_chr_tag <- paste0("Chr", target_snp_chr)
  
  overlapping_genes <- gff_all[gff_all$V1 == gff_chr_tag & gff_all$V3 == "gene" &
                               gff_all$V4 <= b_end & gff_all$V5 >= b_start, ]
  
  if(nrow(overlapping_genes) > 0) {
    block_start_mb <- min(c(overlapping_genes$V4, b_start)) / 1e6
    block_end_mb   <- max(c(overlapping_genes$V5, b_end)) / 1e6
  } else {
    block_start_mb <- b_start / 1e6
    block_end_mb   <- b_end / 1e6
  }


############## USE THIS SECTION TO EDIT FOR OTHER SPECIES ############

  # --- 2. COORDINATES & FILTERING --- You can change this manually for any other species, just add the chromosome lengths in ascending order
  rlengths <- c(43270923, 35937250, 36413819, 35502694, 29958434, 31248787,
                    29697621, 28443022, 23011544, 23207287, 29021106, 27531856)
  chr_offsets <- c(0, cumsum(as.numeric(rlengths))[-12]) # change here for the total chromosome number, if the chromosme number is 23 change it to -23
  names(chr_offsets) <- paste0("Chr", 1:12) # change the chromosome numbers here, if the chromosme number is 23 change it to 1:23
  chr_centers <- chr_offsets + (rlengths / 2)
  
############## dont touch anything below!!!!!!!!!!!!!!!!! ############
  
  df_cmplot$cum_pos <- df_cmplot$pos + chr_offsets[paste0("Chr", df_cmplot$chrom)]
  df_zoom <- df_cmplot[df_cmplot$chrom == target_snp_chr &
                       df_cmplot$pos >= (block_start_mb*1e6) &
                       df_cmplot$pos <= (block_end_mb*1e6), ]
  
  genes_zoom <- gff_all[gff_all$V1 == gff_chr_tag & gff_all$V3 == "gene" &
                        gff_all$V4 <= (block_end_mb*1e6) & gff_all$V5 >= (block_start_mb*1e6), ]
  
  exons_zoom <- NULL
  if(nrow(genes_zoom) > 0) {
    genes_zoom$locus_id <- gsub(".*ID=([^;]+);.*", "\\1", genes_zoom$V9)
    genes_zoom$display_name <- gsub(".*Name=([^;]+);.*", "\\1", genes_zoom$V9)
    exons_zoom <- gff_all[gff_all$V1 == gff_chr_tag & gff_all$V3 == "exon" &
                          gff_all$V4 <= (block_end_mb*1e6) & gff_all$V5 >= (block_start_mb*1e6), ]
    if(nrow(exons_zoom) > 0) {
      exons_zoom$full_parent <- gsub(".*Parent=([^;]+).*", "\\1", exons_zoom$V9)
      exons_zoom$gene_base <- gsub("\\.[0-9]+.*", "", exons_zoom$full_parent)
      exons_zoom$version <- as.numeric(gsub(".*\\.([0-9]+).*", "\\1", exons_zoom$V9))
      exons_zoom <- do.call(rbind, lapply(split(exons_zoom, exons_zoom$gene_base), 
                                          function(sub) sub[sub$version == max(sub$version, na.rm = TRUE), ]))
    }
  }

  # --- 3. LD CALCULATION ---
  ld_mat <- NULL; final_pos <- NULL; snp_names_ld <- NULL; n_snps <- 0
  vcf_data <- data.table::fread(vcf_file, skip = "##", header = TRUE)
  data.table::setnames(vcf_data, 1, "CHROM")
  ld_subset <- vcf_data[CHROM == gff_chr_tag & POS >= floor(block_start_mb*1e6) & POS <= ceiling(block_end_mb*1e6)]
  
  if(nrow(ld_subset) >= 2) {
    sample_cols <- names(ld_subset)[10:ncol(ld_subset)]
    geno_mat <- as.matrix(ld_subset[, ..sample_cols])
    geno_mat[geno_mat %in% c("0/0", "0|0", "0/0/0", ".")] <- "0"
    geno_mat[geno_mat %in% c("0/1", "1/0", "0|1", "1|0")] <- "1"
    geno_mat[geno_mat %in% c("1/1", "1|1")] <- "2"
    geno_mat[!geno_mat %in% c("0","1","2")] <- "0"
    geno_num_mat <- matrix(as.numeric(geno_mat), nrow = nrow(ld_subset))
    vars <- apply(geno_num_mat, 1, stats::var, na.rm = TRUE)
    poly_idx <- which(!is.na(vars) & vars > 1e-8)
    if(length(poly_idx) >= 2) {
      pos_temp <- ld_subset$POS[poly_idx]; sort_order <- order(pos_temp)
      final_pos <- pos_temp[sort_order]
      raw_ld <- stats::cor(t(geno_num_mat[poly_idx, ]), use = "pairwise.complete.obs")^2
      ld_mat <- raw_ld[sort_order, sort_order]; n_snps <- nrow(ld_mat)
      
      # Extract SNP names for the heatmap labels
      snp_names_ld <- df_cmplot$snp[match(final_pos, df_cmplot$pos)]
    }
  }

  # --- 4. FIXED LAYOUT ---
  grDevices::pdf(output_name, width = 12, height = 16)
  graphics::layout(matrix(c(rep(1, 40), 0, 0, rep(2, 36), rep(3, 18), rep(4, 60)), ncol = 1))
  graphics::par(oma = c(5, 1, 1, 0))

  # Panel 1: Manhattan (Custom Axis Step)
  max_p_man <- max(-log10(df_cmplot$p), na.rm = TRUE)
  ylim_manhattan <- c(0, max_p_man * 1.1)
  graphics::par(mar = c(0, 6, 2, 2))
  graphics::plot(df_cmplot$cum_pos, -log10(df_cmplot$p), pch = 19, cex = 0.3, 
                 col = rep(c("grey50", "navy"), 6)[df_cmplot$chrom], xaxt = "n", yaxt = "n", bty = "n", ann = FALSE, ylim = ylim_manhattan) #for arabidopisis col = rep(c("grey50", "navy"), 3), for larger crops, leave it as it is
  
  step_man <- if(max_p_man < 20) 2 else 4
  y_ticks_man <- seq(0, ceiling(ylim_manhattan[2]/step_man)*step_man, by = step_man)
  y_ticks_man <- y_ticks_man[y_ticks_man <= ylim_manhattan[2]]
  graphics::axis(2, at = y_ticks_man, las = 1, lwd = 0.5); graphics::mtext("-log10(p)", side = 2, line = 3.5, cex = 0.8, font = 2)
  graphics::axis(1, at = chr_centers, labels = 1:12, tick = FALSE, cex.axis = 0.8, line = -1)
  graphics::abline(h = 7, col = "red", lty = 2, lwd = 0.8)
  
  zoom_left_cum  <- (block_start_mb * 1e6) + chr_offsets[gff_chr_tag]
  zoom_right_cum <- (block_end_mb * 1e6) + chr_offsets[gff_chr_tag]
  graphics::segments(x0 = c(zoom_left_cum-1e6, zoom_right_cum+1e6), y0 = 0, x1 = c(zoom_left_cum-1e6, zoom_right_cum+1e6), y1 = ylim_manhattan[2], col = "#D3D3D3", lty = "42", lwd = 0.8)
  
  # NEW: Unrotated SNP label near the left line
  if(!is.na(target_snp_id)) {
    graphics::text(zoom_left_cum-1e6, ylim_manhattan[2], labels = target_snp_id, adj = c(1.1, 1), cex = 0.7, font = 2)
  }

  graphics::par(xpd = NA)
  graphics::segments(x0 = c(zoom_left_cum-1e6, zoom_right_cum+1e6), y0 = 0, x1 = c(min(df_cmplot$cum_pos), max(df_cmplot$cum_pos)), y1 = -1, col = "#D3D3D3", lwd = 0.8, lty = "42")
  graphics::par(xpd = FALSE); graphics::box(lwd = 0.5)

  # Panel 2: Zoom (Custom Axis Step)
  max_p_zoom <- max(-log10(df_zoom$p), na.rm = TRUE)
  ylim_zoom <- c(0, max_p_zoom * 1.1)
  graphics::par(mar = c(0, 6, 0, 2))
  graphics::plot(df_zoom$pos/1e6, -log10(df_zoom$p), pch = 19, cex = 1.2, col = "black", xaxt = "n", yaxt = "n", bty = "n", xlim = c(block_start_mb, block_end_mb), ylim = ylim_zoom, ann = FALSE)
  
  step_zoom <- if(max_p_zoom < 20) 2 else 4
  y_ticks_zoom <- seq(0, ceiling(ylim_zoom[2]/step_zoom)*step_zoom, by = step_zoom)
  y_ticks_zoom <- y_ticks_zoom[y_ticks_zoom <= ylim_zoom[2]]
  graphics::axis(2, at = y_ticks_zoom, las = 1, lwd = 0.5); graphics::mtext("-log10(p)", side = 2, line = 3.5, cex = 0.8, font = 2)
  graphics::mtext(paste0("Chromosome ", target_snp_chr), side = 3, line = -1.5, adj = 0.02, font = 2, cex = 0.7)
  graphics::abline(h = 7, col = "red", lty = 2, lwd = 0.8); graphics::box(lwd = 0.5)

  # Panel 3: Gene Track
  graphics::par(mar = c(0, 6, 0, 2))
  graphics::plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlim = c(block_start_mb, block_end_mb), ylim = c(-4, 4), ann = FALSE)
  graphics::axis(2, at = c(1.5, -1.5), labels = c("+", "-"), las = 1, lwd = 0.5); graphics::mtext("Gene track", side = 2, line = 3.5, cex = 0.8, font = 2)
  if(!is.null(genes_zoom) && nrow(genes_zoom) > 0) {
    for(i in 1:nrow(genes_zoom)) {
      is_fwd <- genes_zoom$V7[i] == "+"; is_even <- i %% 2 == 0
      y_mid <- if(is_fwd) ifelse(is_even, 2.6, 1.0) else ifelse(is_even, -2.6, -1.0)
      graphics::lines(c(genes_zoom$V4[i], genes_zoom$V5[i])/1e6, c(y_mid, y_mid), lwd = 1)
      if(!is.null(exons_zoom)) {
        t_exons <- exons_zoom[exons_zoom$gene_base == genes_zoom$locus_id[i], ]
        if(nrow(t_exons) > 0) graphics::rect(t_exons$V4/1e6, y_mid-0.4, t_exons$V5/1e6, y_mid+0.4, col = ifelse(is_fwd, "#2c7bb6", "#d7191c"), border = NA)
      }
      graphics::text((genes_zoom$V4[i]+genes_zoom$V5[i])/(2*1e6), ifelse(is_fwd, y_mid+0.8, y_mid-0.8), labels = genes_zoom$display_name[i], cex = 0.55, font = 1)
    }
  }
  graphics::box(lwd = 0.5)

  # Panel 4: FIXED HEATMAP BOX
  graphics::par(mar = c(0, 6, 0, 2))
  if(!is.null(ld_mat) && n_snps >= 2) {
    plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlim = c(block_start_mb, block_end_mb), ylim = c(0, 110), ann = FALSE, yaxs = "i")
    graphics::mtext(expression(bold(r^2)), side = 2, line = 3.5, cex = 0.8, at = 55)
    plt_ratio <- (graphics::par("pin")[1] / diff(graphics::par("usr")[1:2])) / (graphics::par("pin")[2] / diff(graphics::par("usr")[3:4]))
    
    y_ceiling <- 110; y_hypo <- 105; y_apex <- 5
    idx_x <- seq(block_start_mb, block_end_mb, length.out = n_snps)
    graphics::segments(x0 = final_pos/1e6, y0 = y_ceiling, x1 = idx_x, y1 = y_hypo, lwd = 0.2, col = "grey70")
    
    dx <- (idx_x[2] - idx_x[1]) / 2; dy <- dx * plt_ratio
    col_ramp <- grDevices::colorRampPalette(c("#FFFFE5", "#FED976", "#FD8D3C", "#E31A1C", "#B10026"))(100)
    for(i in 1:(n_snps-1)) for(j in (i+1):n_snps) {
      cx <- (idx_x[i] + idx_x[j]) / 2; cy <- y_hypo - (j - i) * dy
      graphics::polygon(c(cx, cx + dx, cx, cx - dx), c(cy, cy + dy, cy + 2*dy, cy + dy), 
                        col = col_ramp[max(1, round(ld_mat[i,j] * 100))], border = "grey90", lwd = 0.1)
    }
    
    # NEW: Rotated labels on triangle shoulders
    for(k in 1:n_snps) {
      if(!is.na(snp_names_ld[k])) {
        lx <- (idx_x[1] + idx_x[k]) / 2; ly <- y_hypo - (k - 1) * dy
        graphics::text(lx - dx, ly, labels = snp_names_ld[k], srt = 45, adj = c(1, 0.5), cex = 0.4)
        rx <- (idx_x[k] + idx_x[n_snps]) / 2; ry <- y_hypo - (n_snps - k) * dy
        graphics::text(rx + dx, ry, labels = snp_names_ld[k], srt = -45, adj = c(0, 0.5), cex = 0.4)
      }
    }

    leg_x_s <- block_end_mb - (block_end_mb - block_start_mb) * 0.25; leg_x_e <- block_end_mb - (block_end_mb - block_start_mb) * 0.05
    for(k in 1:100) {
      rx_l <- leg_x_s + (k-1) * (leg_x_e - leg_x_s) / 100; rx_r <- leg_x_s + k * (leg_x_e - leg_x_s) / 100
      graphics::rect(rx_l, 8, rx_r, 13, col = col_ramp[k], border = NA)
    }
    graphics::text(c(leg_x_s, (leg_x_s + leg_x_e)/2, leg_x_e), 4, c("0", "0.5", "1"), cex = 0.8, font = 2)
  }
  global_x_at <- pretty(c(block_start_mb, block_end_mb), n = 20)
  graphics::axis(1, at = global_x_at, labels = round(global_x_at, 3), cex.axis = 0.8, lwd = 0.5); graphics::box(lwd = 0.5); grDevices::dev.off()
}