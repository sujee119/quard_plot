# Genetic distance and phylogenetic trees -----------------------------------------

#' Genetic distance between samples
#'
#' Pairwise distances from allele dosages, using for each pair only the
#' variants called in both samples.
#'
#' * `"ibs"` (default): proportion of alleles not shared identical-by-state,
#'   sum(|x_i - x_j|) / (2 n), between 0 (identical) and 1.
#' * `"euclidean"`: root mean squared dosage difference per variant.
#'
#' @param geno A `qp_geno` object (apply [genotype_qc()] first).
#' @param method `"ibs"` or `"euclidean"`.
#' @return A symmetric samples x samples distance matrix.
#' @export
genetic_distance <- function(geno, method = c("ibs", "euclidean")) {
  method <- match.arg(method)
  if (!inherits(geno, "qp_geno")) .stop("`geno` must be a qp_geno object.")
  if (nrow(geno$info) < 2L) .stop("At least two variants are needed for genetic distances.")
  X <- t(geno$dosage)
  M <- !is.na(X)
  storage.mode(M) <- "double"
  X0 <- X
  X0[is.na(X0)] <- 0
  N <- tcrossprod(M)
  A <- tcrossprod(X0^2, M)
  sq <- A + t(A) - 2 * tcrossprod(X0)
  sq[sq < 0] <- 0
  if (method == "euclidean") {
    D <- sqrt(sq / N)
  } else {
    ind <- function(v) {
      m <- (X == v)
      m[is.na(m)] <- FALSE
      storage.mode(m) <- "double"
      m
    }
    I0 <- ind(0)
    I2 <- ind(2)
    opp <- tcrossprod(I0, I2)
    D <- (sq - 2 * (opp + t(opp))) / (2 * N)
    D[D < 0] <- 0
  }
  if (any(N == 0)) {
    .warn(sum(N[upper.tri(N)] == 0), " sample pair(s) share no called variant; their distance is set to the maximum.")
    D[N == 0] <- max(D[N > 0], na.rm = TRUE)
  }
  diag(D) <- 0
  dimnames(D) <- list(geno$samples, geno$samples)
  attr(D, "method") <- method
  attr(D, "n_variants") <- nrow(geno$info)
  D
}

# Neighbour-joining (Saitou & Nei 1987) -> unrooted edge list.
.nj_edges <- function(D) {
  n <- nrow(D)
  Dm <- unname(as.matrix(D))
  nodes <- seq_len(n)
  next_id <- n + 1L
  from <- to <- integer()
  len <- numeric()
  while (length(nodes) > 3L) {
    r <- length(nodes)
    S <- rowSums(Dm)
    Q <- (r - 2) * Dm - outer(S, S, "+")
    diag(Q) <- Inf
    k <- which.min(Q)
    i <- (k - 1L) %% r + 1L
    j <- (k - 1L) %/% r + 1L
    dij <- Dm[i, j]
    li <- 0.5 * dij + (S[i] - S[j]) / (2 * (r - 2))
    lj <- dij - li
    u <- next_id
    next_id <- next_id + 1L
    from <- c(from, u, u)
    to <- c(to, nodes[i], nodes[j])
    len <- c(len, max(li, 0), max(lj, 0))
    du <- 0.5 * (Dm[i, ] + Dm[j, ] - dij)
    keep <- setdiff(seq_len(r), c(i, j))
    Dm <- rbind(cbind(Dm[keep, keep, drop = FALSE], du[keep]), c(du[keep], 0))
    nodes <- c(nodes[keep], u)
  }
  if (length(nodes) == 3L) {
    a <- Dm
    u <- next_id
    l <- 0.5 * c(a[1, 2] + a[1, 3] - a[2, 3], a[1, 2] + a[2, 3] - a[1, 3], a[1, 3] + a[2, 3] - a[1, 2])
    from <- c(from, rep(u, 3))
    to <- c(to, nodes)
    len <- c(len, pmax(l, 0))
  } else if (length(nodes) == 2L) {
    from <- c(from, nodes[1])
    to <- c(to, nodes[2])
    len <- c(len, Dm[1, 2])
  }
  data.frame(from = from, to = to, len = len)
}

# Root an unrooted edge list (midpoint or at its last internal node) and
# return an ape-compatible "phylo" object (tips 1..n, root n + 1).
.edges_to_phylo <- function(E, tips, root = c("midpoint", "none")) {
  root <- match.arg(root)
  n <- length(tips)
  key <- function(v) as.character(v)
  nb <- split(c(E$to, E$from), key(c(E$from, E$to)))
  wt <- split(c(E$len, E$len), key(c(E$from, E$to)))
  walk <- function(src) {
    d <- c()
    par <- c()
    d[key(src)] <- 0
    par[key(src)] <- NA
    stack <- src
    while (length(stack)) {
      v <- stack[length(stack)]
      stack <- stack[-length(stack)]
      for (z in seq_along(nb[[key(v)]])) {
        y <- nb[[key(v)]][z]
        if (is.na(d[key(y)])) {
          d[key(y)] <- d[key(v)] + wt[[key(v)]][z]
          par[key(y)] <- v
          stack <- c(stack, y)
        }
      }
    }
    list(d = d, par = par)
  }
  if (root == "midpoint" && n > 2L) {
    w1 <- walk(1L)
    td <- w1$d[key(seq_len(n))]
    a <- seq_len(n)[which.max(td)]
    wa <- walk(a)
    tb <- wa$d[key(seq_len(n))]
    b <- seq_len(n)[which.max(tb)]
    L <- tb[[key(b)]]
    path <- b
    while (path[length(path)] != a) path <- c(path, wa$par[[key(path[length(path)])]])
    dd <- wa$d[key(path)]
    k <- which(dd >= L / 2 & c(dd[-1], -Inf) <= L / 2)[1]
    v <- path[k]
    u <- path[k + 1L]
    r <- max(c(E$from, E$to)) + 1L
    hit <- which((E$from == u & E$to == v) | (E$from == v & E$to == u))[1]
    E <- E[-hit, , drop = FALSE]
    E <- rbind(E, data.frame(from = c(r, r), to = c(u, v), len = c(L / 2 - wa$d[[key(u)]], wa$d[[key(v)]] - L / 2)))
    nb <- split(c(E$to, E$from), key(c(E$from, E$to)))
    wt <- split(c(E$len, E$len), key(c(E$from, E$to)))
    rt <- r
  } else {
    rt <- max(c(E$from, E$to))
  }
  # orient edges away from the root (breadth-first)
  seen <- rt
  queue <- rt
  par_v <- ch_v <- integer()
  len_v <- numeric()
  while (length(queue)) {
    v <- queue[1]
    queue <- queue[-1]
    for (z in seq_along(nb[[key(v)]])) {
      y <- nb[[key(v)]][z]
      if (!(y %in% seen)) {
        seen <- c(seen, y)
        queue <- c(queue, y)
        par_v <- c(par_v, v)
        ch_v <- c(ch_v, y)
        len_v <- c(len_v, wt[[key(v)]][z])
      }
    }
  }
  internal <- unique(c(rt, par_v))
  internal <- internal[internal > n | internal == rt]
  internal <- internal[!(internal %in% seq_len(n))]
  newid <- stats::setNames(n + seq_along(internal), key(internal))
  map <- function(v) ifelse(v <= n, v, newid[key(v)])
  edge <- cbind(as.integer(map(par_v)), as.integer(map(ch_v)))
  structure(list(edge = edge, edge.length = pmax(len_v, 0), tip.label = tips, Nnode = length(internal)),
            class = "phylo")
}

.hclust_to_phylo <- function(hc) {
  n <- length(hc$order)
  m <- hc$merge
  h <- hc$height / 2
  node_id <- function(k) n + (n - k)
  parent <- child <- integer()
  len <- numeric()
  for (k in seq_len(n - 1L)) {
    for (x in m[k, ]) {
      cid <- if (x < 0) -x else node_id(x)
      ch <- if (x < 0) 0 else h[x]
      parent <- c(parent, node_id(k))
      child <- c(child, cid)
      len <- c(len, h[k] - ch)
    }
  }
  structure(list(edge = cbind(as.integer(parent), as.integer(child)), edge.length = len,
                 tip.label = hc$labels, Nnode = n - 1L), class = "phylo")
}

#' Build a tree from a distance matrix
#'
#' Neighbour-joining (Saitou & Nei 1987; uses the faster `ape::nj()` for
#' large data sets when the ape package is installed) or UPGMA (average
#' linkage). NJ trees are midpoint-rooted for display by default. The result
#' is an ape-compatible `phylo` object, so it can also be used with ape,
#' ggtree or iTOL (via [write_newick()]).
#'
#' @param d Distance matrix (e.g. from [genetic_distance()]) or `dist` object.
#' @param method `"nj"` or `"upgma"`.
#' @param root `"midpoint"` or `"none"` (NJ only).
#' @return A `phylo` object (attribute `distance` keeps the distance matrix).
#' @export
build_tree <- function(d, method = c("nj", "upgma"), root = c("midpoint", "none")) {
  method <- match.arg(method)
  root <- match.arg(root)
  D <- as.matrix(d)
  n <- nrow(D)
  if (n < 3L) .stop("At least three samples are needed for a tree.")
  tips <- rownames(D) %||% paste0("S", seq_len(n))
  if (method == "upgma") {
    hc <- stats::hclust(stats::as.dist(D), method = "average")
    hc$labels <- tips
    tr <- .hclust_to_phylo(hc)
  } else {
    if (n > 400L && requireNamespace("ape", quietly = TRUE)) {
      a <- get("nj", envir = asNamespace("ape"))(stats::as.dist(D))
      E <- data.frame(from = a$edge[, 1], to = a$edge[, 2], len = a$edge.length)
      tr <- .edges_to_phylo(E, a$tip.label, root)
    } else {
      if (n > 800L) message("Neighbour-joining ", n, " samples without the 'ape' package may take a while.")
      tr <- .edges_to_phylo(.nj_edges(D), tips, root)
    }
  }
  attr(tr, "distance") <- D
  attr(tr, "method") <- method
  tr
}

#' Write a tree in Newick format
#'
#' @param tree A `phylo` object from [build_tree()].
#' @param file Optional output file (`.nwk`); the Newick string is returned
#'   invisibly in any case. Opens in FigTree, iTOL, MEGA, ape, ggtree.
#' @param digits Significant digits of branch lengths.
#' @export
write_newick <- function(tree, file = NULL, digits = 6) {
  E <- tree$edge
  n <- length(tree$tip.label)
  len <- tree$edge.length %||% rep(1, nrow(E))
  root <- setdiff(E[, 1], E[, 2])[1]
  esc <- function(s) ifelse(grepl("[][ (),:;']", s), paste0("'", gsub("'", "''", s), "'"), s)
  ord <- root
  i <- 1L
  while (i <= length(ord)) {
    ord <- c(ord, E[E[, 1] == ord[i], 2])
    i <- i + 1L
  }
  str <- character(max(E))
  for (v in rev(ord)) {
    if (v <= n) {
      str[v] <- esc(tree$tip.label[v])
    } else {
      ek <- which(E[, 1] == v)
      str[v] <- paste0("(", paste0(str[E[ek, 2]], ":", signif(len[ek], digits), collapse = ","), ")")
    }
  }
  out <- paste0(str[root], ";")
  if (!is.null(file)) {
    .ensure_dir(file)
    writeLines(out, file)
  }
  invisible(out)
}

# Node coordinates for drawing.
.tree_coords <- function(tree, ladderize = TRUE) {
  E <- tree$edge
  n <- length(tree$tip.label)
  len <- tree$edge.length %||% rep(1, nrow(E))
  N <- max(E)
  root <- setdiff(E[, 1], E[, 2])[1]
  kids <- split(E[, 2], factor(E[, 1], levels = seq_len(N)))
  klen <- split(len, factor(E[, 1], levels = seq_len(N)))
  bfs <- root
  i <- 1L
  while (i <= length(bfs)) {
    bfs <- c(bfs, kids[[bfs[i]]])
    i <- i + 1L
  }
  ntip <- integer(N)
  ntip[seq_len(n)] <- 1L
  for (v in rev(bfs)) if (v > n) ntip[v] <- sum(ntip[kids[[v]]])
  if (ladderize) {
    for (v in bfs[bfs > n]) {
      o <- order(ntip[kids[[v]]])
      kids[[v]] <- kids[[v]][o]
      klen[[v]] <- klen[[v]][o]
    }
  }
  depth <- numeric(N)
  for (v in bfs) {
    if (v > n || v == root) {
      kv <- kids[[v]]
      if (length(kv)) depth[kv] <- depth[v] + klen[[v]]
    }
  }
  # tip order by depth-first traversal
  y <- numeric(N)
  stack <- root
  tip_rank <- 0L
  while (length(stack)) {
    v <- stack[length(stack)]
    stack <- stack[-length(stack)]
    if (v <= n) {
      tip_rank <- tip_rank + 1L
      y[v] <- tip_rank
    } else {
      stack <- c(stack, rev(kids[[v]]))
    }
  }
  for (v in rev(bfs)) if (v > n) y[v] <- mean(range(y[kids[[v]]]))
  list(E = E, len = len, n = n, N = N, root = root, kids = kids, klen = klen, ntip = ntip,
       depth = depth, y = y, bfs = bfs)
}

.group_colours <- function(lv, colors = NULL) {
  if (!is.null(colors)) {
    if (!is.null(names(colors))) return(colors[lv])
    return(stats::setNames(rep_len(colors, length(lv)), lv))
  }
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F3C8D", "#B8860B")
  pal <- if (length(lv) <= length(base)) base[seq_along(lv)] else grDevices::hcl.colors(length(lv), "Dark 3")
  stats::setNames(pal, lv)
}

.as_groups <- function(groups, tips) {
  if (is.null(groups)) return(NULL)
  if (is.character(groups) && length(groups) == 1L && file.exists(groups)) {
    g <- .fread_df(groups, header = TRUE, colClasses = "character")
    groups <- stats::setNames(as.character(g[[2]]), as.character(g[[1]]))
  } else if (is.data.frame(groups)) {
    groups <- stats::setNames(as.character(groups[[2]]), as.character(groups[[1]]))
  } else if (is.factor(groups)) {
    groups <- stats::setNames(as.character(groups), names(groups))
  }
  if (is.null(names(groups))) {
    if (length(groups) != length(tips)) .stop("Unnamed `groups` must have one value per sample.")
    names(groups) <- tips
  }
  out <- unname(groups[tips])
  if (all(is.na(out))) .warn("No sample names in `groups` match the tree tips.")
  out
}

#' Plot a phylogenetic tree
#'
#' Draws a tree from [build_tree()] in a circular (fan), rectangular or
#' unrooted layout, similar to the GAPIT phylogenetic tree output, with the
#' name of every individual (sample) at its leaf. Tips and branches can be
#' coloured by known groups (e.g. subpopulations) or by `k` clusters cut
#' from the distance matrix.
#'
#' The names stay readable: they are written at `label_size` pt and the
#' recommended PDF size grows with the number of individuals (stored in
#' `attr(p, "size_mm")` and used by [save_pdf()] and [hapmap_tree()] when no
#' size is given). If you fix the page size with `size_mm`, the names are
#' made as large as fits on that page.
#'
#' @param tree A `phylo` object from [build_tree()].
#' @param layout `"circular"`, `"rectangular"` or `"unrooted"`.
#' @param groups Optional groups: a named vector (names = sample names), a
#'   two-column data frame (sample, group) or a CSV/TSV file with those two
#'   columns.
#' @param k Optional number of clusters (average-linkage clustering of the
#'   tree's distance matrix) used for colours when `groups` is not given.
#' @param labels Write the sample names at the leaves (default `TRUE`).
#' @param label_size Font size of the sample names in pt (default 6, or what
#'   fits on the page when `size_mm` is given).
#' @param label_color Colour of the names, or `"group"` to use the group colours.
#' @param align_labels Write all names in one line (rectangular) or circle
#'   (circular), joined to their leaf by a dotted line, so that they never
#'   overlap. `FALSE` writes each name next to its leaf.
#' @param tip_points Draw a point at each tip.
#' @param colors Optional group colours.
#' @param branch_color Colour of branches not assigned to a single group.
#' @param branch_width Branch line width.
#' @param scale_bar Add a scale bar for branch lengths.
#' @param title Optional title.
#' @param size_mm Optional page size `c(width, height)` in mm the plot will be
#'   saved at (`NA` = automatic); the names are fitted to it.
#' @param base_size Base font size (pt) of the title and legend.
#' @return A ggplot object with attribute `size_mm` (recommended PDF size).
#'   Save with [save_pdf()].
#' @export
plot_tree <- function(tree, layout = c("circular", "rectangular", "unrooted"), groups = NULL, k = NULL,
                      labels = TRUE, label_size = NULL, label_color = "grey15", align_labels = TRUE,
                      tip_points = TRUE, colors = NULL, branch_color = "grey35", branch_width = 0.3,
                      scale_bar = TRUE, title = NULL, size_mm = NULL, base_size = 8) {
  layout <- match.arg(layout)
  if (!inherits(tree, "phylo")) .stop("`tree` must be a phylo object (see build_tree()).")
  tc <- .tree_coords(tree)
  n <- tc$n
  N <- tc$N
  labels <- !isFALSE(labels)
  grp <- .as_groups(groups, tree$tip.label)
  if (is.null(grp) && !is.null(k)) {
    D <- attr(tree, "distance")
    if (is.null(D)) .stop("`k` needs the distance matrix stored by build_tree().")
    cl <- stats::cutree(stats::hclust(stats::as.dist(D), method = "average"), k = k)
    grp <- paste("Cluster", cl[tree$tip.label])
  }
  node_grp <- rep(NA_character_, N)
  if (!is.null(grp)) {
    node_grp[seq_len(n)] <- grp
    for (v in rev(tc$bfs)) {
      if (v > n) {
        g <- unique(node_grp[tc$kids[[v]]])
        node_grp[v] <- if (length(g) == 1L) g else NA_character_
      }
    }
  }
  lv <- if (!is.null(grp)) sort(unique(stats::na.omit(grp))) else character()
  pal <- c(.group_colours(lv, colors), other = branch_color)
  E <- tc$E
  ecol <- ifelse(is.na(node_grp[E[, 2]]), "other", node_grp[E[, 2]])
  tip_col <- ifelse(is.na(node_grp[seq_len(n)]), "other", node_grp[seq_len(n)])
  maxd <- max(tc$depth)
  if (!is.finite(maxd) || maxd <= 0) maxd <- 1

  # page geometry: names of `fs` pt, 1 pt = 0.3528 mm --------------------------
  pt <- 0.3528
  nch <- max(1L, nchar(tree$tip.label, type = "width"))
  cw <- 0.62                         # average character width / font size
  gap <- 1.2                         # mm between leaf and name
  pad <- 3                           # mm around the tree
  extra <- 24                        # mm for title, legend and scale bar
  page_max <- 5000                   # mm (PDF page size limit)
  sz <- if (is.null(size_mm)) c(NA_real_, NA_real_) else rep(as.numeric(size_mm), length.out = 2L)
  lab_mm <- function(fs) if (labels) gap + cw * nch * fs * pt else 0
  if (layout == "rectangular") {
    W <- if (is.na(sz[1])) 180 else sz[1]
    fs <- label_size %||% if (is.na(sz[2])) 6 else min(8, max(1, (sz[2] - extra) / (n * 1.3 * pt)))
    H <- if (is.na(sz[2])) max(90, n * 1.3 * fs * pt + extra) else sz[2]
    if (H > page_max) {
      fs <- max(1, (page_max - extra) / (n * 1.3 * pt))
      H <- page_max
    }
  } else {
    need_r <- function(fs) max(40, 1.3 * fs * pt * n / (2 * pi))
    if (all(is.na(sz))) {
      fs <- label_size %||% 6
      half <- need_r(fs) + lab_mm(fs) + pad
      if (2 * half > page_max) {
        fs <- max(1, (page_max / 2 - pad - gap) / (pt * (1.3 * n / (2 * pi) + cw * nch)))
        half <- need_r(fs) + lab_mm(fs) + pad
      }
    } else {
      half <- min(c(sz[1] - 4, sz[2] - extra), na.rm = TRUE) / 2
      fs <- label_size %||% min(8, max(1, (half - pad - gap) / (pt * (1.3 * n / (2 * pi) + cw * nch))))
    }
    R <- max(10, half - lab_mm(fs) - pad)
    W <- if (is.na(sz[1])) 2 * half + 4 else sz[1]
    H <- if (is.na(sz[2])) 2 * half + extra else sz[2]
  }
  tsize <- fs / ggplot2::.pt

  conn <- NULL
  if (layout == "rectangular") {
    seg <- data.frame(x = tc$depth[E[, 1]], xend = tc$depth[E[, 2]], y = tc$y[E[, 2]], yend = tc$y[E[, 2]],
                      col = ecol)
    vert <- do.call(rbind, lapply(tc$bfs[tc$bfs > n], function(v) {
      yy <- tc$y[tc$kids[[v]]]
      data.frame(x = tc$depth[v], xend = tc$depth[v], y = min(yy), yend = max(yy),
                 col = ifelse(is.na(node_grp[v]), "other", node_grp[v]))
    }))
    seg <- rbind(seg, vert)
    tx <- tc$depth[seq_len(n)]
    ty <- tc$y[seq_len(n)]
    tree_w <- max(20, W - 6 - lab_mm(fs) - 2)          # mm for the branches
    sx <- maxd / tree_w                                   # data units per mm
    lx <- (if (labels && align_labels) maxd else tx) + gap * sx
    tips <- data.frame(x = tx, y = ty, lab = tree$tip.label, col = tip_col, lx = lx, ly = ty,
                       angle = 0, hjust = 0)
    if (labels && align_labels) {
      far <- tx < maxd - 1e-9 * maxd
      if (any(far)) conn <- data.frame(x = tx[far], xend = maxd, y = ty[far], yend = ty[far])
    }
  } else if (layout == "circular") {
    th <- function(y) 2 * pi * (y - 1) / n
    radial <- data.frame(r0 = tc$depth[E[, 1]], r1 = tc$depth[E[, 2]], a = th(tc$y[E[, 2]]), col = ecol)
    seg <- data.frame(x = radial$r0 * cos(radial$a), y = radial$r0 * sin(radial$a),
                      xend = radial$r1 * cos(radial$a), yend = radial$r1 * sin(radial$a), col = radial$col)
    arcs <- do.call(rbind, lapply(tc$bfs[tc$bfs > n], function(v) {
      aa <- th(range(tc$y[tc$kids[[v]]]))
      m <- max(2L, ceiling((aa[2] - aa[1]) / (pi / 180)))
      a <- seq(aa[1], aa[2], length.out = m + 1L)
      r <- tc$depth[v]
      data.frame(x = r * cos(a[-length(a)]), y = r * sin(a[-length(a)]), xend = r * cos(a[-1]),
                 yend = r * sin(a[-1]), col = ifelse(is.na(node_grp[v]), "other", node_grp[v]))
    }))
    seg <- rbind(seg, arcs)
    s_mm <- maxd / R                                      # data units per mm
    a <- th(tc$y[seq_len(n)])
    deg <- a * 180 / pi
    flip <- deg > 90 & deg < 270
    r_tip <- tc$depth[seq_len(n)]
    r_lab <- (if (labels && align_labels) rep(maxd, n) else r_tip) + gap * s_mm
    tips <- data.frame(x = r_tip * cos(a), y = r_tip * sin(a), lab = tree$tip.label, col = tip_col,
                       lx = r_lab * cos(a), ly = r_lab * sin(a),
                       angle = ifelse(flip, deg + 180, deg), hjust = ifelse(flip, 1, 0))
    if (labels && align_labels) {
      far <- r_tip < maxd - 1e-9 * maxd
      if (any(far)) {
        conn <- data.frame(x = r_tip[far] * cos(a[far]), y = r_tip[far] * sin(a[far]),
                           xend = maxd * cos(a[far]), yend = maxd * sin(a[far]))
      }
    }
  } else {
    pos <- matrix(0, N, 2)
    wedge <- matrix(0, N, 2)
    wedge[tc$root, ] <- c(0, 2 * pi)
    for (v in tc$bfs) {
      kv <- tc$kids[[v]]
      if (!length(kv)) next
      a0 <- wedge[v, 1]
      span <- wedge[v, 2] - wedge[v, 1]
      sizes <- tc$ntip[kv] / sum(tc$ntip[kv])
      ends <- a0 + cumsum(sizes) * span
      starts <- c(a0, ends[-length(ends)])
      for (z in seq_along(kv)) {
        wedge[kv[z], ] <- c(starts[z], ends[z])
        phi <- (starts[z] + ends[z]) / 2
        pos[kv[z], ] <- pos[v, ] + tc$klen[[v]][z] * c(cos(phi), sin(phi))
      }
    }
    seg <- data.frame(x = pos[E[, 1], 1], y = pos[E[, 1], 2], xend = pos[E[, 2], 1], yend = pos[E[, 2], 2],
                      col = ecol)
    rmax <- max(sqrt(rowSums(pos^2)), 1e-12)
    s_mm <- rmax / R
    phi <- (wedge[seq_len(n), 1] + wedge[seq_len(n), 2]) / 2
    deg <- (phi * 180 / pi) %% 360
    flip <- deg > 90 & deg < 270
    tips <- data.frame(x = pos[seq_len(n), 1], y = pos[seq_len(n), 2], lab = tree$tip.label, col = tip_col,
                       angle = ifelse(flip, deg + 180, deg), hjust = ifelse(flip, 1, 0))
    tips$lx <- tips$x + gap * s_mm * cos(phi)
    tips$ly <- tips$y + gap * s_mm * sin(phi)
  }

  p <- ggplot2::ggplot()
  if (!is.null(conn)) {
    p <- p + ggplot2::geom_segment(data = conn, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                                                             yend = .data$yend),
                                   colour = "grey75", linewidth = 0.12, linetype = "dotted")
  }
  p <- p + ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                                                          yend = .data$yend, colour = .data$col),
                                 linewidth = branch_width, lineend = "round")
  if (tip_points) {
    p <- p + ggplot2::geom_point(data = tips, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$col),
                                 size = max(0.25, min(1.6, 0.9 * fs * pt * 2)))
  }
  if (labels) {
    if (identical(label_color, "group")) {
      p <- p + ggplot2::geom_text(data = tips, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$lab,
                                                            angle = .data$angle, hjust = .data$hjust,
                                                            colour = .data$col),
                                  size = tsize, vjust = 0.5, show.legend = FALSE)
    } else {
      p <- p + ggplot2::geom_text(data = tips, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$lab,
                                                            angle = .data$angle, hjust = .data$hjust),
                                  colour = label_color, size = tsize, vjust = 0.5)
    }
  }
  brk <- lv
  p <- p + ggplot2::scale_colour_manual(values = pal, breaks = brk, name = NULL,
                                        guide = if (length(brk)) "legend" else "none")
  sb <- signif(maxd / 5, 1)
  if (layout == "rectangular") {
    xmax <- maxd + (lab_mm(fs) + 2) * sx
    ylo <- if (scale_bar) -2.5 else 0.3
    if (scale_bar && maxd > 0) {
      p <- p + ggplot2::annotate("segment", x = 0, xend = sb, y = -1, yend = -1, linewidth = 0.4) +
        ggplot2::annotate("text", x = sb / 2, y = -1, label = format(sb), vjust = 1.5, size = .txt(base_size, 0.75))
    }
    p <- p + ggplot2::scale_x_continuous(limits = c(-maxd * 0.01, xmax), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(ylo, n + 0.7), expand = c(0, 0)) +
      ggplot2::coord_cartesian(clip = "off")
  } else {
    lim <- half * (if (layout == "circular") maxd / R else rmax / R)
    if (scale_bar && maxd > 0) {
      x0 <- -lim * 0.98
      y0 <- -lim * 0.97
      p <- p + ggplot2::annotate("segment", x = x0, xend = x0 + sb, y = y0, yend = y0, linewidth = 0.4) +
        ggplot2::annotate("text", x = x0 + sb / 2, y = y0, label = format(sb), vjust = -0.6,
                          size = .txt(base_size, 0.75))
    }
    p <- p + ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE, clip = "off")
  }
  p <- p + ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = base_size),
                   plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
                   plot.margin = ggplot2::margin(4, 4, 4, 4)) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.5, linewidth = 1)))
  ttl <- title %||% sprintf("%s tree of %s samples (%s distance%s)",
                            if (identical(attr(tree, "method"), "upgma")) "UPGMA" else "Neighbour-joining",
                            .fmt(n), toupper(attr(attr(tree, "distance"), "method") %||% "genetic"),
                            if (!is.null(attr(attr(tree, "distance"), "n_variants")))
                              paste0(", ", .fmt(attr(attr(tree, "distance"), "n_variants")), " variants") else "")
  p <- p + ggplot2::ggtitle(ttl)
  attr(p, "groups") <- if (!is.null(grp)) stats::setNames(grp, tree$tip.label) else NULL
  attr(p, "size_mm") <- round(c(W, H))
  attr(p, "label_size") <- fs
  p
}

#' Phylogenetic tree from a HapMap or VCF file
#'
#' Whole-genome wrapper similar to the GAPIT phylogenetic tree: reads an
#' evenly spaced subset of variants (`max_snps`), applies genotype QC,
#' computes genetic distances, builds a neighbour-joining (or UPGMA) tree,
#' and writes the plot to PDF and the tree to a Newick file (`.nwk`) that
#' opens in FigTree, iTOL or MEGA.
#'
#' @param file HapMap or VCF file (or a `qp_geno` object).
#' @param output PDF file for the tree.
#' @param max_snps Maximum number of variants used (evenly spaced along the
#'   file); 5,000-20,000 is usually plenty for relationships.
#' @param min_maf,max_missing Variant filters (see [genotype_qc()]).
#' @param distance `"ibs"` or `"euclidean"` (see [genetic_distance()]).
#' @param method `"nj"` or `"upgma"` (see [build_tree()]).
#' @param layout `"circular"`, `"rectangular"` or `"unrooted"`.
#' @param groups,k Group colours (see [plot_tree()]).
#' @param newick Also write `<output>.nwk`.
#' @param width,height PDF size in mm. Default: chosen from the number of
#'   individuals so that every name is readable (see [plot_tree()]); a fixed
#'   size makes the names as large as fits.
#' @param verbose Print progress.
#' @param ... Further arguments for [plot_tree()].
#' @return Invisibly, a list with `tree`, `distance`, `plot` and `files`.
#' @export
hapmap_tree <- function(file, output = "tree.pdf", max_snps = 10000, min_maf = 0.05, max_missing = 0.2,
                        distance = c("ibs", "euclidean"), method = c("nj", "upgma"),
                        layout = c("circular", "rectangular", "unrooted"), groups = NULL, k = NULL,
                        newick = TRUE, width = NULL, height = NULL, verbose = TRUE, ...) {
  distance <- match.arg(distance)
  method <- match.arg(method)
  layout <- match.arg(layout)
  if (!grepl("\\.pdf$", output, ignore.case = TRUE)) .stop("`output` must end in .pdf.")
  t0 <- proc.time()[["elapsed"]]
  if (inherits(file, "qp_geno")) {
    g <- file
  } else {
    n_lines <- .count_data_lines(file)
    every <- max(1L, floor(n_lines / max_snps))
    .msg(verbose, "Reading every ", every, " of ", .fmt(n_lines), " variant(s) (max_snps = ", .fmt(max_snps), ").")
    g <- read_genotypes(file, every = every, verbose = verbose)
  }
  g <- genotype_qc(g, min_maf = min_maf, max_missing = max_missing, snps_only = FALSE, verbose = verbose)
  empty <- colSums(!is.na(g$dosage)) == 0
  if (any(empty)) {
    .warn(sum(empty), " sample(s) have no called genotype and were left out of the tree: ",
          paste(utils::head(g$samples[empty], 10), collapse = ", "), if (sum(empty) > 10) ", ..." else "", ".")
    g <- .geno_keep_samples(g, !empty)
  }
  D <- genetic_distance(g, distance)
  tr <- build_tree(D, method)
  fixed <- if (is.null(width) && is.null(height)) NULL else c(width %||% NA_real_, height %||% NA_real_)
  p <- plot_tree(tr, layout = layout, groups = groups, k = k, size_mm = fixed, ...)
  save_pdf(p, output, width = width %||% attr(p, "size_mm")[1], height = height %||% attr(p, "size_mm")[2])
  files <- c(pdf = output)
  if (newick) {
    nw <- sub("\\.pdf$", ".nwk", output, ignore.case = TRUE)
    write_newick(tr, nw)
    files <- c(files, newick = nw)
  }
  .msg(verbose, "Tree of ", .fmt(length(tr$tip.label)), " samples from ", .fmt(nrow(g$info)), " variants written to ",
       paste(files, collapse = " and "), " in ", round(proc.time()[["elapsed"]] - t0, 1), " s.")
  invisible(list(tree = tr, distance = D, plot = p, files = files))
}

# Number of data lines in a VCF or HapMap file.
.count_data_lines <- function(file) {
  con <- file(file, "r")
  on.exit(close(con))
  n <- 0
  repeat {
    x <- readLines(con, n = 50000L, warn = FALSE)
    if (!length(x)) break
    n <- n + sum(nzchar(x) & !startsWith(x, "#"))
  }
  if (.is_hapmap_file(file)) n <- n - 1
  n
}
