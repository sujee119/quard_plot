# Internal helpers --------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

.qp_version <- "2.0.0"

.stop <- function(...) stop(paste0(...), call. = FALSE)

.warn <- function(...) warning(paste0(...), call. = FALSE)

.msg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(paste0(...))
  invisible(NULL)
}

.fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# Check that `path` is one existing, readable data file. Detects Git LFS
# pointer files (a common reason for cryptic "subscript out of bounds"
# errors when a repository was cloned without git-lfs).
.check_file <- function(path, what = "Input") {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    .stop("`", what, "` must be a single file path (character string).")
  }
  if (!file.exists(path)) .stop(what, " file not found: '", path, "'.")
  if (dir.exists(path)) .stop(what, " path '", path, "' is a directory, not a file.")
  if (file.size(path) == 0) .stop(what, " file '", path, "' is empty.")
  first <- tryCatch(suppressWarnings(readLines(path, n = 1L, warn = FALSE)),
                    error = function(e) "")
  if (length(first) && startsWith(first[1], "version https://git-lfs.github.com/spec/")) {
    .stop(what, " file '", path, "' is a Git LFS pointer, not the data itself. ",
          "Download the real file (for example with `git lfs pull`) or start with ",
          "the small simulated data set: simulate_example_data().")
  }
  invisible(path)
}

.is_compressed <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 6L)
  if (length(magic) < 3L) return(FALSE)
  gz <- magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)
  bz <- rawToChar(magic[1:3]) == "BZh"
  xz <- length(magic) >= 6L && all(magic[1:6] == as.raw(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)))
  gz || bz || xz
}

# Decompress gzip/bzip2/xz input to a temporary file (no extra packages).
.decompress_to_temp <- function(path) {
  tmp <- tempfile(fileext = ".txt")
  inp <- gzfile(path, "rb")
  out <- file(tmp, "wb")
  on.exit({
    close(inp)
    close(out)
  })
  repeat {
    b <- readBin(inp, "raw", 1e7)
    if (!length(b)) break
    writeBin(b, out)
  }
  tmp
}

# data.table::fread for plain or compressed text files, returning a data.frame.
.fread_df <- function(path, ...) {
  src <- path
  if (.is_compressed(path)) {
    src <- .decompress_to_temp(path)
    on.exit(unlink(src), add = TRUE)
  }
  as.data.frame(data.table::fread(src, data.table = FALSE, showProgress = FALSE, ...),
                stringsAsFactors = FALSE)
}

# Find the first column whose simplified name matches one of `candidates`.
.find_col <- function(nms, candidates) {
  simp <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
  low <- simp(nms)
  for (cand in simp(candidates)) {
    hit <- which(low == cand)
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

.check_number <- function(x, what, lower = -Inf, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower) {
    .stop("`", what, "` must be a single number", if (is.finite(lower)) paste0(" >= ", lower) else "", ".")
  }
  invisible(x)
}

# Validate and normalise a region specification.
.as_region <- function(region) {
  if (inherits(region, "qp_region")) return(region)
  if (is.list(region) && all(c("chr", "start", "end") %in% names(region))) {
    out <- list(chr = normalize_chr(region$chr), start = as.numeric(region$start),
                end = as.numeric(region$end), source = region$source %||% "user")
  } else {
    .stop("`region` must be a list with elements chr, start and end ",
          "(for example the output of define_region()).")
  }
  if (!is.finite(out$start) || !is.finite(out$end) || out$end < out$start) {
    .stop("Region end must not be smaller than region start.")
  }
  structure(out, class = "qp_region")
}

#' @export
print.qp_region <- function(x, ...) {
  core <- x$core %||% c(x$start, x$end)
  cat(sprintf("<qp_region> chr %s:%s-%s (%s kb; %s)\n", x$chr, .fmt(round(core[1])),
              .fmt(round(core[2])), format(round((core[2] - core[1]) / 1000, 1), nsmall = 1),
              x$source %||% "user"))
  if (!identical(round(core), round(c(x$start, x$end)))) {
    cat(sprintf("  plot range %s-%s (with margin)\n", .fmt(round(x$start)), .fmt(round(x$end))))
  }
  invisible(x)
}

# Legend placed inside the panel, compatible with ggplot2 < 3.5 and >= 3.5.
.legend_inside <- function(x, y, just = c(1, 1)) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(x, y),
                   legend.justification.inside = just)
  } else {
    ggplot2::theme(legend.position = c(x, y), legend.justification = just)
  }
}

# Create the directory of an output file if needed.
.ensure_dir <- function(file) {
  d <- dirname(file)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(file)
}

# File extension in lower case ("" when there is none).
.file_ext <- function(path) {
  b <- basename(path)
  ifelse(grepl("\\.[A-Za-z0-9]+$", b), tolower(sub("^.*\\.([A-Za-z0-9]+)$", "\\1", b)), "")
}

# Package used to write Excel files, or NA when it is not installed.
.xlsx_engine <- function() {
  if (requireNamespace("writexl", quietly = TRUE)) return("writexl")
  NA_character_
}

.xlsx_help <- paste0("Excel (.xlsx) files need the small 'writexl' package. Install it once with ",
                     "install.packages(\"writexl\"), or use .csv files (they open in Excel too).")

# Write a named list of data frames to one .xlsx file (one sheet each).
.write_xlsx <- function(sheets, file) {
  eng <- .xlsx_engine()
  if (is.na(eng)) .stop(.xlsx_help)
  sheets <- sheets[!vapply(sheets, is.null, NA)]
  nm <- substr(gsub("[\\[\\]:*?/\\\\]", "_", names(sheets), perl = TRUE), 1, 31)
  names(sheets) <- make.unique(nm, sep = "_")
  sheets <- lapply(sheets, function(d) {
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    rownames(d) <- NULL
    d
  })
  .ensure_dir(file)
  writexl::write_xlsx(sheets, file)
  invisible(file)
}

#' Save a table as CSV, tab-separated text or Excel
#'
#' The file type follows the file name: `.csv` (comma-separated; opens in
#' Excel), `.xlsx` (Excel workbook; needs the `writexl` package) or anything
#' else, e.g. `.txt` or `.tsv` (tab-separated text). For `.xlsx`, `x` can
#' also be a named list of tables: each becomes a sheet. If `writexl` is not
#' installed, CSV files are written instead (with a warning that says how to
#' install it).
#'
#' @param x A data frame (or, for `.xlsx`, a named list of data frames).
#' @param file Output file name.
#' @param sheet Sheet name when `x` is a single table written to `.xlsx`.
#' @return The file name(s) written (invisibly).
#' @examples
#' f <- tempfile(fileext = ".csv")
#' write_table(data.frame(gene = c("A", "B"), length = c(1200, 3400)), f)
#' @export
write_table <- function(x, file, sheet = "Sheet1") {
  if (!is.character(file) || length(file) != 1L || !nzchar(file)) .stop("`file` must be a file name.")
  ext <- .file_ext(file)
  .ensure_dir(file)
  if (ext == "xlsx") {
    sheets <- if (is.data.frame(x)) stats::setNames(list(x), sheet) else x
    if (!is.list(sheets) || is.null(names(sheets))) .stop("For .xlsx give a data frame or a named list of data frames.")
    if (is.na(.xlsx_engine())) {
      sheets <- sheets[!vapply(sheets, is.null, NA)]
      base <- sub("\\.xlsx$", "", file, ignore.case = TRUE)
      files <- if (length(sheets) == 1L) paste0(base, ".csv") else paste0(base, "_", names(sheets), ".csv")
      for (k in seq_along(sheets)) data.table::fwrite(sheets[[k]], files[k], na = "", quote = "auto")
      .warn(.xlsx_help, " Written as CSV instead: ", paste(files, collapse = ", "), ".")
      return(invisible(files))
    }
    .write_xlsx(sheets, file)
  } else {
    if (!is.data.frame(x)) .stop("Only .xlsx files can hold several tables: give one data frame.")
    data.table::fwrite(x, file, sep = if (ext == "csv") "," else "\t", na = "", quote = "auto")
  }
  invisible(file)
}
