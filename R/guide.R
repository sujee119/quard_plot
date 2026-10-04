# Step-by-step guide ------------------------------------------------------------

#' Copy the step-by-step guide to your folder
#'
#' Copies the guide `quard_plot_examples.R` - every function of quardplot,
#' explained step by step in plain language, with the structure your input
#' files must have - into `dir` (default: the working directory) and opens
#' it, so you can run it line by line with your own files.
#'
#' @param dir Folder to copy the guide to.
#' @param open Open the copy (in RStudio or your default editor).
#' @param overwrite Replace an existing copy with the same name.
#' @return The path of the copy (invisibly).
#' @examples
#' \dontrun{
#' quard_guide()            # copies quard_plot_examples.R to the working directory
#' }
#' @export
quard_guide <- function(dir = ".", open = interactive(), overwrite = FALSE) {
  src <- system.file("guide", "quard_plot_examples.R", package = "quardplot")
  if (!nzchar(src)) .stop("The guide was not found in the installed quardplot package.")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dest <- file.path(dir, "quard_plot_examples.R")
  if (file.exists(dest) && !overwrite) {
    message("'", dest, "' already exists and was not replaced (use overwrite = TRUE to replace it).")
  } else {
    file.copy(src, dest, overwrite = TRUE)
    message("The guide was copied to ", normalizePath(dest), ".")
  }
  if (isTRUE(open)) utils::file.edit(dest)
  invisible(dest)
}
