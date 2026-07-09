## -----------------------------------------------------------------------------
##
## Script Name: R/utils.R
##
## Purpose of Script: General utility functions for project management, 
##    specifically focused on maintaining data lineage and naming consistency 
##    across the bioinformatics pipeline.
##
## Author: Alyssa N. Easton
##
## Date Created: 2026-05-13
##
## Copyright (c) Alyssa N. Easton, 2026
##
## -----------------------------------------------------------------------------
##
## Notes:
##   - Contains get_prefix() for dynamic file naming based on the calling script.
##   - Designed to be sourced at the beginning of RMarkdown notebooks.
##
## -----------------------------------------------------------------------------

# Load required packages
library(knitr)  # Required for get_prefix() logic

#' Generate a short filename prefix based on the current RMarkdown file
#'
#' By default, this extracts the leading numeric prefix from the current
#' notebook filename. For example:
#'   "06_A_regress_mediate.Rmd" -> "06_"
#'   "11_tables.Rmd"            -> "11_"
#'
#' If the filename does not begin with digits, the full filename stem is used:
#'   "cytokine_plot_final.Rmd" -> "cytokine_plot_final_"
#'
#' Works when knitting and, when possible, in interactive RStudio sessions.
#'
#' @param fallback Prefix to use if no current file can be detected.
#' @param numeric_only If TRUE, use only the leading numeric prefix when present.
#'
#' @return A string prefix ending in "_".
get_prefix <- function(fallback = "interactive_", numeric_only = TRUE) {
  
  # 1. Check if the file is currently being knit/rendered
  current_file <- knitr::current_input()
  
  # 2. If not knitting, try to get the path from the RStudio editor
  if (is.null(current_file) || current_file == "") {
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      
      current_path <- rstudioapi::getSourceEditorContext()$path
      
      if (!is.null(current_path) && current_path != "") {
        current_file <- basename(current_path)
      }
    }
  }
  
  # 3. Final fallback if file is unsaved or not in RStudio
  if (is.null(current_file) || current_file == "") {
    return(fallback)
  }
  
  # 4. Remove extension to create filename stem
  file_stem <- tools::file_path_sans_ext(basename(current_file))
  
  # 5. Extract leading numeric prefix, if requested and available
  if (numeric_only) {
    numeric_prefix <- regmatches(file_stem, regexpr("^[0-9]+", file_stem))
    
    if (length(numeric_prefix) == 1 && numeric_prefix != "") {
      return(paste0(numeric_prefix, "_"))
    }
  }
  
  # 6. Otherwise use the full file stem
  paste0(file_stem, "_")
}


#' Standardized save function that enforces script-based naming
#' @param object The R object to save
#' @param base_name Descriptive name for the file (e.g., "cleaned_meta")
#' @param dir Directory to save in (relative to project root)
#' @param ext File extension (default is ".rds")
save_output <- function(object, base_name, dir, ext = ".rds") {
  # Get dynamic prefix from the current notebook
  prefix <- get_prefix()
  
  # Construct full file path
  full_name <- paste0(prefix, base_name, ext)
  full_path <- file.path(dir, full_name)
  
  # Ensure directory exists before saving
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  message(paste("Saving object to:", full_path))
  
  if (ext == ".rds") {
    saveRDS(object, full_path)
  } else if (ext == ".csv") {
    write.csv(object, full_path, row.names = FALSE)
  } else {
    stop("Unsupported extension. Use .rds or .csv.")
  }
}

# ensure_dir
# read_intermediate
# save_intermediate