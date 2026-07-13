# PROJECT_ROOT <- Sys.getenv(
#   "VNS_ROOT",
#   unset = normalizePath(here::here(), mustWork = FALSE)
# )

PROJECT_ROOT <- "/proj/gibbons/2026_easton_vns"
OLD_ROOT <- "/proj/gibbons/2024_easton_vns"

paths <- list(
  root         = PROJECT_ROOT,
  raw          = file.path(PROJECT_ROOT, "data", "raw"),
  metaphlan    = file.path(OLD_ROOT, "data", "raw", "metaphlan"),
  intermediate = file.path(PROJECT_ROOT, "data", "intermediate"),
  final_data   = file.path(PROJECT_ROOT, "data", "final"),
  figures      = file.path(PROJECT_ROOT, "results", "figures"),
  tables       = file.path(PROJECT_ROOT, "results", "tables"),
  analysis     = file.path(PROJECT_ROOT, "analysis"),
  funcs        = file.path(PROJECT_ROOT, "R")
)
