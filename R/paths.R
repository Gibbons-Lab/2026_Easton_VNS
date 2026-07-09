PROJECT_ROOT <- Sys.getenv(
  "VNS_ROOT",
  unset = normalizePath(here::here(), mustWork = FALSE)
)

paths <- list(
  root         = PROJECT_ROOT,
  raw          = file.path(PROJECT_ROOT, "data", "raw"),
  metadata     = file.path(PROJECT_ROOT, "data", "raw", "metadata"),
  metaphlan    = file.path(PROJECT_ROOT, "data", "raw", "metaphlan"),
  intermediate = file.path(PROJECT_ROOT, "data", "intermediate"),
  final_data   = file.path(PROJECT_ROOT, "data", "final"),
  figures      = file.path(PROJECT_ROOT, "results", "figures"),
  tables       = file.path(PROJECT_ROOT, "results", "tables")
)