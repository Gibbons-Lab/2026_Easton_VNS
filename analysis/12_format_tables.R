# table_format.R
setwd(file.path(Sys.getenv("VNS_ROOT"), "data"))  # stored VNS_ROOT in ~/.Renviron
options(expressions = 6000)  # Increase allowed recursion depth
options(digits = 10) ## Formats output to 10 digits
library(magrittr)
library(tidyverse)
library(moments)
library(broom)
library(broom.helpers)
library(HDMT)
library(rlang)
library(glue)
library(janitor)
library(openxlsx)
library(purrr)


# open up data to make tables
data <- readRDS("/proj/gibbons/2024_easton_vns/data/1_preprocessed/species_october.rds")
lm2a <- readRDS("/proj/gibbons/2024_easton_vns/data/lmer/A_microbe_2var.rds")
lm3a <- readRDS("/proj/gibbons/2024_easton_vns/data/lmer/A_microbe_3var.rds")
lmer2b <- readRDS("/proj/gibbons/2024_easton_vns/data/lmer/B_microbe_2var.rds")
lmer3b <- readRDS("/proj/gibbons/2024_easton_vns/data/lmer/B_microbe_3var.rds")

taxmap_to_merge <- data$taxonomy_map %>%
  dplyr::select(c("Species", "Genus", "Family", "Order", "Class", "Phylum", "full_name")) %>%
  rename("y" = "Species")

results <- list(
  "lm2a" = lm2a$lmer_res,
  "lm3a" = lm3a$lmer_res,
  "lmer2b" = lmer2b$lmer_res,
  "lmer3b" = lmer3b$lmer_res
)

for (res in names(results)) {
  results[[res]]$y <- gsub("s__", "", results[[res]]$y)     # cut off "s__" from taxon names
  results[[res]]$term <- gsub("\\(Intercept)", "Intercept", results[[res]]$term)  # get rid of parentheses
  results[[res]]$term <- gsub("blast_and_vns", "interaction", results[[res]]$term)  # get rid of parentheses
  
  results[[res]] <- results[[res]] %>%
    dplyr::select(-signif) %>%
    pivot_wider(
      names_from = "term",
      values_from = c("estimate", "se", "p", "fdr"),
      names_glue = "{term}_{.value}"
    ) %>% 
    dplyr::select(
      any_of(c("y", "r2", "r2_marginal", "r2_conditional")),
      any_of(matches("^Intercept_")),
      matches("^blast_"),
      matches("^vns_"),
      matches("^interaction_"),
      everything()
    ) %>%
    left_join(taxmap_to_merge, by = "y")
}



# save everything to excel spreaadsheets
wb <- createWorkbook()

addWorksheet(wb, "A_lm_3var")
addWorksheet(wb, "A_lm_2var")
addWorksheet(wb, "B_lmer_3var")
addWorksheet(wb, "B_lmer_2var")

writeData(wb, "A_lm_3var", results$lm3a)
writeData(wb, "A_lm_2var", results$lm2a)
writeData(wb, "B_lmer_3var", results$lmer3b)
writeData(wb, "B_lmer_2var", results$lmer2b)

saveWorkbook(
  wb,
  file = "/proj/gibbons/2024_easton_vns/data/final/clean_diffabund.xlsx",
  overwrite = TRUE
  )



# head(results$lm2a$lmer_res$y)
# head(results$lm2a$lmer_res)
# colnames(results$lm2a$lmer_res)
# 
# results$lm2a %>%
#   dplyr::select(-signif) %>%
#   pivot_wider(
#     names_from = "term",
#     values_from = c("estimate", "se", "p", "fdr"),
#     names_glue = "{term}_{.value}"
#     )

#   pivot_wider(
#     names_from = covariate,
#     values_from = c("estimate", "se", "lower", "upper", "score_stat", "pval"),
#     names_glue = "{covariate}_{.value}"


meta <- readRDS("/proj/gibbons/2024_easton_vns/data/data/1_preprocessed/01_preprocess_metadata_cleaned_meta.rds")

wb <- createWorkbook()

addWorksheet(wb, "cleaned_metadata")
writeData(wb, "cleaned_metadata", meta)

saveWorkbook(
  wb,
  file = "/proj/gibbons/2024_easton_vns/data/final/clean_meta.xlsx",
  overwrite = TRUE
)