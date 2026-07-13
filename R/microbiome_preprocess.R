return.presence <- function(data, taxa_cols) {
  # Extract abundance data and determine presence (non-zero) for each taxon
  presence_df <- as.data.frame(data) %>%
    dplyr::select(all_of(taxa_cols)) %>%
    mutate(across(all_of(taxa_cols), ~ . != 0))  # Convert to TRUE (non-zero) / FALSE (zero)
  
  # returns named vector of prevalence for each taxon
  return(colSums(presence_df) / nrow(presence_df))
}

presence_by_category <- function(data, category_col) {
  data %>%
    pivot_longer(
      cols = -all_of(category_col),
      names_to = "microbe",
      values_to = "present"
    ) %>%
    group_by(microbe, .data[[category_col]]) %>%
    summarise(presence = mean(present, na.rm=TRUE)) %>%
    pivot_wider(
      names_from = all_of(category_col),
      values_from = presence
    )
}

# filter_and_clr <- function(
#     bug_info,
#     counts,
#     total_presence      = 0.3,
#     batch_presence      = 0.1,
#     impute_method       = "GBM",
#     filt_before_clr     = TRUE
# ) {
#   
#   filt_bugs <- bug_info %>% filter(
#     pres_total     >= total_presence,
#     pres_batch1    >= batch_presence,
#     pres_batch3    >= batch_presence) %>%
#     pull(microbe)
#   
#   if (filt_before_clr) {
#     filt_counts <- counts %>%
#       dplyr::select(any_of(filt_bugs))
#   } else {
#     filt_counts <- counts
#   }
#   
#   imputed_counts <- zCompositions::cmultRepl(filt_counts, label = 0, method = impute_method, z.delete = FALSE)
#   gm <- exp(rowMeans(log(imputed_counts)))
#   clr_mat <- log(imputed_counts) - log(gm)
#   clr_df_raw <- data.frame(clr_mat)
#   
#   res <- list(
#     "clr_df_raw" = clr_df_raw,
#     "filt_bugs" = filt_bugs
#   )
#   
#   return(res)
# }

filter_and_clr <- function(
    bug_info,
    counts,
    total_presence = 0.3,
    batch_presence = 0.1,
    batch_mode = c("all_batches", "min_batches"),
    min_batches = 2,
    batch_cols = NULL,
    impute_method = "GBM",
    filt_before_clr = TRUE,
    scale_clr = FALSE
) {
  batch_mode <- match.arg(batch_mode)
  
  if (is.null(batch_cols)) {
    batch_cols <- grep("^pres_batch", names(bug_info), value = TRUE)
  }
  
  if (batch_mode == "all_batches") {
    filt_bugs <- bug_info %>%
      dplyr::filter(
        pres_total >= total_presence,
        dplyr::if_all(dplyr::all_of(batch_cols), ~ .x >= batch_presence)
      ) %>%
      dplyr::pull(microbe)
    
  } else if (batch_mode == "min_batches") {
    filt_bugs <- bug_info %>%
      dplyr::filter(
        pres_total >= total_presence,
        rowSums(dplyr::across(dplyr::all_of(batch_cols), ~ .x >= batch_presence), na.rm = TRUE) >= min_batches
      ) %>%
      dplyr::pull(microbe)
  }
  
  filt_counts <- if (filt_before_clr) {
    counts %>% dplyr::select(dplyr::any_of(filt_bugs))
  } else {
    counts
  }
  
  imputed_counts <- zCompositions::cmultRepl(
    filt_counts, label = 0, method = impute_method, z.delete = FALSE
  )
  
  gm <- exp(rowMeans(log(imputed_counts)))
  clr_mat <- log(imputed_counts) - log(gm)
  clr_df <- as.data.frame(clr_mat)
  
  if (scale_clr) {
    clr_df <- as.data.frame(scale(clr_df))
  }
  
  list(
    clr_df_raw = clr_df,
    filt_bugs = filt_bugs
  )
}

join_bugs_and_div <- function(
    raw_clr_df,
    raw_alpha_df
) {
  
  bugs_only <- colnames(raw_clr_df)
  
  adiv_features <- raw_alpha_df %>%
    transmute(shannon, richness, berger_parker) %>%
    rownames_to_column("seq_name")
  
  bug_features <- raw_clr_df %>%
    rownames_to_column("seq_name")
  
  raw_merge <- full_join(adiv_features, bug_features, by = "seq_name")
  
  scaled_merge <- raw_merge %>%
    mutate(
      shannon = scale(shannon),
      richness = scale(richness),
      berger_parker = scale(log(berger_parker)),
      across(all_of(bugs_only), ~ scale(.))
    )
  
  features <- list(
    "raw_feat" = raw_merge,
    "norm_feat" = scaled_merge
  )
  
  return(features)
}