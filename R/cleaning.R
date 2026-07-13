########################  clean vectors of names  ########################
clean_vec <- function(name_vec) {
  tmp_df <- data.frame(matrix(NA, nrow=5, ncol=length(name_vec)))
  colnames(tmp_df) <- name_vec
  tmp_df_clean <- clean_names(tmp_df, case = "none")
  colnames(tmp_df_clean)
}

clean_meta <- function(meta) {
  meta %>%
    mutate(Round_main    = droplevels(Round_main)) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%  # same as meta <- meta[, colSums(is.na(meta)) != nrow(meta)]
    relocate(c("cohort", "inj_trt_group"), .before = 3) %>%
    relocate(c("animal_code", "seq_name"), .before = 1) %>%
    dplyr::select(-c("cohort", "file_name", "feces_tp"))
}

clean_bug_names <- function(name_vector) {
  clean_vector <- name_vector %>%
    gsub(" sp\\. ", "_sp_", .) %>%
    gsub("([0-9]+)_([0-9]+)", "\\1-\\2", .) %>%
    gsub("_", " ", .) %>%
    gsub("unclassified", "uncl.", .) %>%
    gsub("shannon", "Shannon Diversity", .) %>%
    gsub("berger parker", "Berger-Parker Index", .) %>%
    gsub("richness", "Richness (# Species)", .)
  
  # italicized <- clean_vector %>% paste0("*", ., "*")
  
  return(clean_vector)
}

clean_out_names <- function(name_vector) {
  outcome_naming <- c(
    `EtOH_pref_24h_final`      = "Preference\n(ethanol/water)",
    `startle_PPI`              = "Pre-pulse\ninhibition (AUC)",
    `startle_habituation`      = "Startle\nHabituation (block 5/block 1)",
    `EtOH_24h_intakeave_final` = "Daily alcohol\nintake (g/kg)",
    `OFB_Cdistance`            = "Distance\ntraveled in open field (m)",
    `OFB_Clatency`             = "Delay to first\ncenter entry in open field (s)",
    `PD_choice_filt`           = "Percent choice\nof risky option (AUC)",
    `acute_cyto_PC1`           = "Cytokine PC1",
    `acute_cyto_PC2`           = "Cytokine PC2",
    `acute_cyto_PC3`           = "Cytokine PC3",
    `acute_cyto_PC4`           = "Cytokine PC4",
    `acute_cyto_PC5`           = "Cytokine PC5"
  )
  
  out <- unname(outcome_naming[name_vector])
  out[is.na(out)] <- name_vector[is.na(out)]
  out
}

# Create the wrapper function
format_mediator_labels <- function(x) {
  # 1. Run your existing cleaning
  cleaned <- clean_bug_names(x)
  
  # 2. Apply conditional formatting
  ifelse(x %in% c("shannon", "richness", "berger_parker"), 
         cleaned, 
         paste0("*", cleaned, "*"))
}