
# Assumes count_table contains non-negative count-like estimated read counts (rep(seq_along(x), x) expects counts)
rarefy <- function(count_table) {
  # just use minimum n_reads
  n <- min(rowSums(count_table, na.rm = TRUE))
  counts_rarefied <- t(apply(count_table, 1, function(x) {
    if (sum(x) < n) {
      # not enough reads, return original counts
      return(x)
    } else {
      sampled <- sample(rep(seq_along(x), x), n)
      tab <- tabulate(sampled, nbins = length(x))
      return(tab)
    }
  }))
  colnames(counts_rarefied) <- colnames(count_table)
  rownames(counts_rarefied) <- rownames(count_table)
  
  zero_cols <- colnames(counts_rarefied)[colSums(counts_rarefied) < 1]
  counts_rarefied[, !(colnames(counts_rarefied) %in% zero_cols)]
}

relabund.betadiv.analyze <- function(abund, meta_filt, id_col = "seq_name", perm_formula, condition, n_perm) {
  # generate Bray-Curtis Dissimilarity Matrix
  filt_abund <- abund %>%
    filter(.data[[id_col]] %in% meta_filt[[id_col]]) %>%
    column_to_rownames("seq_name") %>%
    as.matrix()
  
  # do some rarefaction (since this isn't a count table, we can't take a random sample, but oh well)
  min_nonzero <- apply(filt_abund, 1, function(x) min(x[x > 0], na.rm = TRUE))
  cutoff <- max(min_nonzero)
  filt_abund[filt_abund < cutoff] <- 0
  filt_abund <- filt_abund[, colSums(filt_abund != 0) > 0, drop = FALSE]
  
  bc <- vegdist(filt_abund, method = "bray")
  
  if (!identical(rownames(meta_filt), labels(bc))) {
    meta_filt <- meta_filt[match(labels(bc), rownames(meta_filt)), , drop = FALSE]
  }
  
  # Do traditional PCoA, outputting variance plot
  pcoa <- vegan::wcmdscale(bc, eig = TRUE)
  pcoa_var <- var.explain(pcoa = pcoa)
  plt_var <- var.part.plot(
    pcoa_var, "Variance in Bray-Curtis Beta Diversity Explained by PCo",
    num_pc = 8) + theme_classic()
  
  # Do conditioned PCoA (corrected for some nuisance variable)
  cpcoa_formula <- as.formula(paste("bc ~", condition))
  cpcoa <- vegan::capscale(cpcoa_formula, data = meta_filt)
  
  # extract eigenvalues for CCA and CA
  eig <- cpcoa$CA$eig                             # PCs orthogonal to batch (conditioned)
  prop_var <- eig / sum(eig)                      # proportion of non-batch-attributed variance explained by each PC
  
  # do PERMANOVA
  perm_formula <- as.formula(paste("bc ~", perm_formula))
  
  perm_marg <- adonis2(
    perm_formula,
    data = meta_filt,
    permutations = n_perm,
    by = "margin"
  ) %>%
    tidy() %>% mutate(across(c(SumOfSqs, R2, statistic), ~ round(.x, 3)))
  
  perm_terms <- adonis2(
    perm_formula,         # bc ~ X + Y + Z
    data = meta_filt,
    permutations = n_perm,
    by = "terms"
  ) %>%
    tidy() %>% mutate(across(c(SumOfSqs, R2, statistic), ~ round(.x, 3)))
  
  list(
    meta_filt   = meta_filt,
    bc_dist_mat = bc,
    pcoa        = pcoa,
    pcoa_var    = pcoa_var,
    pcoa_var_plt= plt_var,
    cpcoa       = cpcoa,
    cpcoa_var   = prop_var,
    perm_marg   = perm_marg,
    perm_terms  = perm_terms
  )
}


# rarefied count matrix
count.betadiv.analyze <- function(
    rarefied_counts,                 # rarefied count matrix
    meta_filt,                       # metadata with a row (id_col) that matches count matrix row names/col names
    id_col = "seq_name",             # name of column containing the ID
    perm_formula,                    # RHS formula for PERMANOVA (i.e. "Round_main + blast + vns + blast_and_vns")
    condition,                       # if doing a conditional PCoA (residualized), which variable are you removing?
    n_perm = 999                     # number of permutations
) {
  rownames(meta_filt) <- meta_filt[[id_col]]
  # generate Bray-Curtis Dissimilarity Matrix
  bc <- vegdist(rarefied_counts, method = "bray")
  
  if (!identical(rownames(meta_filt), labels(bc))) {
    meta_filt <- meta_filt[match(labels(bc), rownames(meta_filt)), , drop = FALSE]
  }
  
  # Do traditional PCoA, outputting variance plot
  pcoa <- vegan::wcmdscale(bc, eig = TRUE)
  pcoa_var <- var.explain(pcoa = pcoa)
  plt_var <- var.part.plot(
    pcoa_var, "Variance in Bray-Curtis Beta Diversity Explained by PCo",
    num_pc = 8) + theme_classic()
  
  pcoa_var <- pcoa_var %>%
    mutate(axis = paste0("Dim", axis))
  
  pcoa_var_vector <- pcoa_var %>% pull(var_explained)
  names(pcoa_var_vector) <- pcoa_var %>% pull(axis)
  
  # Do conditioned PCoA (corrected for some nuisance variable)
  cpcoa_formula <- as.formula(paste("bc ~", condition))
  cpcoa <- vegan::capscale(cpcoa_formula, data = meta_filt)
  
  # extract eigenvalues for CCA and CA
  eig <- cpcoa$CA$eig                             # PCs orthogonal to batch (conditioned)
  prop_var <- eig / sum(eig)                      # proportion of non-batch-attributed variance explained by each PC
  
  # do PERMANOVA
  perm_formula <- as.formula(paste("bc ~", perm_formula))
  
  perm_marg <- adonis2(
    perm_formula,
    data = meta_filt,
    permutations = n_perm,
    by = "margin"
  ) %>%
    tidy() %>% mutate(across(c(SumOfSqs, R2, statistic), ~ round(.x, 3)))
  
  perm_terms <- adonis2(
    perm_formula,         # bc ~ X + Y + Z
    data = meta_filt,
    permutations = n_perm,
    by = "terms"
  ) %>%
    tidy() %>% mutate(across(c(SumOfSqs, R2, statistic), ~ round(.x, 3)))
  
  list(
    meta_filt   = meta_filt,
    bc_dist_mat = bc,
    pcoa        = pcoa,
    pcoa_var    = pcoa_var_vector,
    pcoa_var_plt= plt_var,
    cpcoa       = cpcoa,
    cpcoa_var   = prop_var,
    perm_marg   = perm_marg,
    perm_terms  = perm_terms
  )
}

alpha_div <- function(raref_counts) {
  shannon    <- vegan::diversity(raref_counts, index = "shannon", MARGIN = 1) # row-wise
  n_feats    <- rowSums(raref_counts >= 1)
  simpson    <- vegan::diversity(raref_counts, index = "simpson", MARGIN = 1)
  invsimpson <- vegan::diversity(raref_counts, index = "invsimpson")
  richness   <- rowSums(raref_counts > 0)
  pielou     <- shannon / log(richness)
  berg_park  <- apply(raref_counts, 1, function(x) {
    max(x) / sum(x)
  })
  
  data.frame(
    shannon       = shannon,
    simpson       = simpson,
    invsimpson    = invsimpson,
    richness      = richness,
    pielou        = pielou,
    berger_parker = berg_park
  )
}

########################## PCoA Variance Explained  ############################
# From a vegan pcoa, calculate variance explained by each PC, output in df
var.explain <- function(pcoa) {
  
  ev <- pcoa$eig
  total_var <- sum(abs(ev))
  var_explained <- ev / total_var           # don't take absolute value of eigenvalue in numerator
  cumulative_var <- cumsum(var_explained)
  
  pcoa_var <-
    data.frame(
      axis = seq_along(var_explained),
      var_explained,
      cumulative_var
    )
  
  return(pcoa_var)
}

var.part.plot <- function(var_df, title_char, num_pc = 5) {
  
  top_pcs <- var_df$var_explained[1:num_pc] %>%
    data.frame(PC = paste0("PC", 1:num_pc), var = .)
  
  plt <- ggplot(top_pcs, aes(x = PC, y = var)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = title_char,
         x = "Principal Coordinate",
         y = "Variance Explained") +
    theme_bw()
  return(plt)
}