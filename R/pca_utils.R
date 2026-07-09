
# wrapper function to run PCA, calculate variance explained,
# keep PCs with R2 > threshold, join PCs back to metadata

run_pca <- function(meta, block, var_threshold = 0.05) {
  # 1. Prep metadata
  meta0 <- meta %>%
    mutate(
      Round_main = droplevels(Round_main),
      SampleID  = as.character(seq_name)
    )
  
  # 2. Extract & scale the variables matrix
  mat <- meta0 %>%
    dplyr::select(matches(block)) %>%
    as.matrix() %>%
    scale()
  rownames(mat) <- meta0$SampleID
  
  # 3. Run PCA
  pca_obj <- prcomp(mat, center = TRUE)
  
  # 4. Variance explained table
  vars <- summary(pca_obj)$importance %>%
    as_tibble(rownames = "measure") %>%
    pivot_longer(-measure, names_to = "PC", values_to = "value") %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    dplyr::select(PC, `Proportion of Variance`, `Cumulative Proportion`)
  
  # 5. Choose which PCs to keep
  keep_pcs <- vars %>%
    filter(`Proportion of Variance` >= var_threshold) %>%
    pull(PC)
  
  # 6. Augment metadata with the kept PCs
  aug <- broom::augment(pca_obj, data = as.data.frame(mat)) %>%
    dplyr::select(.rownames, starts_with(".fittedPC")) %>%
    rename_with(~ stringr::str_remove(.x, "^\\.fitted"), starts_with(".fitted")) %>%
    dplyr::select(.rownames, all_of(keep_pcs)) %>%
    rename(SampleID = .rownames)
  
  meta_enriched <- meta0 %>%
    left_join(aug, by = "SampleID")
  
  # 7. Collect the loadings matrix
  loadings <- as_tibble(pca_obj$rotation, rownames = "Feature")
  
  # 8. Return everything in a list
  list(
    meta     = meta_enriched,
    pca      = pca_obj,
    vars     = vars,
    loadings = loadings
  )
}

# barplot of variance explained by PC (defaults to top 10 PCs)
plot_pca_variance <- function(res, num = 10) {
  n_pc_var <- res$vars %>%
    # extract first `num` rows of the vars tibble
    slice_head(n = num) %>%
    # make PC into an ordered factor
    mutate(PC = factor(PC, levels = paste0("PC", seq_len(num))),
           prop_var = `Proportion of Variance`)
  
  ggplot(n_pc_var, aes(x = PC, y = prop_var)) +
    geom_col(fill = "steelblue") +
    ylab("Proportion of Variance Explained") +
    xlab("Principal Component") +
    theme_minimal()
}

# Vector plot of PCA feature loadings over PC1 and PC2
plot_pca_vectors <- function(res, i_vectors = 10) {
  # res$loadings has Feature, PC1, PC2, …
  load_df <- res$loadings %>%
    as.data.frame() %>%
    dplyr::select(Feature, PC1, PC2) %>%
    mutate(length = sqrt(PC1^2 + PC2^2)) %>%
    arrange(desc(length)) %>%
    slice_head(n = i_vectors)    # we keep only the i longest vectors in PC1/PC2 space
  
  # Colors & linetypes
  feats <- load_df$Feature
  # vector_colors <- RColorBrewer::brewer.pal(max(i_vectors, 3), "Set3")   # this is only for i ≤ 12
  # __________________ Set up color/linetype coding ________________________
  set1 <- RColorBrewer::brewer.pal(9, "Set1")
  set2 <- RColorBrewer::brewer.pal(8, "Set2")
  set3 <- RColorBrewer::brewer.pal(12, "Set3")
  dark2 <- RColorBrewer::brewer.pal(8, "Dark2")
  paired <- RColorBrewer::brewer.pal(12, "Paired")
  
  # If we need more colors, add a gradient
  all_palette <- unique(c(set1, set2, set3, dark2, paired))
  if (i_vectors > length(all_palette)) {
    extra_colors <- colorRampPalette(all_palette)(i_vectors - length(all_palette))
    vector_colors <- c(all_palette, extra_colors)
  } else {
    vector_colors <- all_palette[1:i_vectors]
  }
  
  # set the names of color vector as PCA feature names
  names(vector_colors) <- feats
  
  # define line types
  line_types   <- setNames(rep(c("solid", "22"), length.out = i_vectors), feats)
  # _____________________________________________________________________________________
  
  max_comp <- max(load_df$length)
  
  ggplot(load_df, aes(x = 0, y = 0, xend = PC1, yend = PC2)) +
    geom_segment(aes(color = Feature, linetype = Feature),
                 linewidth = 0.8) +
    geom_point(aes(x = PC1, y = PC2, color = Feature), size = 2) +
    coord_fixed() +
    xlim(c(-max_comp, max_comp)) +
    ylim(c(-max_comp, max_comp)) +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    labs(x = "PC1 loading", y = "PC2 loading") +
    scale_color_manual(values = vector_colors) +
    scale_linetype_manual(values = line_types)
}

plot_pca_biplot <- function(res, group_col, group_colors, i_vectors = 10) {
  # points + confidence ellipses (1 st.dev from group centroid)
  p <- ggplot(res$meta, aes_string("PC1", "PC2", color = group_col)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(level = 0.68, type = "norm", linetype = "dashed") +
    scale_color_manual(values = group_colors) +
    theme_minimal() + 
    theme(legend.title = element_blank())
  
  if (i_vectors <= 0) {
    return(p)
  }
  
  # arrows from loadings
  vec <- plot_pca_vectors(res, i_vectors)$data  # if you expose load_df
  p + geom_segment(data = vec,
                   aes(x = 0, y = 0, xend = PC1, yend = PC2, color = "black"),
                   arrow = arrow(length = unit(0.2, "cm")))
}

# Fixed marginal boxplot function
marginal_boxplot <- function(res, PC = "PC1", vertical = TRUE, group_var = "inj_trt_group", group_colors) {
  pc_df <- as.data.frame(res$meta)
  
  # Make sure PC column exists in the data
  if(!(PC %in% colnames(pc_df))) {
    stop(paste("PC column", PC, "not found in data"))
  }
  
  # Create the plot
  plt <- ggplot(data = pc_df, aes(x = .data[[group_var]], y = .data[[PC]], fill = .data[[group_var]])) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    scale_fill_manual(values = group_colors) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  # Apply coord_flip conditionally
  if(vertical) {
    plt <- plt + coord_flip() +
      theme(axis.text.y = element_blank())
  }
  
  return(plt)
}
