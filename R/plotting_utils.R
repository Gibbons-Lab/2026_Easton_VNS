## -----------------------------------------------------------------------------
##
## Script Name: R/plotting_utils.R
##
## Purpose of Script: 
##
## Author: Alyssa N. Easton
##
## Date Created: 2026-05-13
##
## Copyright (c) Alyssa N. Easton, 2026
##
## -----------------------------------------------------------------------------

#' Standardized function to save a list of plots to a multi-page PDF
#' using notebook-based naming prefixes
#' @param plot_list A list of ggplot or grid-based plot objects
#' @param base_name Descriptive name for the file (e.g., "glo_micro_plots")
#' @param dir Directory to save in (e.g., "./figs")
#' @param width Width of the PDF pages in inches (default 5.5)
#' @param height Height of the PDF pages in inches (default 6.0)
save_plots_pdf <- function(plot_list, base_name, dir, width = 5.5, height = 6.0) {
  
  # 1. Enforce that input is a list
  if (!is.list(plot_list)) {
    stop("plot_list must be a list of plot objects.")
  }
  
  # 2. Get dynamic prefix from the current notebook
  prefix <- get_prefix()
  
  # 3. Construct full file path (automatically appends .pdf if missing)
  ext <- if (grepl("\\.pdf$", base_name, ignore.case = TRUE)) "" else ".pdf"
  full_name <- paste0(prefix, base_name, ext)
  full_path <- file.path(dir, full_name)
  
  # 4. Ensure directory exists before saving
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  message(paste("Saving multi-page PDF to:", full_path))
  
  # 5. Open PDF device and guarantee dev.off() executes via on.exit
  pdf(full_path, width = width, height = height)
  on.exit({
    if (names(dev.cur()) == "pdf") {
      dev.off()
    }
  })
  
  # 6. Walk through and print each plot to the PDF
  purrr::walk(plot_list, print)
}

save_heatmap_pdf <- function(ht_object, base_name, dir, width = 10, height = 8, ...) {
  # Get dynamic prefix from the current notebook
  prefix <- get_prefix()
  
  # Construct full file path with .pdf extension
  full_name <- paste0(prefix, base_name, ".pdf")
  full_path <- file.path(dir, full_name)
  
  # Ensure directory exists before saving
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  message(paste("Saving heatmap PDF to:", full_path))
  
  # Open the PDF device with specified dimensions
  pdf(full_path, width = width, height = height)
  
  # Draw the ComplexHeatmap object
  # We pass '...' so you can pass extra arguments like main_heatmap directly
  ComplexHeatmap::draw(ht_object, ...)
  
  # Close the device to actually write the file
  dev.off()
}

save_pages <- function(pages, file, width = 7, height = 7) {
  
  chunk_plots <- function(x, n_per_page) {
    split(x, ceiling(seq_along(x) / n_per_page))
  }
  
  pdf(file, width = width, height = height, useDingbats = FALSE)
  
  purrr::walk(pages, function(pg) {
    
    # CASE 1: one plot per page
    if (pg$layout == "single") {
      purrr::walk(pg$plots, print)
      
      # CASE 2: all plots on one page
    } else if (pg$layout == "grid") {
      p <- patchwork::wrap_plots(pg$plots, ncol = rlang::`%||%`(pg$ncol, 2)) +
        # patchwork::plot_layout(padding = unit(6, "pt")) +
        patchwork::plot_annotation(title = pg$title)
      
      p <- p & theme(plot.margin = margin(10, 15, 10, 15))
      
      print(p)
      
      # ✅ CASE 3: paginated grid (THIS is what you want)
    } else if (pg$layout == "grid_paginated") {
      
      chunks <- chunk_plots(pg$plots, rlang::`%||%`(pg$n_per_page, 4))
      
      purrr::imap(chunks, function(chunk, i) {
        
        p <- patchwork::wrap_plots(chunk, ncol = rlang::`%||%`(pg$ncol, 2)) +
          # patchwork::plot_layout(padding = unit(6, "pt")) +
          patchwork::plot_annotation(
            title = paste0(pg$title, " (page ", i, ")")
          )
        p <- p & theme(plot.margin = margin(10, 15, 10, 15))
        
        print(p)
      })
    }
  })
  
  dev.off()
}

box_easier <- function(
    data,
    dv,
    y_lab,
    italicize_y_lab = TRUE,
    ylab_size   = 25,
    groups      = "inj_trt_group",
    palette     = injtrt_pal,
    x_positions = c(4,10,20,26),
    box_width   = 4,
    axis_length = 30,
    show_legend = FALSE
) {
  # --- Rescale positions and widths to proportions ---
  pos_map <- x_positions / axis_length
  names(pos_map) <- levels(factor(data[[groups]]))  # match group order
  
  # tick_map <- c(
  #   pos_map[1] + ((pos_map[2] - pos_map[1]) / 2),
  #   pos_map[3] + ((pos_map[4] - pos_map[3]) / 2)
  # )
  box_width <- box_width / axis_length              # rescale width
  error_width <- box_width/2
  
  # assign numeric xpos to each row
  data$xpos <- pos_map[as.character(data[[groups]])]
  
  # calculate base ymax, increase height by ~1/4
  dv_vals <- data[[dv]]
  y_range <- range(dv_vals, na.rm=TRUE)
  base_height <- y_range[2] - y_range[1]
  new_height <- base_height * 1.3
  
  new_y_range <- c(
    (y_range[1] - (0.1*base_height)),
    (y_range[2] + (0.2*base_height))
  )
  
  x_lab_nice <- ("  Sham        VNS   ")
  
  if(italicize_y_lab) {
    y_lab_nice <- paste0("*", clean_bug_names(y_lab), "*")
  } else {
    y_lab_nice <- y_lab
  }
  
  p <- ggplot(data, aes(x = xpos, y = .data[[dv]], group = .data[[groups]], fill = .data[[groups]])) +
    geom_boxplot(
      width    = box_width,
      color    = "grey",
      fill     = "transparent",
      outliers = FALSE,
      coef     = 0
    ) +
    stat_summary(
      geom = "errorbar",
      fun.data = "mean_cl_boot",
      width = error_width
    ) +
    geom_point(
      position = position_jitter(width = box_width/4, height = 0), # jitter relative to box width
      shape    = 22,
      size     = 5,
      aes(fill = .data[[groups]]),
      color    = "black",
      stroke   = 1.2
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = pos_map, labels = names(pos_map)) +
    scale_y_continuous(limits = new_y_range) +
    scale_fill_manual(values = palette) +
    coord_fixed(
      ratio = 1,
      clip = "off"
    ) +
    labs(
      x = x_lab_nice,
      y = y_lab_nice,
      fill = NULL
    ) +
    theme_classic() +
    theme(
      legend.position      = if(!show_legend) "none" else "right",
      legend.justification = "left",
      legend.direction     = "vertical",
      legend.background    = element_blank(),
      legend.box.margin    = margin(0, 0, 10, 0),
      plot.margin          = margin(80, 10, 10, 10),
      plot.background      = element_blank(),
      panel.background     = element_blank(),
      
      aspect.ratio         = 1, # square plotting area
      axis.title.y         = element_markdown(size = ylab_size),
      axis.title.x         = element_text(size = 30),
      axis.text.x          = element_blank(),                # element_text(size = 20),
      axis.text.y          = element_text(size = 30),
      
      # just make all ticks and remove later
      axis.ticks.x         = element_blank(),
      axis.ticks.length.y  = unit(0.2, "cm"),
      # axis.ticks.length         = unit(0.2, "cm"),
      
      axis.line                 = element_line(colour = "black", linewidth = 1)
    )
  
  p
}

plot_adjusted_scatter <- function(data, 
                                  x_var, 
                                  y_var, 
                                  covariate_formula, 
                                  residualize_x = FALSE, 
                                  residualize_y = FALSE, 
                                  color_var = "inj_trt_group", 
                                  color_pal = injtrt_colors, 
                                  font_size = 21,
                                  method = "lmer") { # New argument: "lmer" or "lm"
  
  # 1. Setup Data
  # We filter NAs on the outcome, but na.exclude will handle covariate NAs later
  data <- data %>% filter(!is.na(.data[[y_var]]))
  plot_data <- data
  
  # Helper function to fit the correct model type
  fit_model <- function(f, d) {
    if (method == "lmer") {
      # REML=FALSE is standard for lmer comparisons; na.exclude keeps row alignment
      lme4::lmer(f, data = d, REML = FALSE, na.action = na.exclude)
    } else if (method == "lm") {
      # Standard linear model
      stats::lm(f, data = d, na.action = na.exclude)
    } else {
      stop("Method must be 'lmer' or 'lm'")
    }
  }
  
  # 2. Handle Y (outcome)
  if (residualize_y) {
    form_y <- as.formula(paste(y_var, covariate_formula))
    fit_y  <- fit_model(form_y, data)
    
    # na.exclude ensures this matches plot_data length perfectly
    plot_data$y_plot_val <- resid(fit_y) 
    
    # Logic to break label if too long (Fixed to use y_var dynamically)
    y_name_clean <- clean_out_names(y_var)
    y_prefix <- ifelse(nchar(sub("\n.*", "", y_name_clean)) < 15, "Residual ", "Residual\n")
    
  } else {
    plot_data$y_plot_val <- plot_data[[y_var]]
    y_prefix <- ""
  }
  
  # 3. Handle X (Predictor)
  if (residualize_x) {
    form_x <- as.formula(paste(x_var, covariate_formula))
    fit_x  <- fit_model(form_x, data)
    
    plot_data$x_plot_val <- resid(fit_x)
    x_prefix <- "Residual \n"
    
  } else {
    plot_data$x_plot_val <- plot_data[[x_var]]
    x_prefix <- ""
  }
  
  # 4. Generate Labels
  # Check if x_var is in your special list for non-italic names
  is_metric <- x_var %in% c("shannon", "richness", "berger_parker")
  
  clean_x <- gsub("_", " ", clean_bug_names(x_var))
  
  if (is_metric) {
    x_lab_nice <- paste0(x_prefix, " ", clean_x)
  } else {
    x_lab_nice <- paste0(x_prefix, " *", clean_x, "*")
  }
  
  y_lab_nice <- paste0(y_prefix, clean_out_names(y_var))
  
  # 5. Plotting
  ggplot(plot_data, aes(x = x_plot_val, y = y_plot_val)) +
    geom_point(aes(color = .data[[color_var]]), size = 4, alpha = 0.7) +
    
    scale_color_manual(values = color_pal) +
    
    # We use "lm" for the trend line regardless of the residualization method
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    
    labs(x = x_lab_nice, y = y_lab_nice) +
    
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = font_size),
      axis.title.y = element_text(size = font_size),
      axis.title.x = element_markdown(size = font_size)
    )
}

make_manual_legend <- function(palette,
                               labels = c("Sham", "Blast"),
                               box_size = unit(6, "mm"),
                               text_size = 20,
                               hgap = unit(3, "mm"),
                               vgap = unit(2, "mm"),
                               right_pad = unit(3, "mm")) {
  
  n <- length(labels)
  
  # Create text grobs
  text_grobs <- lapply(labels, function(lab) {
    textGrob(lab,
             x = 0,
             y = 0.5,
             just = c("left", "center"),
             gp = gpar(fontsize = text_size))
  })
  
  # Use longest label to define column width (robust)
  longest_label <- labels[which.max(nchar(labels))]
  ref_text <- textGrob(longest_label,
                       gp = gpar(fontsize = text_size))
  
  text_col_width <- grobWidth(ref_text) + right_pad
  
  # Layout: alternating box row + spacer row
  heights <- rep(unit.c(box_size, vgap), length.out = 2*n - 1)
  
  gt <- gtable(
    widths  = unit.c(box_size, hgap, text_col_width),
    heights = heights
  )
  
  row_index <- seq(1, 2*n - 1, by = 2)
  
  for (i in seq_along(labels)) {
    
    rect <- rectGrob(
      width  = box_size,
      height = box_size,
      gp = gpar(fill = palette[labels[i]],
                col = "black",
                lwd = 4.0)
    )
    
    gt <- gtable_add_grob(gt, rect, t = row_index[i], l = 1)
    gt <- gtable_add_grob(gt, text_grobs[[i]], t = row_index[i], l = 3)
  }
  
  gt
}

pdf_split <- function(
    hit_df, data_df,
    sub_height, sub_width,
    n_cols,
    plot_function,
    filename
) {
  
  n_plts <- ifelse(
    plot_function == "box_easier", yes = length(hit_df),
    no = nrow(hit_df)
  )
  
  n_rows <- ceiling(n_plts / n_cols)
  
  pdf_width <- sub_width * n_cols
  pdf_height <- sub_height * n_rows
  
  plot_vector <- vector("list", n_plts)
  
  if (plot_function == "residout_bug") {
    for (i in seq_len(nrow(hit_df))) {
      plot_vector[[i]] <- plot_adjusted_scatter(
        data = data_df,
        x_var = hit_df[[i, "microbe"]],   # e.g., "bug_name"
        y_var = hit_df[[i, "outcome"]],   # e.g., "outcome_name"
        covariate_formula = "~ blast + vns + blast_and_vns + (1|Round_main)",
        residualize_x = FALSE,
        residualize_y = TRUE
      )
      # residout_bug(
      #   bug_name     = hit_df[[i, "microbe"]],
      #   data         = data_df,
      #   outcome_name = hit_df[[i, "outcome"]]
      #   )
    }
    
  } else if (plot_function == "residout_residbug") {
    for (i in seq_len(nrow(hit_df))) {
      
      plot_vector[[i]] <- plot_adjusted_scatter(
        data = data_df,
        x_var = hit_df[[i, "microbe"]],   # 
        y_var = hit_df[[i, "outcome"]],   # e.g., "outcome_name"
        covariate_formula = "~ blast + vns + blast_and_vns + (1|Round_main)",
        residualize_x = TRUE,
        residualize_y = TRUE
      )
      
      # plot_vector[[i]] <- residout_residbug(
      #   bug_name     = hit_df[[i, "mediator"]],
      #   data         = data_df,
      #   outcome_name = hit_df[[i, "outcome"]]
      # )
    }
    
  } else if (plot_function == "box_easier") {
    for (i in 1:length(hit_df)) {  # hit_df can be a vector here
      
      plot_vector[[i]] <- box_easier(
        data_df,
        hit_df[i],
        hit_df[i]
      )
    }
    
  }
  
  pdf(filename, width = pdf_width, height = pdf_height)
  grid.arrange(grobs = plot_vector, ncol = n_cols)
  dev.off()
}

a_pdf_split <- function(
    hit_df, data_df,
    sub_height, sub_width,
    n_cols,
    plot_function,
    filename
) {
  
  n_plts <- ifelse(
    plot_function == "box_easier", yes = length(hit_df),
    no = nrow(hit_df)
  )
  
  n_rows <- ceiling(n_plts / n_cols)
  
  pdf_width <- sub_width * n_cols
  pdf_height <- sub_height * n_rows
  
  plot_vector <- vector("list", n_plts)
  
  if (plot_function == "residout_bug") {
    for (i in seq_len(nrow(hit_df))) {
      plot_vector[[i]] <- plot_adjusted_scatter(
        data = data_df,
        x_var = hit_df[[i, "microbe"]],   # e.g., "bug_name"
        y_var = hit_df[[i, "outcome"]],   # e.g., "outcome_name"
        covariate_formula = "~ blast + vns + blast_and_vns + Round_main",
        residualize_x = FALSE,
        residualize_y = TRUE,
        method = "lm"
      )
    }
    
  } else if (plot_function == "residout_residbug") {
    for (i in seq_len(nrow(hit_df))) {
      
      plot_vector[[i]] <- plot_adjusted_scatter(
        data = data_df,
        x_var = hit_df[[i, "microbe"]],   # 
        y_var = hit_df[[i, "outcome"]],   # e.g., "outcome_name"
        covariate_formula = "~ blast + vns + blast_and_vns + Round_main",
        residualize_x = TRUE,
        residualize_y = TRUE,
        method = "lm"
      )
    }
    
  } else if (plot_function == "box_easier") {
    for (i in 1:length(hit_df)) {  # hit_df can be a vector here
      
      plot_vector[[i]] <- box_easier(
        data_df,
        hit_df[i],
        hit_df[i]
      )
    }
  }
  
  pdf(filename, width = pdf_width, height = pdf_height)
  grid.arrange(grobs = plot_vector, ncol = n_cols)
  dev.off()
}

make_position_map <- function(data, groups, x_positions, axis_length) {
  levs <- levels(factor(data[[groups]]))
  
  pos_map <- x_positions / axis_length
  names(pos_map) <- levs
  
  data$xpos <- pos_map[as.character(data[[groups]])]
  
  list(
    data = data,
    pos_map = pos_map,
    box_width = diff(range(pos_map)) * (x_positions[2] - x_positions[1]) / axis_length,
    error_width = (x_positions[2] - x_positions[1]) / axis_length / 2
  )
}


make_pairs_df <- function(
    fit,
    data,
    response,
    factorA = "injury",
    factorB = "treatment",
    comparison_mode = c("all", "reduced"),
    factorA_level_for_factorB = NULL
) {
  
  comparison_mode <- match.arg(comparison_mode)
  
  levA <- levels(factor(data[[factorA]]))
  
  # default = first level of factorA
  if (is.null(factorA_level_for_factorB)) {
    factorA_level_for_factorB <- levA[1]
  }
  
  form_A <- as.formula(
    paste("pairwise ~", factorA, "|", factorB)
  )
  
  form_B <- as.formula(
    paste("pairwise ~", factorB, "|", factorA)
  )
  
  # ---------------------------------------------------------
  # Simple effects
  # ---------------------------------------------------------
  
  res_A <- emmeans(fit, form_A)
  res_B <- emmeans(fit, form_B)
  
  # ---------------------------------------------------------
  # factorA simple effects
  # ---------------------------------------------------------
  
  pairs_A <- tidy(res_A$contrasts) %>%
    tidyr::separate(
      contrast,
      into = c("g1", "g2"),
      sep = " - "
    ) %>%
    dplyr::mutate(
      group1 = paste(g1, !!rlang::sym(factorB)),
      group2 = paste(g2, !!rlang::sym(factorB)),
      effect_type = "factorA"
    )
  
  # ---------------------------------------------------------
  # factorB simple effects
  # ---------------------------------------------------------
  
  pairs_B <- tidy(res_B$contrasts) %>%
    tidyr::separate(
      contrast,
      into = c("g1", "g2"),
      sep = " - "
    ) %>%
    dplyr::mutate(
      group1 = paste(!!rlang::sym(factorA), g1),
      group2 = paste(!!rlang::sym(factorA), g2),
      effect_type = "factorB"
    )
  
  # keep only one factorA level for factorB comparisons
  if (comparison_mode == "reduced") {
    
    pairs_B <- pairs_B %>%
      dplyr::filter(
        .data[[factorA]] == factorA_level_for_factorB
      )
  }
  
  # ---------------------------------------------------------
  # Combine
  # ---------------------------------------------------------
  
  pairs <- dplyr::bind_rows(pairs_A, pairs_B)
  
  # ---------------------------------------------------------
  # Multiple testing correction
  # ---------------------------------------------------------
  
  pairs <- pairs %>%
    dplyr::mutate(
      p_adj = p.adjust(p.value, method = "bonferroni"),
      stars = dplyr::case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ "ns"
      ),
      label = ifelse(p_adj < 0.05, stars, "ns")
    )
  
  # ---------------------------------------------------------
  # y positions
  # ---------------------------------------------------------
  
  y_max <- max(data[[response]], na.rm = TRUE)
  y_min <- min(data[[response]], na.rm = TRUE)
  y_range <- y_max - y_min
  
  n_pairs <- nrow(pairs)
  
  pairs$y.position <- y_max + (
    pmax(1, seq_len(n_pairs) - 1) * 0.2 * y_range # compares against 1, returns 1 if zero
  )
  
  pairs %>%
    dplyr::select(
      group1,
      group2,
      y.position,
      label
    )
}

make_anova_pairs <- function(
    data,
    dv,
    factorA,
    factorB,
    pos_map,
    fit_fun = aov,
    comparison_mode = "all",
    factorA_level_for_factorB = NULL
) {
  
  f <- reformulate(
    termlabels = paste0(factorA, "*", factorB),
    response = dv
  )
  
  fit <- fit_fun(f, data = data)
  
  pairs_df <- make_pairs_df(
    fit = fit,
    data = data,
    response = dv,
    factorA = factorA,
    factorB = factorB,
    comparison_mode = comparison_mode,
    factorA_level_for_factorB = factorA_level_for_factorB
  ) %>%
    dplyr::mutate(
      group1 = gsub(" ", " + ", group1),
      group2 = gsub(" ", " + ", group2),
      xmin = pos_map[group1],
      xmax = pos_map[group2]
    )
  
  pairs_df
}

make_y_range <- function(x, expand = 0.2) {
  r <- range(x, na.rm = TRUE)
  h <- diff(r)
  c(r[1] - expand * h, r[2] + expand * h)
}

plot_box_with_stats <- function(
    data,
    dv,
    groups,
    pos_map,
    pairs_df,
    palette,
    y_lab,
    x_lab = NULL,
    italicize_y_lab = FALSE,
    ylab_size = 25,
    box_width,
    error_width,
    show_legend = FALSE,
    y_expand = 0.2
) {
  
  y_range <- make_y_range(data[[dv]], expand = y_expand)
  
  if (italicize_y_lab) {
    y_lab <- paste0("*", clean_bug_names(y_lab), "*")
  }
  
  cluster_centers <- c(
    Sham = mean(pos_map[c(levels(data[[groups]])[[1]], levels(data[[groups]])[[2]])]),   # adjust to your levels
    VNS  = mean(pos_map[c(levels(data[[groups]])[[3]], levels(data[[groups]])[[4]])])
  )
  
  ggplot(
    data,
    aes(
      x = xpos,
      y = .data[[dv]],
      group = .data[[groups]],
      fill = .data[[groups]]
    )
  ) +
    geom_boxplot(
      width = box_width,
      color = "grey",
      fill = "transparent",
      outliers = FALSE,
      coef = 0
    ) +
    stat_summary(
      geom = "errorbar",
      fun.data = "mean_cl_boot",
      width = error_width
    ) +
    geom_point(
      position = position_jitter(width = box_width / 4, height = 0),
      shape = 22,
      size = 5,
      color = "black",
      stroke = 1.2
    ) +
    stat_pvalue_manual(
      pairs_df,
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      label = "label",
      bracket_size = 1.0,
      size = 5,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = cluster_centers,
      labels = names(cluster_centers)
    ) +
    scale_fill_manual(values = palette) +
    coord_fixed(ratio = 1, clip = "off", ylim = y_range) +
    labs(x = x_lab, y = y_lab, fill = NULL) +
    theme_classic() +
    theme(
      legend.position = if (!show_legend) "none" else "right",
      plot.margin = margin(80, 10, 10, 10),
      aspect.ratio = 1,
      axis.title.y = element_markdown(size = ylab_size),
      axis.text.x = element_text(size = 25),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 25),
      axis.line = element_line(colour = "black", linewidth = 1),
      axis.ticks.x = element_line(color = "black", linewidth = 0.7)
    )
}

# labels
get_label <- function(dv, label_df) {
  out <- as.character(label_df$label[match(dv, label_df$dv)])
  ifelse(is.na(out), dv, out)
}

plot_multiple_dv <- function(
    data,
    dvs,
    groups,
    factorA,
    factorB,
    x_positions,
    axis_length,
    palette,
    label_df = NULL,
    italicize = FALSE,
    ylab_size = 30,
    show_legend = FALSE,
    fit_fun = aov,
    comparison_mode = "all",
    factorA_level_for_factorB = NULL,
    return_plots = TRUE
) {
  # --- positions (shared across all DVs) ---
  pos <- make_position_map(
    data = data,
    groups = groups,
    x_positions = x_positions,
    axis_length = axis_length
  )
  
  data <- pos$data
  
  # --- loop over DVs ---
  plot_list <- lapply(dvs, function(dv) {
    
    # label handling
    y_lab <- if (!is.null(label_df)) {
      get_label(dv, label_df)
    } else {
      dv
    }
    
    # stats
    pairs_df <- make_anova_pairs(
      data = data,
      dv = dv,
      factorA = factorA,
      factorB = factorB,
      pos_map = pos$pos_map,
      fit_fun = fit_fun,
      comparison_mode = comparison_mode,
      factorA_level_for_factorB = factorA_level_for_factorB
    )
    
    # plot
    plot_box_with_stats(
      data = data,
      dv = dv,
      groups = groups,
      pos_map = pos$pos_map,
      pairs_df = pairs_df,
      palette = palette,
      y_lab = y_lab,
      italicize_y_lab = italicize,
      ylab_size = ylab_size,
      box_width = pos$box_width,
      error_width = pos$error_width,
      show_legend = show_legend
    )
  })
  
  names(plot_list) <- dvs
  
  if (return_plots) {
    return(plot_list)
  } else {
    invisible(plot_list)
  }
}

generate_adjusted_plots <- function(
    data, 
    x_vars, 
    y_vars, 
    covariate_formula = "~ blast + vns + blast_and_vns + Round_main", 
    residualize_x = TRUE, 
    residualize_y = TRUE, 
    method = "lm",
    color_var = "inj_trt_group", 
    color_pal = injtrt_colors, # Assumes this exists in your environment
    font_size = 21
) {
  
  # 1. Input Validation
  if (length(x_vars) != length(y_vars)) {
    stop("Error: x_vars and y_vars must be the same length to form pairs.")
  }
  
  # 2. Iterate over the pairs to generate the list of plots
  plot_list <- lapply(seq_along(x_vars), function(i) {
    x_var <- x_vars[i]
    y_var <- y_vars[i]
    
    # Isolate data and filter NAs for this specific outcome pair
    pair_data <- data %>% filter(!is.na(.data[[y_var]]))
    
    # Helper to fit the correct model type
    fit_model <- function(f, d) {
      if (method == "lmer") {
        lme4::lmer(f, data = d, REML = FALSE, na.action = na.exclude)
      } else if (method == "lm") {
        stats::lm(f, data = d, na.action = na.exclude)
      } else {
        stop("Method must be 'lmer' or 'lm'")
      }
    }
    
    # --- Y (Outcome) Handling ---
    if (residualize_y) {
      form_y <- as.formula(paste(y_var, covariate_formula))
      fit_y  <- fit_model(form_y, pair_data)
      pair_data$y_plot_val <- resid(fit_y)
      
      y_name_clean <- clean_out_names(y_var) # Assumes helper function exists
      y_prefix <- ifelse(nchar(sub("\n.*", "", y_name_clean)) < 15, "Residual ", "Residual\n")
    } else {
      pair_data$y_plot_val <- pair_data[[y_var]]
      y_name_clean <- clean_out_names(y_var)
      y_prefix <- ""
    }
    
    # --- X (Predictor) Handling ---
    if (residualize_x) {
      form_x <- as.formula(paste(x_var, covariate_formula))
      fit_x  <- fit_model(form_x, pair_data)
      pair_data$x_plot_val <- resid(fit_x)
      x_prefix <- "Residual \n"
    } else {
      pair_data$x_plot_val <- pair_data[[x_var]]
      x_prefix <- ""
    }
    
    # --- Generate Labels ---
    is_metric <- x_var %in% c("shannon", "richness", "berger_parker")
    clean_x <- gsub("_", " ", clean_bug_names(x_var)) # Assumes helper function exists
    
    if (is_metric) {
      x_lab_nice <- paste0(x_prefix, " ", clean_x)
    } else {
      x_lab_nice <- paste0(x_prefix, " *", clean_x, "*")
    }
    y_lab_nice <- paste0(y_prefix, y_name_clean)
    
    # --- Plotting ---
    p <- ggplot(pair_data, aes(x = x_plot_val, y = y_plot_val)) +
      geom_point(aes(color = .data[[color_var]]), size = 5, alpha = 0.7) +
      scale_color_manual(values = color_pal) +
      
      geom_smooth(method = "lm",se = FALSE, color = "black", linetype = "dashed") +
      # Group-specific lines of fit matching the point colors
      geom_smooth(
        aes(color = .data[[color_var]]), 
        method = "lm", 
        se = FALSE, 
        alpha = 0.3 
      ) +
      
      labs(x = x_lab_nice, y = y_lab_nice) +
      theme_bw() +
      theme(
        aspect.ratio = 1, # Forces the plot panel to be a square
        legend.position = "none",
        axis.text = element_text(size = font_size),
        axis.title.y = element_text(size = font_size),
        axis.title.x = ggtext::element_markdown(size = font_size)
      )
    
    return(p)
  })
  
  # Name the list elements for easy access later 
  # e.g., my_plots[["shannon_vs_outcome1"]]
  names(plot_list) <- paste0(x_vars, "_vs_", y_vars)
  
  return(plot_list)
}

# 1) Main biplot
biplot_main <- function(pc_df, prop_var,
                        x_axis, y_axis,
                        color_var, color_palette, 
                        shape_var, shape_palette,
                        point_size = 3,
                        legend_text_size = 7) {
  ggplot() +
    geom_point(
      data = pc_df,
      aes(x = .data[[x_axis]],
          y = .data[[y_axis]],
          color = .data[[color_var]],
          shape = .data[[shape_var]]),
      size = point_size, alpha = 0.8
    ) +
    scale_shape_manual(values = shape_palette) +
    stat_ellipse(
      data = pc_df, aes(x = .data[[x_axis]], y = .data[[y_axis]], color = .data[[color_var]]),# Group),
      level = 0.68, type = "norm", linetype = "dashed"
    ) +
    scale_color_manual(values = color_palette) +
    # arrows, if load_df is supplied (see old confo "merging dataframes by keys" if you want this code)
    theme_bw() +        # theme_bw has outlines around plot area by default
    labs(
      x = paste0(x_axis, " (", round(prop_var[[x_axis]] * 100, 1), "%)"),
      y = paste0(y_axis, " (", round(prop_var[[y_axis]] * 100, 1), "%)")
    ) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size),
      axis.title = element_text(size = 20),
      # axis.text  = element_blank(),
      # axis.ticks = element_blank(),         # remove all axis ticks
      # panel.grid.minor = element_blank(),   # remove gridlines
      # panel.grid.major = element_blank()
    )
}

extract.legend <- function(ggplot_obj) {
  g <- ggplotGrob(ggplot_obj)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

# 2) Marginal boxplot for PC1 (horizontal)
marginal_box_x <- function(pc_df, x_axis, color_var, color_palette) {
  # reverse order of group variable levels
  pc_df <- pc_df %>%
    mutate(!!color_var := fct_rev(.data[[color_var]]))
  
  ggplot(pc_df, aes(x = .data[[x_axis]], y = .data[[color_var]], fill = .data[[color_var]])) +
    geom_boxplot(alpha = 0.8, ) +
    scale_fill_manual(values = color_palette) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # panel.grid = element_blank()
    )
}

# 3) Marginal boxplot for PC2 (vertical)
marginal_box_y <- function(pc_df, y_axis, color_var, color_palette) {
  ggplot(pc_df, aes(
    x = .data[[color_var]],
    y = .data[[y_axis]],
    fill = .data[[color_var]]
  )) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = color_palette) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # panel.grid = element_blank()
    )
}


combined.biplot <- function(pc_df, prop_var,
                            x_axis, y_axis,
                            color_var, color_palette, 
                            shape_var, shape_palette,
                            rel_widths = c(4, 1),
                            point_size = 3,
                            legend_text_size = 7) {
  
  biplot <- biplot_main(pc_df, prop_var, x_axis, y_axis, color_var,
                        color_palette, shape_var, shape_palette, point_size)
  
  # 2. Extract legend
  leg <- extract.legend(biplot)
  
  # 3. Extract x and y axes of main biplot
  build_main <- ggplot_build(biplot)
  x_limits <- build_main$layout$panel_scales_x[[1]]$range$range
  y_limits <- build_main$layout$panel_scales_y[[1]]$range$range
  
  # 4. Remove legends from all panels
  biplot <- biplot + theme(legend.position = "none")
  
  xplt <- 
    marginal_box_x(
      pc_df, x_axis,
      color_var, color_palette
    ) +
    theme(legend.position = "none") +
    scale_x_continuous(limits = x_limits)
  
  yplt <-
    marginal_box_y(
      pc_df, y_axis,
      color_var, color_palette
    ) +
    theme(legend.position = "none") +
    scale_y_continuous(limits = y_limits)
  
  ### align p_x & p_main on the vertical edge
  xplt_align     <- cowplot::align_plots(xplt, biplot, align = "v", axis = "l")[[1]]
  ### align p_main & p_y on the horizontal edge
  yplt_align     <- cowplot::align_plots(biplot, yplt, align = "h", axis = "b")[[2]]
  
  # 5. Arrange the grid:
  #         [ xplt_align ][     leg    ]      
  #         [ biplot     ][ yplt_align ]
  
  top_row    <- plot_grid(xplt_align,        leg, ncol = 2, rel_widths = rel_widths, align = "hv")
  bottom_row <- plot_grid(biplot,     yplt_align, ncol = 2, rel_widths = rel_widths)
  
  final_plot <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = rev(rel_widths))  # defaults to main 4/5, legend 1/5
  
  return(final_plot)
}

extract_two_way_anova_table <- function(data, dvs, factorA, factorB, fit_fun = aov) {
  one_dv <- function(dv) {
    d <- data %>%
      dplyr::filter(!is.na(.data[[dv]])) %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(c(factorA, factorB)), as.factor)
      )
    
    f <- stats::reformulate(paste0(factorA, "*", factorB), response = dv)
    fit <- fit_fun(f, data = d)
    
    an <- broom::tidy(fit) %>%
      dplyr::filter(term != "Residuals")
    
    get_term <- function(term_name) {
      row <- an %>% dplyr::filter(term == term_name)
      if (nrow(row) == 0) {
        return(tibble::tibble(statistic = NA_real_, p.value = NA_real_))
      }
      row %>% dplyr::select(statistic, p.value)
    }
    
    int_term <- paste0(factorA, ":", factorB)
    a_term <- factorA
    b_term <- factorB
    
    tibble::tibble(
      dv = dv,
      interaction_F = get_term(int_term)$statistic[[1]],
      interaction_p = get_term(int_term)$p.value[[1]],
      injury_F = get_term(a_term)$statistic[[1]],
      injury_p = get_term(a_term)$p.value[[1]],
      treatment_F = get_term(b_term)$statistic[[1]],
      treatment_p = get_term(b_term)$p.value[[1]]
    )
  }
  
  purrr::map_dfr(dvs, one_dv)
}