# mediation_utils.R

safe_pi <- function(x, floor = 1e-3) {
  if (is.finite(x)) max(floor, x) else floor
}

safe_eval_thr <- function(thr, value) {
  # Try original value first, then rounded versions, mimicking authors' fallback idea
  vals <- c(value, round(value, 10:0))
  for (v in vals) {
    out <- tryCatch(thr(v), error = function(e) NA_real_)
    if (is.finite(out)) return(out)
  }
  NA_real_
}

# compute FDR threshold using sorted-statistic search
find_discrete_dact_threshold <- function(thr, Ts) {
  sorted_index <- order(Ts)
  lower_bound <- 1
  upper_bound <- length(Ts)
  
  if (length(Ts) < 2) return(NA_real_)
  
  while (upper_bound - lower_bound > 1) {
    mid_index <- floor((lower_bound + upper_bound) / 2)
    mid_value <- Ts[sorted_index[mid_index]]
    
    x <- safe_eval_thr(thr, mid_value)
    if (!is.finite(x)) return(NA_real_)
    
    if (x < 0) {
      lower_bound <- mid_index
    } else {
      upper_bound <- mid_index
    }
  }
  
  Ts[sorted_index[lower_bound]]
}

# compute FDR threshold using a continuous root-finding step (uniroot)
find_root_dact_threshold <- function(thr, interval = c(0, 1), tol = 1e-10) {
  f0 <- tryCatch(thr(interval[1]), error = function(e) NA_real_)
  f1 <- tryCatch(thr(interval[2]), error = function(e) NA_real_)
  
  if (!is.finite(f0) || !is.finite(f1)) return(NA_real_)
  if (f0 == 0) return(interval[1])
  if (f1 == 0) return(interval[2])
  if (sign(f0) == sign(f1)) return(NA_real_)
  
  tryCatch(
    uniroot(thr, interval = interval, tol = tol)$root,
    error = function(e) NA_real_
  )
}


# My own wrapper function to run M-DACT given sets of p-values
# with the option to compute FDR threshold via uniroot or sorted search
one_mdact <- function(all_p_df, pm, py, fdr, diagnostic_root = TRUE) {
  p_m <- setNames(all_p_df[[pm]], all_p_df[["microbe"]])
  p_y <- setNames(all_p_df[[py]], all_p_df[["microbe"]])
  
  keep <- is.finite(p_m) & is.finite(p_y) &
    p_m >= 0 & p_m <= 1 &
    p_y >= 0 & p_y <= 1
  
  if (any(!keep)) {
    warning(sum(!keep), " mediators dropped because p-values were missing/non-finite/out of [0,1].")
  }
  
  p_m <- p_m[keep]
  p_y <- p_y[keep]
  
  p_matrix <- cbind(p_m, p_y)
  
  estws <- tryCatch(
    HDMT::null_estimation(p_matrix),
    error = function(e) {
      warning(sprintf("null_estimation() failed: %s", e$message))
      return(NULL)
    }
  )
  
  if (is.null(estws)) {
    return(tibble(
      exposure = pm,
      mediator = names(p_m),
      outcome = py,
      p_m = p_m,
      p_y = p_y,
      p_composite = NA_real_,
      dact_stat = NA_real_,
      dact_thr = NA_real_,
      dact_thr_discrete = NA_real_,
      dact_thr_root = NA_real_,
      selected = FALSE,
      hit_rule = NA_character_,
      fdr_target = fdr
    ))
  }
  
  p_composite <- tryCatch(
    MDACT_pvalues(p_m, p_y, estws),
    error = function(e) {
      warning(sprintf("MDACT_pvalues() failed: %s", e$message))
      rep(NA_real_, length(p_m))
    }
  )
  
  result <- balancing_DACT_control_DR_adjust(
    p.M = p_m,
    p.Y = p_y,
    estws = estws,
    significance_upper = fdr,
    control.method = "FDR",
    diagnostic_root = diagnostic_root
  )
  
  tibble(
    exposure = pm,
    mediator = names(p_m),
    outcome = py,
    p_m = as.numeric(p_m),
    p_y = as.numeric(p_y),
    p_composite = as.numeric(p_composite),
    dact_stat = as.numeric(result$dact_list),
    dact_thr = result$dact_thr,
    dact_thr_discrete = result$dact_thr_discrete,
    dact_thr_root = result$dact_thr_root,
    selected = seq_along(p_m) %in% result$hits,
    hit_rule = result$hit_rule,
    fdr_target = fdr
  )
}


# Summary wrapper for one_mdact
one_mdact_thr <- function(all_p_df, pm, py, fdr, diagnostic_root = TRUE) {
  one_mdact(
    all_p_df = all_p_df,
    pm = pm,
    py = py,
    fdr = fdr,
    diagnostic_root = diagnostic_root
  ) %>%
    distinct(
      exposure, outcome, fdr_target,
      dact_thr, dact_thr_discrete, dact_thr_root,
      hit_rule
    )
}


# Prep bootstrap data splits
prep_group_bootstrap <- function(data, group_var) {
  split_data  <- split(data, data[[group_var]])
  group_sizes <- vapply(split_data, nrow, integer(1))
  list(split_data = split_data, group_sizes = group_sizes)
}

# Stratified resampling for bootstraps
resample_groups <- function(split_data, group_sizes) {
  bind_rows(Map(
    function(d, n) d[sample.int(n, replace = TRUE), , drop = FALSE],
    split_data,
    group_sizes
  ))
}

# stratified bootstrap using lm()
strat_boot_lm <- function(data, exposure, mediator, outcome, n_boot = 1000) {
  
  # Create 4-level treatment cell for stratification
  data$treat_cell <- interaction(data$blast, data$vns, drop = TRUE)
  
  # Pre-split once
  boot_prep   <- prep_group_bootstrap(data, "treat_cell")
  split_data  <- boot_prep$split_data
  group_sizes <- boot_prep$group_sizes
  
  # Model formulas
  form_m <- as.formula(
    paste0(mediator, " ~ blast + vns + blast_and_vns + Round_main")
  )
  form_y <- as.formula(
    paste0(outcome, " ~ ", mediator,
           " + blast + vns + blast_and_vns + Round_main")
  )
  
  # Storage matrix
  boot_mat <- matrix(NA_real_, n_boot, 6)
  colnames(boot_mat) <- c(
    "ACME_cond0", "ACME_cond1",
    "ADE_cond0",  "ADE_cond1",
    "Total_cond0","Total_cond1"
  )
  
  for (i in seq_len(n_boot)) {
    
    boot_data <- resample_groups(split_data, group_sizes)
    
    fit_m <- lm(form_m, boot_data)
    fit_y <- lm(form_y, boot_data)
    
    fe_m <- coef(fit_m)
    fe_y <- coef(fit_y)
    
    # Check required coefficients exist
    needed_m <- c(exposure, "blast_and_vns")
    needed_y <- c(exposure, mediator, "blast_and_vns")
    
    if (!all(needed_m %in% names(fe_m)) ||
        !all(needed_y %in% names(fe_y))) {
      next
    }
    
    b <- fe_y[mediator]
    
    a_cond0 <- fe_m[exposure]
    a_cond1 <- fe_m[exposure] + fe_m["blast_and_vns"]
    
    c_cond0 <- fe_y[exposure]
    c_cond1 <- fe_y[exposure] + fe_y["blast_and_vns"]
    
    ACME_0  <- a_cond0 * b
    ACME_1  <- a_cond1 * b
    
    ADE_0   <- c_cond0
    ADE_1   <- c_cond1
    
    boot_mat[i, ] <- c(
      ACME_0,            # ACME cond 0
      ACME_1,
      ADE_0,             # ADE cond 0
      ADE_1,
      ACME_0 + ADE_0,    # total cond 0
      ACME_1 + ADE_1
    )
  }
  
  boot_df <- as_tibble(boot_mat) |> drop_na()
  
  boot_df |>
    summarise(
      across(
        everything(),
        list(
          est   = mean,
          lower = ~ quantile(.x, 0.025),
          upper = ~ quantile(.x, 0.975)
        ),
        .names = "{.col}_{.fn}"
      )
    ) |>
    mutate(
      exposure = exposure,
      mediator = mediator,
      outcome  = outcome
    )
}

# stratified bootstrap using lmer()
strat_boot_lmer <- function(data, exposure, mediator, outcome, n_boot = 1000) {
  # Create 4-level treatment cell for stratification
  data$treat_cell <- interaction(data$blast, data$vns, drop = TRUE)
  data$strat <- interaction(data$treat_cell, data$Round_main, drop = TRUE)
  
  # Pre-split once
  boot_prep   <- prep_group_bootstrap(data, "strat")  # need to split on BOTH "treat_cell" AND "batch"
  split_data  <- boot_prep$split_data
  group_sizes <- boot_prep$group_sizes
  
  # formulas
  form_m <- as.formula(paste0(mediator, " ~ blast + vns + blast_and_vns + (1 | Round_main)"))
  form_y <- as.formula(paste0(outcome, " ~ ", mediator, " + blast + vns + blast_and_vns + (1 | Round_main)"))
  
  # Storage matrix
  boot_mat <- matrix(NA_real_, n_boot, 6)
  colnames(boot_mat) <- c(
    "ACME_cond0", "ACME_cond1",
    "ADE_cond0",  "ADE_cond1",
    "Total_cond0","Total_cond1"
  )
  
  for (i in seq_len(n_boot)) {
    
    boot_data <- resample_groups(split_data, group_sizes) # need to resample within "treat_cell" and "batch" combinations 
    
    fit_m <- lme4::lmer(form_m, boot_data, REML = FALSE)
    fit_y <- lme4::lmer(form_y, boot_data, REML = FALSE)
    
    fe_m <- lme4::fixef(fit_m)
    fe_y <- lme4::fixef(fit_y)
    
    # Check required coefficients exist
    needed_m <- c(exposure, "blast_and_vns")
    needed_y <- c(exposure, mediator, "blast_and_vns")
    
    if (!all(needed_m %in% names(fe_m)) ||
        !all(needed_y %in% names(fe_y))) {
      next
    }
    
    b <- fe_y[mediator]
    
    a_cond0 <- fe_m[exposure]
    a_cond1 <- fe_m[exposure] + fe_m["blast_and_vns"]
    
    c_cond0 <- fe_y[exposure]
    c_cond1 <- fe_y[exposure] + fe_y["blast_and_vns"]
    
    ACME_0  <- a_cond0 * b
    ACME_1  <- a_cond1 * b
    
    ADE_0   <- c_cond0
    ADE_1   <- c_cond1
    
    boot_mat[i, ] <- c(
      ACME_0,            # ACME cond 0
      ACME_1,
      ADE_0,             # ADE cond 0
      ADE_1,
      ACME_0 + ADE_0,    # total cond 0
      ACME_1 + ADE_1
    )
  }
  
  boot_df <- as_tibble(boot_mat) |> drop_na()
  
  boot_df |>
    summarise(
      across(
        everything(),
        list(
          est   = mean,
          lower = ~ quantile(.x, 0.025),
          upper = ~ quantile(.x, 0.975)
        ),
        .names = "{.col}_{.fn}"
      )
    ) |>
    mutate(
      exposure = exposure,
      mediator = mediator,
      outcome  = outcome
    )
}
