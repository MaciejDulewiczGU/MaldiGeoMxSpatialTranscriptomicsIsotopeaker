# MALDI Spectra Processing Module (with "private" helpers)
#
# Helper functions have a '.' prefix. User-facing API functions do not.
#
# -----------------------------------------------------------------------------
# DEPENDENCIES (import at script or package level)
# -----------------------------------------------------------------------------
{
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(xlsx)      # for single-sheet .xlsx
  library(openxlsx)  # for multi-sheet .xlsx (optional)
  library(signal)    # for sgolay smoothing
  library(ggplot2)   # for analysis/plotting
  library(minpack.lm) # for nlsLM fit
}
# -----------------------------------------------------------------------------

# ----------------------------- 1. DATA INGEST --------------------------------

# "Private" helper
.load_file <- function(file) {
  data <- read.table(file, header = FALSE)
  colnames(data) <- c("mz_value", "intensity")
  return(data)
}

# User-facing
load_spectra_files <- function(folder_path, pattern = "\\.txt$") {
  files <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  names <- basename(files) %>% sub(pattern, "", .)
  data  <- lapply(files, .load_file)
  names(data) <- names
  data
}

.create_mass_spectra_from_list_custom <- function(data_list, file_names, search_spectra, range_for_spectra) {
  mass_spectra_list <- list()
  for (i in seq_along(data_list)) {
    file_name <- file_names[i]
    filtered_data <- data_list[[i]] %>%
      dplyr::filter(mz_value >= (search_spectra - range_for_spectra) & 
                      mz_value <= (search_spectra + range_for_spectra))
    if (nrow(filtered_data) == 0) {
      warning(paste("There is not files and spectra:", file_name, "between range m/z."))
      next
    }
    spectrum <- .create_spectrum(
      mz_vec = filtered_data$mz_value,
      intensity_vec = filtered_data$intensity,
      file_name = file_name
    )
    mass_spectra_list[[file_name]] <- spectrum
  }
  return(mass_spectra_list)
}

.create_spectrum <- function(mz_vec, intensity_vec, file_name="") {
  list(
    mz = mz_vec,
    intensity = intensity_vec,
    metadata = list(
      filename = file_name,
      baseline = NULL,
      noise = NULL
    )
  )
}

# User-facing (shortcut)
create_spectra <- function(data_list, file_names, target_mz, range = 25) {
  .create_mass_spectra_from_list_custom(
    data_list         = data_list,
    file_names        = file_names,
    search_spectra    = target_mz,
    range_for_spectra = range
  )
}

# ---------------------------- 2. PRE-PROCESSING ------------------------------

.scale_single_spectrum <- function(spectrum) {
  max_intensity <- max(spectrum$intensity, na.rm = TRUE)
  if (max_intensity > 0) {
    spectrum$intensity <- (spectrum$intensity / max_intensity) * 100
  }
  return(spectrum)
}

# User-facing
scale_spectra_pct <- function(spectra) {
  lapply(spectra, .scale_single_spectrum)
}

# Smoothing
.smooth_sgolay <- function(y, polyOrder=3, windowSize=11) {
  if(windowSize <= polyOrder) stop("Savitzky-Golay error: window size must exceed polynomial order.")
  signal::sgolayfilt(y, p=polyOrder, n=windowSize)
}

.smooth_spectra <- function(spectra_list, method = "sgolay", ...) {
  lapply(spectra_list, function(spectrum) {
    spectrum$intensity <- switch(method,
                                 "sgolay" = .smooth_sgolay(spectrum$intensity, ...),
                                 spectrum$intensity  # default: no smoothing
    )
    spectrum
  })
}

# User-facing
smooth_spectra_sgolay <- function(spectra, polyOrder = 3, windowSize = 49) {
  .smooth_spectra(
    spectra    = spectra,
    method     = "sgolay",
    polyOrder  = polyOrder,
    windowSize = windowSize
  )
}

# Baseline
.remove_baseline <- function(spectrum, baseline, valid = NULL, clip_negative = TRUE) {
  if (is.null(valid)) {
    corrected <- spectrum$intensity - baseline
    if (clip_negative) corrected[corrected < 0] <- 0
    spectrum$intensity <- corrected
  } else {
    corrected <- spectrum$intensity[valid] - baseline
    if (clip_negative) corrected[corrected < 0] <- 0
    spectrum$intensity[valid] <- corrected
  }
  spectrum
}

.remove_baseline_spectra <- function(spectra_list, method = "typical", clip_negative = TRUE, ...) {
  filter_spectrum <- function(mz, intensity) {
    valid <- !is.na(mz) & !is.na(intensity) & is.finite(mz) & is.finite(intensity)
    list(mz = mz[valid], intensity = intensity[valid], valid = valid)
  }
  lapply(spectra_list, function(spectrum) {
    if (is.null(spectrum$mz) || is.null(spectrum$intensity)) {
      warning("Spectrum does not contain mz or intensity. Returning unchanged.")
      return(spectrum)
    }
    mz <- spectrum$mz; intensity <- spectrum$intensity
    f <- filter_spectrum(mz, intensity)
    mz <- f$mz; intensity <- f$intensity; valid <- f$valid
    if (length(mz) < 3) {
      warning("Too few valid points (< 3). Returning unchanged.")
      return(spectrum)
    }
    baseline <- switch(tolower(method),
                       als = {
                         b_mat <- .als_baseline(matrix(intensity, nrow = 1), ...)
                         if (is.matrix(b_mat)) b_mat[1, ] else as.numeric(b_mat)
                       },
                       stop("Unknown baseline method: ", method)
    )
    .remove_baseline(spectrum, baseline, valid = valid, clip_negative = clip_negative)
  })
}


.als_baseline <- function(y, lambda = 6, p = 0.05, maxit = 20, tol = 1e-6) {
  y <- as.matrix(y)
  n <- ncol(y); m <- nrow(y)
  D <- diff(diag(n), differences = 2)
  DtD <- crossprod(D)
  lam <- 10^lambda
  Baseline <- matrix(0.0, m, n)
  for (row in seq_len(m)) {
    yrow <- y[row, ]
    w <- rep(1, n)
    z_old <- yrow
    for (iter in seq_len(maxit)) {
      z <- solve(diag(w) + lam * DtD, w * yrow)
      r <- yrow - z
      w <- ifelse(r > 0, p, 1 - p)
      if (sqrt(mean((z - z_old)^2)) < tol) break
      z_old <- z
    }
    Baseline[row, ] <- z
  }
  Baseline
}

# User-facing
remove_baseline_als <- function(spectra, lambda = 5, p = 0.025, maxit = 100) {
  .remove_baseline_spectra(
    spectra = spectra,
    method  = "als",
    lambda  = lambda,
    p       = p,
    maxit   = maxit
  )
}

# ----------------------------- 3. ANALYSIS -----------------------------------

.analyze_spectra_simple <- function(
    spectra,
    peak_window   = 18,
    threshold_rel = 1e-4,
    scale_percent = FALSE
) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  plots_list   <- list()
  results_list <- list()
  asym_gaussian <- function(x, A, mu, sigma_left, sigma_right) {
    ifelse(x < mu,
           A * exp(-0.5 * ((x - mu) / sigma_left)^2),
           A * exp(-0.5 * ((x - mu) / sigma_right)^2))
  }
  for (s in spectra) {
    filename <- s$metadata$filename %||% "spectrum"
    mz  <- s$mz
    raw_int <- s$intensity
    if (length(mz) == 0 || length(raw_int) == 0) next
    if (scale_percent) {
      scale_fac <- max(raw_int)
      int <- raw_int / scale_fac * 100
      ylab <- "Intensity [%]"
    } else {
      int <- raw_int
      ylab <- "Intensity"
    }
    df <- data.frame(mz = mz, intensity = int, row.names = NULL)
    peaks <- df %>% dplyr::filter(intensity > lag(intensity, default = first(intensity)),
                                  intensity > lead(intensity, default = last(intensity)))
    main_peak <- peaks %>% dplyr::filter(intensity == max(intensity)) %>% slice(1)
    if (nrow(main_peak) == 0) next
    fit_rng <- c(main_peak$mz - peak_window, main_peak$mz + peak_window)
    fit_df  <- df %>% dplyr::filter(mz >= fit_rng[1], mz <= fit_rng[2]) %>% mutate(y = pmax(intensity, 0))
    start <- list(A = max(fit_df$y), mu = main_peak$mz, sigma_left = 0.11, sigma_right = 0.11)
    fit_mod <- tryCatch(nlsLM(y ~ asym_gaussian(mz, A, mu, sigma_left, sigma_right),
                              data = fit_df, start = start,
                              lower = c(0, main_peak$mz - 1, 0, 0),
                              upper = c(Inf, main_peak$mz + 1, Inf, Inf),
                              control = nls.lm.control(maxiter = 500)), error = function(e) NULL)
    if (is.null(fit_mod)) next
    pars <- coef(fit_mod); A_fit <- pars["A"]; mu_fit <- pars["mu"]
    pred_y <- predict(fit_mod)
    half_max <- A_fit / 2
    fit_df$pred <- pred_y
    left_idx  <- which(fit_df$mz < mu_fit & fit_df$pred <= half_max)
    right_idx <- which(fit_df$mz > mu_fit & fit_df$pred <= half_max)
    if (length(left_idx) == 0 || length(right_idx) == 0) next
    fwhm_left_mz <- fit_df$mz[max(left_idx)]
    fwhm_right_mz <- fit_df$mz[min(right_idx)]
    FWHM_total <- fwhm_right_mz - fwhm_left_mz
    thresh <- max(int) * threshold_rel
    above <- which(df$intensity >= thresh)
    if (!length(above)) next
    cl <- cumsum(c(1, diff(above) > 1))
    mu_idx <- which.min(abs(df$mz - mu_fit)); cl_mu <- cl[match(mu_idx, above)]
    cluster <- above[cl == cl_mu]
    peak_start <- df$mz[cluster[1]]; peak_end <- df$mz[cluster[length(cluster)]]
    peak_width <- peak_end - peak_start
    res <- fit_df$y - pred_y; R2 <- 1 - sum(res^2) / sum((fit_df$y - mean(fit_df$y))^2)
    mz_diff <- abs(main_peak$mz - mu_fit)
    results_list[[filename]] <- data.frame(
      peak_intensity = main_peak$intensity * if (scale_percent) scale_fac/100 else 1,
      mz_difference = round(mz_diff, 3),
      filename = filename,
      FWHM_total = round(FWHM_total, 3),
      peak_start = round(peak_start, 3),
      peak_end = round(peak_end, 3),
      peak_width = round(peak_width, 3),
      R2 = round(R2, 3),
      row.names = NULL
    )
    fit_curve <- data.frame(mz = fit_df$mz, y = pred_y, row.names = NULL)
    fwhm_pts  <- data.frame(mz = c(fwhm_left_mz, fwhm_right_mz), y = half_max, row.names = NULL)
    p <- ggplot(df, aes(mz, intensity)) +
      geom_line(colour = "blue") +
      geom_line(data = fit_curve, aes(mz, y), colour = "red") +
      geom_hline(yintercept = 0, colour = "brown", linetype = "dotted") +
      geom_hline(yintercept = half_max, colour = "red", linetype = "dashed") +
      geom_vline(xintercept = mu_fit, colour = "purple", linetype = "dashed") +
      geom_vline(xintercept = peak_start, colour = "darkgreen", linetype = "dashed") +
      geom_vline(xintercept = peak_end, colour = "darkgreen", linetype = "dashed") +
      geom_vline(xintercept = fwhm_left_mz,  colour = "darkorange", linetype = "dotdash") +
      geom_vline(xintercept = fwhm_right_mz, colour = "darkorange", linetype = "dotdash") +
      geom_point(data = fwhm_pts, aes(mz, y), colour = "black", size = 2) +
      annotate("segment", x = fwhm_left_mz, xend = fwhm_right_mz, y = half_max, yend = half_max, colour = "darkblue") +
      annotate("text", x = fwhm_right_mz, y = half_max*1.05, label = paste("FWHM =", round(FWHM_total, 3)), hjust = -0.3, colour = "darkblue", size = 3) +
      labs(title = paste("Spectrum –", filename), x = "m/z", y = ylab) +
      theme_minimal()
    plots_list[[filename]] <- p
  }
  list(plots = plots_list, df_results = results_list)
}

# User-facing
analyze_spectra_metrics <- function(spectra, scale_percent = TRUE) {
  .analyze_spectra_simple(
    spectra       = spectra,
    scale_percent = scale_percent
  )
}

# ----------------------------- 4. UTILITIES ----------------------------------

combine_metrics <- function(metrics_list, name_col = "name") {
  dplyr::bind_rows(metrics_list, .id = name_col)
}

plot_spectra_grid <- function(spectra, ncol = 2, title_prefix = "Spectrum") {
  plots <- lapply(seq_along(spectra), function(i) {
    if (inherits(spectra[[i]], c("ggplot", "plotly"))) {
      spectra[[i]]
    } else {
      plot_spectrum(spectra[[i]], main = paste0(title_prefix, " ", i))
    }
  })
  do.call(grid.arrange, c(plots, ncol = ncol))
}

export_metrics_to_excel <- function(metrics, file_name = "results.xlsx") {
  if (is.list(metrics) && !is.data.frame(metrics)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' needed for multiple sheets. Install it or pass a data frame instead.")
    }
    openxlsx::write.xlsx(metrics, file_name)
  } else if (is.data.frame(metrics)) {
    xlsx::write.xlsx(metrics, file_name)
  } else {
    stop("'metrics' must be a data frame or a named list of data frames.")
  }
}

# ----------------------------- 5. EXAMPLES -----------------------------------
# Example usage:
# raw   <- load_spectra_files("PATH/TO/SPECTRA/FOLDER")
# spec  <- create_spectra(raw, names(raw), target_mz = 4525)
# spec  <- scale_spectra_pct(spec)
# spec  <- smooth_spectra_sgolay(spec)
# spec  <- remove_baseline_als(spec)
# res   <- analyze_spectra_metrics(spec)
# export_metrics_to_excel(res$df_results, "multi_sheet_results.xlsx")
