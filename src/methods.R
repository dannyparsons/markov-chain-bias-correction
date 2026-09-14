# Threshold adjustment
thresh <- function(obs, est, obs_thresh = 0.85) {
  id <- !is.na(obs) & !is.na(est)
  p_obs <- mean(obs[id] > obs_thresh, na.rm = TRUE)
  quantile(est[id], 1 - p_obs)[[1]]
}

s_loci <- function(obs, est, obs_thresh = 0.85, est_thresh = NULL) {
  if (is.null(est_thresh)) est_thresh <- thresh(obs, est, obs_thresh)
  id <- !is.na(obs) & !is.na(est)
  s_obs <- mean(obs[id & obs > obs_thresh]) - obs_thresh
  s_est <- mean(est[id & est > est_thresh]) - est_thresh
  s <- s_obs / s_est
  s
}

# Local Intensity Scaling
loci <- function(obs, est, obs_thresh = 0.85) {
  est_thresh <- thresh(obs, est, obs_thresh)
  s <- s_loci(obs, est, obs_thresh, est_thresh)
  est_loci <- obs_thresh + s * (est - est_thresh)
  est_loci <- pmax(est_loci, 0)
  est_loci
}

safe_coef <- function(fit, name) {
  if (is.null(fit)) return(NA)
  coef(fit)[[name]]
}

# Quantile Mapping
fit_gamma <- function(obs, obs_thresh = 0.85, 
                      method = c("mle", "mme", "qme", "mge", "mse")) {
  obs <- obs[obs > obs_thresh] - obs_thresh
  obs <- as.numeric(na.omit(obs))
  if(length(obs) > 5) fitdistrplus::fitdist(obs, distr = "gamma", method = method)
  else NULL
}

fit_empirical <- function(obs, obs_thresh = 0.85) {
  obs <- obs[obs > obs_thresh]
  if (sum(!is.na(obs)) == 0) return(NULL)
  ecdf(obs)
}

qm_gamma <- function(est, est_thresh, est_shape, est_rate, obs_shape, obs_rate,
                     obs_thresh) {
  if_else(est <= est_thresh, 0,
          qgamma(pgamma(est - est_thresh, shape = est_shape,
                        rate = est_rate),
                 shape = obs_shape, rate = obs_rate)
          + obs_thresh)
  
}

qm_empirical <- function(est, est_thresh, obs, est_ecdf) {
  if (is.null(est_ecdf)) return(NA)
  if_else(est <= est_thresh, 0,
          quantile(obs, est_ecdf(est), na.rm = TRUE, names = FALSE))
}

quantile_mapping <- function(obs, est, obs_thresh = 0.85, 
                             qm_method = c("gamma", "empirical"), 
                             gamma_method = c("mle", "mme", "qme", 
                                              "mge", "mse")) {
  est_thresh <- thresh(obs, est, obs_thresh)
  qm_method <- match.arg(qm_method)
  if (qm_method == "gamma") {
    obs_gamma <- fit_gamma(obs, obs_thresh, method = gamma_method)
    est_gamma <- fit_gamma(est, est_thresh, method = gamma_method)
    obs_shape <- obs_gamma[["estimate"]][["shape"]]
    obs_rate <- obs_gamma[["estimate"]][["rate"]]
    est_shape <- est_gamma[["estimate"]][["shape"]]
    est_rate <- est_gamma[["estimate"]][["rate"]]
    est_qm <- if_else(est <= est_thresh, 0,
                      qgamma(pgamma(est - est_thresh, shape = est_shape, 
                                    rate = est_rate),
                             shape = obs_shape, rate = obs_rate) + obs_thresh)
  } else if (qm_method == "empirical") {
    obs_rain <- obs[obs > obs_thresh]
    est_ecdf <- fit_empirical(est, est_thresh)
    est_qm <- if_else(est <= est_thresh, 0,
                      as.numeric(quantile(obs_rain, est_ecdf(est), 
                                          na.rm = TRUE)))
  }
  return(est_qm)
}

# Calculate logical wet/dry on x using t_w and t_d
# y0 is a logical wet/dry value on day before start of x series, if available
conditional_wd <- function(x, t_w, t_d, t0, y0) {
  n <- length(x)
  y <- logical(n)
  # TODO Split out the case of thresholds as a column and single number
  # Probably inefficient to do this:
  if (length(t_w) == 1) t_w <- rep(t_w, n)
  if (length(t_d) == 1) t_d <- rep(t_d, n)
  if (length(t0) == 1) t0 <- rep(t0, n)
  # Set first value based only on y0 (if available) and x[1]
  if (!missing(y0) && !is.na(y0)) {
    y[1] <- if(y0) x[1] > t_w[1] else x[1] > t_d[1]
  } else y[1] <- x[1] > t0[1]
  for (i in 2:n) {
    # If rain on previous day is NA
    # then use today's rain only with overall threshold
    if (is.na(y[i - 1])) y[i] <- x[i] > t0[i]
    else y[i] <- if(y[i - 1]) x[i] > t_w[i] else x[i] > t_d[i]
  }
  y
}

markov_thresholds <- function(data, obs_col = "obs", est_col = "est", 
                              date_col = "date", season_col, station_col,
                              obs_thr = 0.85, tol = 1e-2, max_it = 20, 
                              n_conv = 5, damping = 0.4, loci = TRUE,
                              qm_empirical = TRUE, qm_gamma = TRUE,
                              gamma_method = c("mle", "mme", "qme", "mge", 
                                               "mse")) {
  
  data[["year"]] <- lubridate::year(data[[date_col]])
  if (missing(season_col)) {
    data[["season"]] <- factor(lubridate::month(data[[date_col]]))
  } else {
    data <- data %>% rename(season = all_of(season_col))
  }
  season_col <- "season"
  data[["day"]] <- lubridate::day(data[[date_col]])
  
  # Insert NAs where obs/est is NA to ensure that only days with both
  # obs and est values are included in the calibration
  data[[obs_col]] <- if_else(is.na(data[[est_col]]), NA_real_, data[[obs_col]])
  data[[est_col]] <- if_else(is.na(data[[obs_col]]), NA_real_, data[[est_col]])

  obs <- data[[obs_col]]
  est <- data[[est_col]]
  # logical wet/dry from observations
  # Calculate lags before subsetting to ensure correct previous days
  data[["obs_wd"]] <- obs > obs_thr
  data[["obs_wd_prev"]] <- dplyr::lag(data[["obs_wd"]])
  
  results <- tibble(
    station = character(),
    season = character(),
    t0 = numeric(),
    t_w = numeric(),
    t_d = numeric(),
    p0_obs = numeric(),
    p0_est = numeric(),
    p_w_obs = numeric(),
    p_w_est = numeric(),
    p_d_obs = numeric(),
    p_d_est = numeric(),
    converged = character(),
    iterations = integer(),
    n_days = integer(),
    s_all = numeric(),
    s_wet = numeric(),
    s_dry = numeric(),
    obs_all = list(),
    obs_wet = list(),
    obs_dry = list(),
    ecdf_est_all = list(),
    ecdf_est_wet = list(),
    ecdf_est_dry = list(),
    gamma_obs_all = list(),
    gamma_est_all = list(),
    gamma_obs_wet = list(),
    gamma_obs_dry = list(),
    gamma_est_wet = list(),
    gamma_est_dry = list()
  )
  if (!missing(station_col)) {
    stations <- unique(data[[station_col]])
    data <- data %>% rename(station = all_of(station_col))
    station_col <- "station"
  } else stations <- NA
  
  for (s in stations) {
    for (m in levels(data[[season_col]])) {
      if (!is.na(s)) data_m <- dplyr::filter(data, station == s & season == m)
      else data_m <- dplyr::filter(data, season == m)
      est_m <- data_m[[est_col]]
      obs_wd <- data_m[["obs_wd"]]
      obs_wd_prev <- data_m[["obs_wd_prev"]]
      
      # Target probabilities
      p_obs <- mean(obs_wd, na.rm = TRUE)
      p_w_obs <- mean(obs_wd[obs_wd_prev], na.rm = TRUE)
      p_d_obs <- mean(obs_wd[!obs_wd_prev], na.rm = TRUE)
      
      # Initial threshold based on p_obs
      t0 <- quantile(est_m, probs = 1 - p_obs, na.rm = TRUE, 
                     names = FALSE)
      t_d <- t0
      t_w <- t0
      # State of convergence, either: "tol", "neg diff" or "no"
      converged <- "no"
      n_same <- 0
      # Iterate to converge to target probabilities
      for (i in 1:max_it) {
        res <- data_m %>%
          group_by(year, season) %>%
          group_modify(~{
            y0 <- .x$obs_wd_prev[1]
            if (is.na(y0)) y0 <- FALSE
            y  <- conditional_wd(.x[[est_col]], t_w, t_d, t0)
            # dplyr::lag() fails if y is all NAs so do it manually
            if(all(is.na(y))) lag_y <- c(y0, rep(NA, length(y) - 1))
            else lag_y <- dplyr::lag(y, default = y0)
            tibble(.x, est_wd = y, est_wd_prev = lag_y)
          }) %>%
          ungroup()
        # Update probabilities
        est_wd <- res$est_wd
        est_wd_prev <- res$est_wd_prev
        p_est_new <- mean(est_wd, na.rm = TRUE)
        p_w_est_new <- mean(est_wd[est_wd_prev], na.rm = TRUE)
        p_d_est_new <- mean(est_wd[!est_wd_prev], na.rm = TRUE)
        # Calculate difference between previous probabilities 
        # if not first iteration
        if (i > 1) {
          diff_p <- abs(p_est_new - p_est)
          diff_p_w <- abs(p_w_est_new - p_w_est)
          diff_p_d <- abs(p_d_est_new - p_d_est)
        }
        p_est <- p_est_new
        p_w_est <- p_w_est_new
        p_d_est <- p_d_est_new
        # Compare with target probabilities and stop if sufficiently converged
        tol_w <- 1 / max(sum(res$est_wd_prev, na.rm = TRUE), 1)
        tol_d <- 1 / max(sum(!res$est_wd_prev, na.rm = TRUE), 1)
        # Tolerance to check if probabilities have changed from previous iteration
        tol_diff_w <- tol_w / 2
        tol_diff_d <- tol_d / 2
        if (i > 1 && diff_p_w < tol_diff_w && diff_p_d < tol_diff_d) {
          n_same <- n_same + 1
        }
        if (is.na(p_w_est) || is.na(p_w_obs) || 
            is.na(p_d_est) || is.na(p_d_obs)) break
        # Stop if probabilities are within tolerance of targets
        else if ((abs(p_w_est - p_w_obs) < max(tol, tol_w) && 
                  abs(p_d_est - p_d_obs) < max(tol, tol_d))) {
          converged <- "within tol"
          break
          # Stop if difference in probabilities is within negligible difference
          # after n_conv iterations
        } else if (n_same >= n_conv) {
          converged <- "neg diff"
          break
        }
        
        # Update thresholds based on current wet/dry and target probabilities
        if (i != max_it) {
          t_w_new <- quantile(res[[est_col]][res$est_wd_prev],
                              probs = 1 - p_w_obs, na.rm = TRUE, names = FALSE)
          t_w <- (1 - damping) * t_w + damping * t_w_new
          t_d_new <- quantile(res[[est_col]][!res$est_wd_prev],
                              probs = 1 - p_d_obs, na.rm = TRUE, names = FALSE)
          t_d <- (1 - damping) * t_d + damping * t_d_new
        }
      }
      s_wet <- NA
      s_dry <- NA
      obs_prev_wet <- data_m[[obs_col]][obs_wd_prev]
      obs_prev_dry <- data_m[[obs_col]][!obs_wd_prev]
      est_prev_wet <- data_m[[est_col]][est_wd_prev]
      est_prev_dry <- data_m[[est_col]][!est_wd_prev]
      if (loci) {
        s_obs_all <- mean(data_m[[obs_col]][data_m[[obs_col]] > obs_thr], 
                          na.rm = TRUE) - obs_thr
        s_est_all <- mean(data_m[[est_col]][data_m[[est_col]] > t0], 
                          na.rm = TRUE) - t0
        s_all <- s_obs_all / s_est_all
        
        s_obs_wet <- mean(obs_prev_wet[obs_prev_wet > obs_thr], 
                          na.rm = TRUE) - obs_thr
        s_obs_dry <- mean(obs_prev_dry[obs_prev_dry > obs_thr], 
                          na.rm = TRUE) - obs_thr
        s_est_wet <- mean(est_prev_wet[est_prev_wet > t_w], 
                          na.rm = TRUE) - t_w
        s_est_dry <- mean(est_prev_dry[est_prev_dry > t_d], 
                          na.rm = TRUE) - t_d
        
        s_wet <- s_obs_wet / s_est_wet
        s_dry <- s_obs_dry / s_est_dry
      }
      if (qm_empirical) {
        obs_all <- data_m[[obs_col]][data_m[[obs_col]] > obs_thr]
        obs_wet <- obs_prev_wet[obs_prev_wet > obs_thr]
        obs_dry <- obs_prev_dry[obs_prev_dry > obs_thr]
        ecdf_est_all <- fit_empirical(data_m[[est_col]], t0)
        ecdf_est_wet <- fit_empirical(est_prev_wet, t_w)
        ecdf_est_dry <- fit_empirical(est_prev_dry, t_d)
      }
      if (qm_gamma) {
        gamma_obs_all <- fit_gamma(data_m[[obs_col]], obs_thresh = obs_thr, 
                                   method = gamma_method)
        gamma_est_all <- fit_gamma(data_m[[est_col]], obs_thresh = t0, 
                                   method = gamma_method)
        
        gamma_obs_wet <- fit_gamma(obs_prev_wet, obs_thresh = obs_thr, 
                                   method = gamma_method)
        gamma_obs_dry <- fit_gamma(obs_prev_dry, obs_thresh = obs_thr, 
                                   method = gamma_method)
        
        gamma_est_wet <- fit_gamma(est_prev_wet, obs_thresh = t_w, 
                                   method = gamma_method)
        gamma_est_dry <- fit_gamma(est_prev_dry, obs_thresh = t_d, 
                                   method = gamma_method)
      }
      
      results <- results %>% 
        tibble::add_row(station = s, season = m, t0 = t0, t_w = t_w, t_d = t_d,
                        p0_obs = p_obs, p0_est = p_est, p_w_obs = p_w_obs, 
                        p_w_est = p_w_est, p_d_obs = p_d_obs, 
                        p_d_est = p_d_est, converged = converged, 
                        iterations = i, n_days = nrow(data_m), 
                        s_all = s_all, s_wet = s_wet, s_dry = s_dry,
                        obs_all = list(obs_all), obs_wet = list(obs_wet), 
                        obs_dry = list(obs_dry),
                        ecdf_est_all = list(ecdf_est_all),
                        ecdf_est_wet = list(ecdf_est_wet),
                        ecdf_est_dry = list(ecdf_est_dry),
                        gamma_obs_all = list(gamma_obs_all),
                        gamma_est_all = list(gamma_est_all),
                        gamma_obs_wet = list(gamma_obs_wet),
                        gamma_obs_dry = list(gamma_obs_dry),
                        gamma_est_wet = list(gamma_est_wet),
                        gamma_est_dry = list(gamma_est_dry))
    }
  }
  if (all(is.na(stations))) results$station <- NULL
  if (all(is.na(results$s_all))) results$s_all <- NULL
  if (all(is.na(results$s_wet))) results$s_wet <- NULL
  if (all(is.na(results$s_dry))) results$s_dry <- NULL
  results
}

markov_loci <- function(data, obs_col = "obs", est_col = "est",
                        date_col = "date", season_col, station_col,
                        obs_thr = 0.85, tol = 1e-2, max_it = 20, n_conv = 5,
                        damping = 0.4, blocks) {
  data <- data %>% 
    rename(date = all_of(date_col),
           season = all_of(season_col),
           station = all_of(station_col))
  result_list <- list()
  stations <- unique(data$station)
  c <- 1
  for (st in stations) {
    data_st <- data %>% filter(station == st)
    for (b in 1:(length(blocks) - 1)) {
      # Calibrate parameters on all data except current block, apply to current block
      # data to calibrate parameters
      data_cal <- data_st %>%
        filter(!(date >= blocks[b] & date < blocks[b + 1]))
      # data to apply bias correction to
      data_apply <- data_st %>% 
        filter(date >= blocks[b] & date < blocks[b + 1])
      
      # Calculate thresholds on calibration data
      m_thresh <- markov_thresholds(data_cal, obs_col = obs_col, est_col = est_col, 
                                    date_col = date_col, season_col = season_col, 
                                    station_col = station_col, obs_thr = obs_thr, 
                                    tol = tol, max_it = max_it, n_conv = n_conv, 
                                    damping = damping)
      data_apply <- data_apply %>%
        dplyr::left_join(m_thresh %>% filter(station == st), 
                         by = c("station", "season"))
      # Extract Empirical parameters for QM - Empirical
      obs_all <- data_apply$obs_all
      obs_wet <- data_apply$obs_wet
      obs_dry <- data_apply$obs_dry
      ecdf_est_all <- data_apply$ecdf_est_all
      ecdf_est_wet <- data_apply$ecdf_est_wet
      ecdf_est_dry <- data_apply$ecdf_est_dry
      
      # Extract Gamma parameters for QM - Gamma
      shape_obs_all <- map_dbl(data_apply$gamma_obs_all, ~ safe_coef(.x, "shape"))
      rate_obs_all <- map_dbl(data_apply$gamma_obs_all, ~ safe_coef(.x, "rate"))
      shape_est_all <- map_dbl(data_apply$gamma_est_all, ~ safe_coef(.x, "shape"))
      rate_est_all <- map_dbl(data_apply$gamma_est_all, ~ safe_coef(.x, "rate"))
      shape_obs_wet <- map_dbl(data_apply$gamma_obs_wet, ~ safe_coef(.x, "shape"))
      rate_obs_wet <- map_dbl(data_apply$gamma_obs_wet, ~ safe_coef(.x, "rate"))
      shape_obs_dry <- map_dbl(data_apply$gamma_obs_dry, ~ safe_coef(.x, "shape"))
      rate_obs_dry <- map_dbl(data_apply$gamma_obs_dry, ~ safe_coef(.x, "rate"))
      shape_est_wet <- map_dbl(data_apply$gamma_est_wet, ~ safe_coef(.x, "shape"))
      rate_est_wet <- map_dbl(data_apply$gamma_est_wet, ~ safe_coef(.x, "rate"))
      shape_est_dry <- map_dbl(data_apply$gamma_est_dry, ~ safe_coef(.x, "shape"))
      rate_est_dry <- map_dbl(data_apply$gamma_est_dry, ~ safe_coef(.x, "rate"))
      
      # Extract parameters
      n <- nrow(data_apply)
      est <- data_apply[[est_col]]
      t0 <- data_apply[["t0"]]
      t_w <- data_apply[["t_w"]]
      t_d <- data_apply[["t_d"]]
      
      # Extract LOCI parameters
      s <- data_apply[["s_all"]]
      s_wet <- data_apply[["s_wet"]]
      s_dry <- data_apply[["s_dry"]]
      
      # Standard BC methods
      # TODO These could be done via mutate, loop shouldn't be needed
      est_loci <- numeric(n)
      est_qm_empirical <- numeric(n)
      est_qm_gamma <- numeric(n)
      est_loci[1] <- if_else(est[1] < t0[1], 0, obs_thr + s[1] * (est[1] - t0[1]))
      est_qm_empirical[1] <- qm_empirical(est[1], t0[1], obs_all[[1]], ecdf_est_all[[1]])
      est_qm_gamma[1] <- qm_gamma(est[1], t0[1], shape_est_all[1], rate_est_all[1], 
                                     shape_obs_all[1], rate_obs_all[1], obs_thr)
      
      # Markov chain enhanced BC methods
      est_loci_mk <- numeric(n)
      est_loci_mk[1] <- if_else(est[1] < t0[1], 0, obs_thr + s[1] * (est[1] - t0[1]))
      est_qm_empirical_mk <- numeric(n)
      est_qm_empirical_mk[1] <- qm_empirical(est[1], t0[1], obs_all[[1]], ecdf_est_all[[1]])
      est_qm_gamma_mk <- numeric(n)
      est_qm_gamma_mk[1] <- qm_gamma(est[1], t0[1], shape_est_all[1], rate_est_all[1], 
                                  shape_obs_all[1], rate_obs_all[1], obs_thr)
      
      for (i in 2:n) {
        est_loci[i] <- if_else(est[i] < t0[i], 0, obs_thr + s[i] * (est[i] - t0[i]))
        est_qm_empirical[i] <- qm_empirical(est[i], t0[i], obs_all[[i]], ecdf_est_all[[i]])
        est_qm_gamma[i] <- qm_gamma(est[i], t0[i], shape_est_all[i], rate_est_all[i],
                                    shape_obs_all[i], rate_obs_all[i], obs_thr)
        if (is.na(est_loci_mk[i - 1])) {
          est_loci_mk[i] <- if_else(est[i] < t0[i], 0, obs_thr + s[i] * (est[i] - t0[i]))
          est_qm_empirical_mk[i] <- qm_empirical(est[i], t0[i], obs_all[[i]], ecdf_est_all[[i]])
          est_qm_gamma_mk[i] <- qm_gamma(est[i], t0[i], shape_est_all[i], rate_est_all[i], 
                                      shape_obs_all[i], rate_obs_all[i], obs_thr)
        } else if (est_loci_mk[i - 1] > obs_thr) {
          est_loci_mk[i] <- if_else(est[i] < t_w[i], 0, obs_thr + s_wet[i] * (est[i] - t_w[i]))
          est_qm_empirical_mk[i] <- qm_empirical(est[i], t_w[i], obs_wet[[i]], ecdf_est_wet[[i]])
          est_qm_gamma_mk[i] <- qm_gamma(est[i], t_w[i], shape_est_wet[i], rate_est_wet[i], 
                                      shape_obs_wet[i], rate_obs_wet[i], obs_thr)
        } else {
          est_loci_mk[i] <- if_else(est[i] < t_d[i], 0, obs_thr + s_dry[i] * (est[i] - t_d[i]))
          est_qm_empirical_mk[i] <- qm_empirical(est[i], t_d[i], obs_dry[[i]], ecdf_est_dry[[i]])
          est_qm_gamma_mk[i] <- qm_gamma(est[i], t_d[i], shape_est_dry[i], rate_est_dry[i], 
                                      shape_obs_dry[i], rate_obs_dry[i], obs_thr)
        }
      }
      data_apply$est_loci <- est_loci
      data_apply$est_loci_mk <- est_loci_mk
      data_apply$est_qm_empirical <- est_qm_empirical
      data_apply$est_qm_empirical_mk <- est_qm_empirical_mk
      data_apply$est_qm_gamma <- est_qm_gamma
      data_apply$est_qm_gamma_mk <- est_qm_gamma_mk
      result_list[[c]] <- data_apply
      c <- c + 1
    }
  }
  bind_rows(result_list)
}

## ----------------------------------------------------------------------

# Initial work on possible stochastic methods (not used)
fit_rain_prob <- function(data, obs_col = "obs", est_col = "est", 
                          date_col = "date", season_col, station_col,
                          obs_thr = 0.85) {
  data <- data %>% 
    rename(obs = all_of(obs_col),
           est = all_of(est_col)) %>%
    mutate(obs_wd = obs > obs_thr,
           obs_wd_prev = dplyr::lag(obs_wd),
           est_wd = est > obs_thr,
           est_wd_prev = dplyr::lag(est_wd)) %>%
    # needed for testing models to ensure same data used but can remove in practice
    filter(!is.na(obs_wd_prev))

  if (!missing(station_col)) {
    stations <- unique(data[[station_col]])
    data <- data %>% rename(station = all_of(station_col))
    station_col <- "station"
  } else stations <- NA
  
  data_season <- data %>%
    filter(!is.na(obs_wd_prev)) %>%
    group_by(station, season, obs_wd_prev) %>%
    summarise(p = mean(obs_wd, na.rm = TRUE)) %>%
    pivot_wider(names_from = obs_wd_prev, values_from = p, names_prefix = "p") %>%
    rename(pd_target = pFALSE, pw_target = pTRUE)
  
  results <- tibble(
    station = character(),
    p0 = numeric(),
    p0_w = numeric(),
    p0_d = numeric(),
    fit_all = list(),
    fit_markov_obs = list(),
    fit_markov_est = list(),
    delta_params = list()
  )
  
  for (s in stations) {
    data_s <- dplyr::filter(data, station == s)
    data_season_s <- dplyr::filter(data_season, station == s)
    data_s <- left_join(data_s, data_season_s, by = c("station", "season"))
    p0 <- data_s %>%
      filter(est == 0) %>%
      summarise(p0 = mean(obs_wd, na.rm = TRUE)) %>%
      pull(p0)
    p0_w <- data_s %>%
      filter(est == 0 & obs_wd_prev) %>%
      summarise(p0_w = mean(obs_wd, na.rm = TRUE)) %>%
      pull(p0_w)
    p0_d <- data_s %>%
      filter(est == 0 & !obs_wd_prev) %>%
      summarise(p0_d = mean(obs_wd, na.rm = TRUE)) %>%
      pull(p0_d)
    
    data_s_non_zero <- dplyr::filter(data_s, est > 0)
    fit_all <- glm(obs_wd ~ est + season, data = data_s_non_zero, 
                   family = binomial, na.action = na.exclude)
    fit_markov_obs <- glm(obs_wd ~ est + season + obs_wd_prev, 
                          data = data_s_non_zero, family = binomial,
                          na.action = na.exclude)
    
    # Markov adjustment to model fit_all
    data_s_non_zero$p_base <- predict(fit_all, type = "response")
    
    data_s_non_zero_season <- data_s_non_zero %>% 
      group_by(season) %>%
      summarise(delta_w = solve_delta(p_base[est_wd_prev], first(pw_target)),
                delta_d = solve_delta(p_base[!est_wd_prev], first(pd_target))
                )
    results <- results %>% 
      tibble::add_row(station = s, p0 = p0, p0_w = p0_w, p0_d = p0_d, 
                      fit_all = list(fit_all), fit_markov_obs = list(fit_markov_obs),
                      delta_params = list(data_s_non_zero_season))
  }
  results
}

# logit    <- function(p) log(p / (1 - p)) (logit = qlogis)
# logistic <- function(x) 1 / (1 + exp(-x)) (logistic = plogis)

## Solve: mean(logistic(logit(p) + delta)) = target
solve_delta <- function(p_sub, target, tol = 1e-10, maxiter = 50, 
                        bounds = c(-20, 20)) {
  p_sub <- p_sub[!is.na(p_sub)]
  if (length(p_sub) == 0) return(NA)
  
  fn <- function(delta) mean(plogis(qlogis(p_sub) + delta)) - target
  delta <- uniroot(fn, interval = bounds, tol = 1e-8)$root
  
  return(delta)
}