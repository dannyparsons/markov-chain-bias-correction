me_corr_occ_stats <- function(x) {
  Gauge <- x[, "Gauge"]
  LOCI_QM  <- x[, "LOCI_QM"]
  MC    <- x[, "MC"]
  c(me_diff = mean(MC - Gauge, na.rm = TRUE) - mean(LOCI_QM - Gauge, na.rm = TRUE),
    corr_diff = cor(MC, Gauge, use = "complete.obs") - cor(LOCI_QM, Gauge, use = "complete.obs"))
}

nrain_bootstrap_station <- function(dat, statistic, R = 10000, l = 3) {
  
  x <- dat %>%
    dplyr::select(Gauge, LOCI_QM, MC) %>%
    as.matrix()
  
  nrain_boot <- tsboot(
    tseries = x,
    statistic = statistic,
    R = R,
    l = l,
    sim = "fixed"
  )

  observed_diff <- me_corr_occ_stats(x)
  
  me_ci <- boot.ci(nrain_boot, conf = 0.95, type = "perc", index = 1)
  corr_ci <- boot.ci(nrain_boot, conf = 0.95, type = "perc", index = 2)
  
  tibble(
    me_diff = observed_diff["me_diff"],
    me_ci = sprintf("(%.2f, %.2f)", me_ci$percent[4], me_ci$percent[5]),
    corr_diff = observed_diff["corr_diff"],
    corr_ci = sprintf("(%.2f, %.2f)", corr_ci$percent[4], corr_ci$percent[5])
  )
}

ks_bootstrap_station <- function(dat, R = 10000, l = 3, col) {
  
  # Nested to have one row per year
  # so that bootstrap sampling is done on the years
  dat_nested <- dat %>%
    arrange(s_year) %>%
    group_by(s_year) %>%
    nest() %>%
    ungroup()
  
  year_source <- dat_nested$data %>%
    map(~ list(
      gauge = .x[[col]][.x$source == "Gauge"],
      loci_qm = .x[[col]][.x$source == "LOCI/QM"],
      mc = .x[[col]][.x$source == "MC"]
    ))
  
  n_years <- nrow(dat_nested)
  
  year_ts <- 1:n_years
  
  # Statistic function receives the resampled row numbers
  statistic <- function(samp_rows) {
    
    sampled <- year_source[samp_rows]
    gauge <- unlist(map(sampled, "gauge"), use.names = FALSE)
    loci_qm <- unlist(map(sampled, "loci_qm"), use.names = FALSE)
    mc <- unlist(map(sampled, "mc"), use.names = FALSE)
    c(ks_diff = unname(ks.test(mc, gauge)$statistic - ks.test(loci_qm, gauge)$statistic))
  }
  
  # Observed statistic
  observed <- statistic(year_ts)
  
  # bootstrap
  ks_boot <- tsboot(
    tseries = year_ts,
    statistic = statistic,
    R = R,
    l = l,
    sim = "fixed"
  )
  
  ks_ci <- boot.ci(
    ks_boot,
    conf = 0.95,
    type = "perc",
    index = 1
  )
  
  tibble(
    ks_diff = observed["ks_diff"],
    ks_ci = sprintf("(%.3f, %.3f)", ks_ci$percent[4], ks_ci$percent[5])
  )
}


# 1st order MC RMSE -------------------------------------------------------

mc_first_rmse_bootstrap_station <- function(dat, R = 10000) {
  print("station")
  
  # Model formula
  mc_first_occ_formula <- rainday ~
    lag_rainday +
    sin(2 * pi * s_doy / 366) + cos(2 * pi * s_doy / 366) +
    sin(4 * pi * s_doy / 366) + cos(4 * pi * s_doy / 366) +
    sin(6 * pi * s_doy / 366) + cos(6 * pi * s_doy / 366)
  
  # Prediction data
  doy_df_w <- expand.grid(
    lag_rainday = TRUE,
    s_doy = 1:366
  )
  
  doy_df_d <- expand.grid(
    lag_rainday = FALSE,
    s_doy = 1:366
  )
  
  annual_data <- dat %>%
    arrange(s_year, source, s_doy) %>%
    group_by(s_year, source) %>%
    nest() %>%
    ungroup()
  
  gauge_annual <- annual_data %>%
    filter(source == "Gauge") %>%
    arrange(s_year)
  
  loci_qm_annual <- annual_data %>%
    filter(source == "LOCI/QM") %>%
    arrange(s_year)
  
  mc_annual <- annual_data %>%
    filter(source == "MC") %>%
    arrange(s_year)
  
  statistic <- function(dat, samp_rows) {
    gauge_data <- dplyr::bind_rows(gauge_annual$data[samp_rows])
    loci_qm_data <- dplyr::bind_rows(loci_qm_annual$data[samp_rows])
    mc_data <- dplyr::bind_rows(mc_annual$data[samp_rows])
    
    fit_gauge <- glm(
      mc_first_occ_formula,
      data = gauge_data,
      family = binomial()
    )
    fit_loci_qm <- glm(
      mc_first_occ_formula,
      data = loci_qm_data,
      family = binomial()
    )
    fit_mc <- glm(
      mc_first_occ_formula,
      data = mc_data,
      family = binomial()
    )
    
    X_pred_w <- model.matrix(
      delete.response(terms(fit_gauge)),
      data = doy_df_w
    )
    gauge_w <- plogis(
      X_pred_w %*% coef(fit_gauge)
    )
    loci_qm_w <- plogis(
      X_pred_w %*% coef(fit_loci_qm)
    )
    mc_w <- plogis(
      X_pred_w %*% coef(fit_mc)
    )
    X_pred_d <- model.matrix(
      delete.response(terms(fit_gauge)),
      data = doy_df_d
    )
    gauge_d <- plogis(
      X_pred_d %*% coef(fit_gauge)
    )
    loci_qm_d <- plogis(
      X_pred_d %*% coef(fit_loci_qm)
    )
    mc_d <- plogis(
      X_pred_d %*% coef(fit_mc)
    )
    
    rmse_loci_qm_w <- sqrt(
      mean((gauge_w - loci_qm_w)^2)
    )
    rmse_loci_qm_d <- sqrt(
      mean((gauge_d - loci_qm_d)^2)
    )
    rmse_mc_w <- sqrt(
      mean((gauge_w - mc_w)^2)
    )
    rmse_mc_d <- sqrt(
      mean((gauge_d - mc_d)^2)
    )
    
    c(rmse_w_diff = rmse_mc_w - rmse_loci_qm_w,
      rmse_d_diff = rmse_mc_d - rmse_loci_qm_d)
  }
  
  ann_rows <- 1:nrow(gauge_annual)
  observed <- statistic(annual_data, ann_rows)
  
  rmse_boot <- boot(
    data = ann_rows,
    statistic = statistic,
    R = R, 
    parallel = "snow",
    ncpus = 8
  )
  
  rmse_w_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "perc",
    index = 1
  )

  rmse_d_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "perc",
    index = 2
  )
  
  tibble(
    rmse_w_diff = unname(observed["rmse_w_diff"]),
    rmse_w_ci = sprintf(
      "(%.3f, %.3f)",
      rmse_w_ci$percent[4],
      rmse_w_ci$percent[5]
    ),
    rmse_d_diff = unname(observed["rmse_d_diff"]),
    rmse_d_ci = sprintf(
      "(%.3f, %.3f)",
      rmse_d_ci$percent[4],
      rmse_d_ci$percent[5]
    ),
    boot = list(rmse_boot)
  )
}

# Rainfall occurrence detection -------------------------------------------

occ_detection_bootstrap_station <- function(dat, R = 10000) {
  print("station")
  annual_data <- dat %>%
    arrange(s_year, source, date) %>%
    group_by(source, s_year) %>%
    nest() %>%
    ungroup()
  
  loci_qm_annual <- annual_data %>%
    filter(source == "LOCI/QM") %>%
    arrange(s_year)
  
  mc_annual <- annual_data %>%
    filter(source == "MC") %>%
    arrange(s_year)
  
  stopifnot(identical(loci_qm_annual$s_year, mc_annual$s_year))
  
  statistic <- function(data, samp_rows) {
    loci_qm_data <- dplyr::bind_rows(loci_qm_annual$data[samp_rows])
    mc_data <- dplyr::bind_rows(mc_annual$data[samp_rows])

    loci_qm_ver <- verification::verify(
      loci_qm_data$rainday_station,
      loci_qm_data$rainday,
      frcst.type = "binary",
      obs.type = "binary"
    )
    mc_ver <- verification::verify(
      mc_data$rainday_station,
      mc_data$rainday,
      frcst.type = "binary",
      obs.type = "binary"
    )

    c(
      pod_diff = mc_ver$POD - loci_qm_ver$POD,
      far_diff = mc_ver$FAR - loci_qm_ver$FAR,
      hss_diff = mc_ver$HSS - loci_qm_ver$HSS
    )
  }
  
  observed <- statistic(loci_qm_annual, seq_len(nrow(loci_qm_annual)))

  occ_boot <- boot::boot(
    data = seq_len(nrow(loci_qm_annual)),
    statistic = statistic,
    R = R, 
    parallel = "snow", 
    ncpus = 8
  )
  
  pod_ci <- boot::boot.ci(
    occ_boot,
    conf = 0.95,
    type = "bca",
    index = 1
  )
  
  far_ci <- boot::boot.ci(
    occ_boot,
    conf = 0.95,
    type = "bca",
    index = 2
  )
  
  hss_ci <- boot::boot.ci(
    occ_boot,
    conf = 0.95,
    type = "bca",
    index = 3
  )
  
  tibble(
    pod_diff = unname(observed["pod_diff"]),
    pod_ci = sprintf(
      "(%.3f, %.3f)",
      pod_ci$bca[4],
      pod_ci$bca[5]
    ),
    far_diff = unname(observed["far_diff"]),
    far_ci = sprintf(
      "(%.3f, %.3f)",
      far_ci$bca[4],
      far_ci$bca[5]
    ),
    hss_diff = unname(observed["hss_diff"]),
    hss_ci = sprintf(
      "(%.3f, %.3f)",
      hss_ci$bca[4],
      hss_ci$bca[5]
    )
  )
}

# RAINFALL AMOUNTS --------------------------------------------------------

# RMSE Month Mean Rainfall ------------------------------------------------

monthly_rmse_bootstrap_station <- function(dat, R = 10000, l = 3) {
  
  monthly_dat <- dat %>%
    dplyr::select(source, month_abb, year, mean_rain) %>%
    pivot_wider(
      names_from = source,
      values_from = mean_rain
    ) %>%
    arrange(year, month_abb)
  
  annual_data <- monthly_dat %>%
    group_by(year) %>%
    nest() %>%
    ungroup()
  
  statistic <- function(samp_rows) {
    
    sampled_data <- dplyr::bind_rows(annual_data$data[samp_rows])
    
    sampled_means <- sampled_data %>%
      group_by(month_abb) %>%
      summarise(
        Gauge = mean(Gauge, na.rm = TRUE),
        LOCI = mean(LOCI, na.rm = TRUE),
        `MC LOCI` = mean(`MC LOCI`, na.rm = TRUE),
        QM = mean(QM, na.rm = TRUE),
        `MC QM` = mean(`MC QM`, na.rm = TRUE),
        .groups = "drop"
      )
    
    gauge <- sampled_means$Gauge
    loci <- sampled_means$LOCI
    mc_loci <- sampled_means$`MC LOCI`
    qm <- sampled_means$QM
    mc_qm <- sampled_means$`MC QM`
    
    rmse_loci <- sqrt(mean((loci - gauge)^2))
    rmse_mc_loci <- sqrt(mean((mc_loci - gauge)^2))
    rmse_qm <- sqrt(mean((qm - gauge)^2))
    rmse_mc_qm <- sqrt(mean((mc_qm - gauge)^2))
    
    c(mc_loci_diff = rmse_mc_loci - rmse_loci,
      mc_qm_diff   = rmse_mc_qm - rmse_qm)
  }
  
  observed <- statistic(seq_len(nrow(annual_data)))
  
  rmse_boot <- tsboot(
    tseries = seq_len(nrow(annual_data)),
    statistic = statistic,
    R = R,
    l = l,
    sim = "fixed"
  )
  
  mc_loci_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "perc",
    index = 1
  )
  mc_qm_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "perc",
    index = 2
  )
  
  tibble(
    mc_loci_diff = unname(observed["mc_loci_diff"]),
    mc_loci_ci = sprintf(
      "(%.2f, %.2f)",
      mc_loci_ci$percent[4],
      mc_loci_ci$percent[5]
    ),
    mc_qm_diff = unname(observed["mc_qm_diff"]),
    mc_qm_ci = sprintf(
      "(%.2f, %.2f)",
      mc_qm_ci$percent[4],
      mc_qm_ci$percent[5]
    )
  )
}


# Annual amounts summaries ------------------------------------------------

me_corr_rsd_amt_stats <- function(x) {
  gauge <- x[, "Gauge"]
  loci <- x[, "LOCI"]
  mc_loci <- x[, "MC LOCI"]
  qm <- x[, "QM"]
  mc_qm <- x[, "MC QM"]
  
  c(mc_loci_me_diff = mean(mc_loci - gauge, na.rm = TRUE) - mean(loci - gauge, na.rm = TRUE),
    mc_loci_corr_diff = cor(mc_loci, gauge, use = "complete.obs") - cor(loci, gauge, use = "complete.obs"),
    mc_loci_rsd_diff = hydroGOF::rSD(mc_loci, gauge) - hydroGOF::rSD(loci, gauge),
    mc_qm_me_diff = mean(mc_qm - gauge, na.rm = TRUE) - mean(qm - gauge, na.rm = TRUE),
    mc_qm_corr_diff = cor(mc_qm, gauge, use = "complete.obs") - cor(qm, gauge, use = "complete.obs"),
    mc_qm_rsd_diff = hydroGOF::rSD(mc_qm, gauge) - hydroGOF::rSD(qm, gauge))
}

annual_amt_bootstrap_station <- function(dat, statistic, R = 10000, l = 3) {
  
  x <- dat %>%
    dplyr::select(Gauge, LOCI, `MC LOCI`, QM, `MC QM`) %>%
    as.matrix()
  
  annual_amt_boot <- tsboot(
    tseries = x,
    statistic = statistic,
    R = R,
    l = l,
    sim = "fixed"
  )
  
  observed_diff <- me_corr_rsd_amt_stats(x)
  
  mc_loci_me_ci <- boot.ci(annual_amt_boot, conf = 0.95, type = "perc", index = 1)
  mc_loci_corr_ci <- boot.ci(annual_amt_boot, conf = 0.95, type = "perc", index = 2)
  mc_loci_rsd_ci <- boot.ci(annual_amt_boot, conf = 0.95, type = "perc", index = 3)
  mc_qm_me_ci <- boot.ci(annual_amt_boot, conf = 0.95, type = "perc", index = 4)
  mc_qm_corr_ci <- boot.ci(annual_amt_boot, conf = 0.95, type = "perc", index = 5)
  mc_qm_rsd_ci <- boot.ci(annual_amt_boot, conf = 0.95, type = "perc", index = 6)
  
  tibble(
    mc_loci_me_diff = observed_diff["mc_loci_me_diff"],
    mc_loci_me_ci = sprintf("(%.2f, %.2f)", mc_loci_me_ci$percent[4], mc_loci_me_ci$percent[5]),
    mc_loci_corr_diff = observed_diff["mc_loci_corr_diff"],
    mc_loci_corr_ci = sprintf("(%.2f, %.2f)", mc_loci_corr_ci$percent[4], mc_loci_corr_ci$percent[5]),
    mc_loci_rsd_diff = observed_diff["mc_loci_rsd_diff"],
    mc_loci_rsd_ci = sprintf("(%.2f, %.2f)", mc_loci_rsd_ci$percent[4], mc_loci_rsd_ci$percent[5]),
    mc_qm_me_diff = observed_diff["mc_qm_me_diff"],
    mc_qm_me_ci = sprintf("(%.2f, %.2f)", mc_qm_me_ci$percent[4], mc_qm_me_ci$percent[5]),
    mc_qm_corr_diff = observed_diff["mc_qm_corr_diff"],
    mc_qm_corr_ci = sprintf("(%.2f, %.2f)", mc_qm_corr_ci$percent[4], mc_qm_corr_ci$percent[5]),
    mc_qm_rsd_diff = observed_diff["mc_qm_rsd_diff"],
    mc_qm_rsd_ci = sprintf("(%.2f, %.2f)", mc_qm_rsd_ci$percent[4], mc_qm_rsd_ci$percent[5])
  )
}

# Zero order MC Amounts RMSE ----------------------------------------------

mc_zero_amt_rmse_bootstrap_station <- function(dat, R = 10000) {
  print("station")
  
  fm <- rr ~ 
    sin(2 * pi * s_doy / 366) + cos(2 * pi * s_doy / 366) +
    sin(4 * pi * s_doy / 366) + cos(4 * pi * s_doy / 366) 
  
  # Prediction data
  doy_df <- tibble(s_doy = 1:366)
  
  years <- sort(unique(dat$s_year))
  sources <- c("Gauge", "LOCI", "MC LOCI", "QM", "MC QM")
  
  annual_data <- dat %>%
    arrange(s_year, source, s_doy) %>%
    group_by(s_year, source) %>%
    nest() %>%
    ungroup() %>%
    tidyr::complete(s_year = years, source = sources,
                    fill = list(data = list(tibble()))) %>%
    mutate(data = map(data, ~ .x %>% filter(rainday))) %>%
    arrange(s_year, source)
  
  gauge_annual <- annual_data %>%
    filter(source == "Gauge") %>%
    arrange(s_year)
  
  loci_annual <- annual_data %>%
    filter(source == "LOCI") %>%
    arrange(s_year)
  
  mc_loci_annual <- annual_data %>%
    filter(source == "MC LOCI") %>%
    arrange(s_year)

  qm_annual <- annual_data %>%
    filter(source == "QM") %>%
    arrange(s_year)
  
  mc_qm_annual <- annual_data %>%
    filter(source == "MC QM") %>%
    arrange(s_year)
  
  stopifnot(
    identical(gauge_annual$s_year, loci_annual$s_year),
    identical(gauge_annual$s_year, mc_loci_annual$s_year),
    identical(gauge_annual$s_year, qm_annual$s_year),
    identical(gauge_annual$s_year, mc_qm_annual$s_year)
  )
  
  statistic <- function(dat, samp_rows) {
    gauge_data <- dplyr::bind_rows(gauge_annual$data[samp_rows])
    loci_data <- dplyr::bind_rows(loci_annual$data[samp_rows])
    mc_loci_data <- dplyr::bind_rows(mc_loci_annual$data[samp_rows])
    qm_data <- dplyr::bind_rows(qm_annual$data[samp_rows])
    mc_qm_data <- dplyr::bind_rows(mc_qm_annual$data[samp_rows])
    
    fit_gauge <- glm(
      fm,
      data = gauge_data,
      family = Gamma(link = "log")
    )
    fit_loci <- glm(
      fm,
      data = loci_data,
      family = Gamma(link = "log")
    )
    fit_mc_loci <- glm(
      fm,
      data = mc_loci_data,
      family = Gamma(link = "log")
    )
    fit_qm <- glm(
      fm,
      data = qm_data,
      family = Gamma(link = "log")
    )
    fit_mc_qm <- glm(
      fm,
      data = mc_qm_data,
      family = Gamma(link = "log")
    )
    
    X_pred <- model.matrix(
      delete.response(terms(fit_gauge)),
      data = doy_df
    )
    gauge <- exp(X_pred %*% coef(fit_gauge))
    loci <- exp(X_pred %*% coef(fit_loci))
    mc_loci <- exp(X_pred %*% coef(fit_mc_loci))
    qm <- exp(X_pred %*% coef(fit_qm))
    mc_qm <- exp(X_pred %*% coef(fit_mc_qm))
    
    rmse_loci <- sqrt(
      mean((gauge - loci)^2)
    )
    rmse_mc_loci <- sqrt(
      mean((gauge - mc_loci)^2)
    )
    rmse_qm <- sqrt(
      mean((gauge - qm)^2)
    )
    rmse_mc_qm <- sqrt(
      mean((gauge - mc_qm)^2)
    )
    
    c(rmse_mc_loci_diff = rmse_mc_loci - rmse_loci,
      rmse_mc_qm_diff = rmse_mc_qm - rmse_qm)
  }
  
  ann_rows <- 1:nrow(gauge_annual)
  observed <- statistic(annual_data, ann_rows)
  
  rmse_boot <- boot(
    data = ann_rows,
    statistic = statistic,
    R = R, 
    parallel = "snow",
    ncpus = 8
  )
  
  rmse_mc_loci_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "bca",
    index = 1
  )
  
  rmse_mc_qm_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "bca",
    index = 2
  )
  
  tibble(
    rmse_mc_loci_diff = unname(observed["rmse_mc_loci_diff"]),
    rmse_mc_loci_ci = sprintf(
      "(%.2f, %.2f)",
      rmse_mc_loci_ci$bca[4],
      rmse_mc_loci_ci$bca[5]
    ),
    rmse_mc_qm_diff = unname(observed["rmse_mc_qm_diff"]),
    rmse_mc_qm_ci = sprintf(
      "(%.2f, %.2f)",
      rmse_mc_qm_ci$bca[4],
      rmse_mc_qm_ci$bca[5]
    ),
    boot = list(rmse_boot)
  )
}

# First order MC Amounts RMSE ---------------------------------------------

mc_first_amt_rmse_bootstrap_station <- function(dat, R = 10000) {
  print("station")
  
  fm <- rr ~ 
    lag_rainday +
    sin(2 * pi * s_doy / 366) + cos(2 * pi * s_doy / 366) +
    sin(4 * pi * s_doy / 366) + cos(4 * pi * s_doy / 366)
  
  # Prediction data
  doy_df_w <- tibble(s_doy = 1:366, lag_rainday = TRUE)
  doy_df_d <- tibble(s_doy = 1:366, lag_rainday = FALSE)
  
  years <- sort(unique(dat$s_year))
  
  sources <- c("Gauge", "LOCI", "MC LOCI", "QM", "MC QM")
  
  annual_data <- dat %>%
    group_by(s_year, source) %>%
    nest() %>%
    ungroup() %>%
    tidyr::complete(
      s_year = years,
      source = sources,
      fill = list(data = list(tibble()))
    )  
  gauge_annual <- annual_data %>%
    filter(source == "Gauge") %>%
    arrange(s_year)
  
  loci_annual <- annual_data %>%
    filter(source == "LOCI") %>%
    arrange(s_year)
  
  mc_loci_annual <- annual_data %>%
    filter(source == "MC LOCI") %>%
    arrange(s_year)
  
  qm_annual <- annual_data %>%
    filter(source == "QM") %>%
    arrange(s_year)
  
  mc_qm_annual <- annual_data %>%
    filter(source == "MC QM") %>%
    arrange(s_year)
  
  stopifnot(
    identical(gauge_annual$s_year, loci_annual$s_year),
    identical(gauge_annual$s_year, mc_loci_annual$s_year),
    identical(gauge_annual$s_year, qm_annual$s_year),
    identical(gauge_annual$s_year, mc_qm_annual$s_year)
  )
  
  statistic <- function(dat, samp_rows) {
    gauge_data <- dplyr::bind_rows(gauge_annual$data[samp_rows])
    loci_data <- dplyr::bind_rows(loci_annual$data[samp_rows])
    mc_loci_data <- dplyr::bind_rows(mc_loci_annual$data[samp_rows])
    qm_data <- dplyr::bind_rows(qm_annual$data[samp_rows])
    mc_qm_data <- dplyr::bind_rows(mc_qm_annual$data[samp_rows])
    
    fit_gauge <- glm(
      fm,
      data = gauge_data,
      family = Gamma(link = "log")
    )
    fit_loci <- glm(
      fm,
      data = loci_data,
      family = Gamma(link = "log")
    )
    fit_mc_loci <- glm(
      fm,
      data = mc_loci_data,
      family = Gamma(link = "log")
    )
    fit_qm <- glm(
      fm,
      data = qm_data,
      family = Gamma(link = "log")
    )
    fit_mc_qm <- glm(
      fm,
      data = mc_qm_data,
      family = Gamma(link = "log")
    )
    
    X_pred_w <- model.matrix(
      delete.response(terms(fit_gauge)),
      data = doy_df_w
    )
    gauge_w <- exp(X_pred_w %*% coef(fit_gauge))
    loci_w <- exp(X_pred_w %*% coef(fit_loci))
    mc_loci_w <- exp(X_pred_w %*% coef(fit_mc_loci))
    qm_w <- exp(X_pred_w %*% coef(fit_qm))
    mc_qm_w <- exp(X_pred_w %*% coef(fit_mc_qm))
    
    X_pred_d <- model.matrix(
      delete.response(terms(fit_gauge)),
      data = doy_df_d
    )
    gauge_d <- exp(X_pred_d %*% coef(fit_gauge))
    loci_d <- exp(X_pred_d %*% coef(fit_loci))
    mc_loci_d <- exp(X_pred_d %*% coef(fit_mc_loci))
    qm_d <- exp(X_pred_d %*% coef(fit_qm))
    mc_qm_d <- exp(X_pred_d %*% coef(fit_mc_qm))
    
    rmse_loci_d <- sqrt(mean((gauge_d - loci_d)^2))
    rmse_mc_loci_d <- sqrt(mean((gauge_d - mc_loci_d)^2))
    rmse_qm_d <- sqrt(mean((gauge_d - qm_d)^2))
    rmse_mc_qm_d <- sqrt(mean((gauge_d - mc_qm_d)^2))
    
    rmse_loci_w <- sqrt(mean((gauge_w - loci_w)^2))
    rmse_mc_loci_w <- sqrt(mean((gauge_w - mc_loci_w)^2))
    rmse_qm_w <- sqrt(mean((gauge_w - qm_w)^2))
    rmse_mc_qm_w <- sqrt(mean((gauge_w - mc_qm_w)^2))
    
    c(rmse_mc_loci_w_diff = rmse_mc_loci_w - rmse_loci_w,
      rmse_mc_qm_w_diff = rmse_mc_qm_w - rmse_qm_w,
      rmse_mc_loci_d_diff = rmse_mc_loci_d - rmse_loci_d,
      rmse_mc_qm_d_diff = rmse_mc_qm_d - rmse_qm_d)
  }
  
  ann_rows <- 1:nrow(gauge_annual)
  observed <- statistic(annual_data, ann_rows)
  
  rmse_boot <- boot(
    data = ann_rows,
    statistic = statistic,
    R = R, 
    parallel = "snow",
    ncpus = 8
  )
  
  rmse_mc_loci_w_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "bca",
    index = 1
  )
  rmse_mc_qm_w_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "bca",
    index = 2
  )
  rmse_mc_loci_d_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "bca",
    index = 3
  )
  rmse_mc_qm_d_ci <- boot.ci(
    rmse_boot,
    conf = 0.95,
    type = "bca",
    index = 4
  )
  
  tibble(
    rmse_mc_loci_w_diff = unname(observed["rmse_mc_loci_w_diff"]),
    rmse_mc_loci_w_ci = sprintf(
      "(%.2f, %.2f)",
      rmse_mc_loci_w_ci$bca[4],
      rmse_mc_loci_w_ci$bca[5]
    ),
    rmse_mc_qm_w_diff = unname(observed["rmse_mc_qm_w_diff"]),
    rmse_mc_qm_w_ci = sprintf(
      "(%.2f, %.2f)",
      rmse_mc_qm_w_ci$bca[4],
      rmse_mc_qm_w_ci$bca[5]
    ),
    rmse_mc_loci_d_diff = unname(observed["rmse_mc_loci_d_diff"]),
    rmse_mc_loci_d_ci = sprintf(
      "(%.2f, %.2f)",
      rmse_mc_loci_d_ci$bca[4],
      rmse_mc_loci_d_ci$bca[5]
    ),
    rmse_mc_qm_d_diff = unname(observed["rmse_mc_qm_d_diff"]),
    rmse_mc_qm_d_ci = sprintf(
      "(%.2f, %.2f)",
      rmse_mc_qm_d_ci$bca[4],
      rmse_mc_qm_d_ci$bca[5]
    ),
    boot = list(rmse_boot)
  )
}
