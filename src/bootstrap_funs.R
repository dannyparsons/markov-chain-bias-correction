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

