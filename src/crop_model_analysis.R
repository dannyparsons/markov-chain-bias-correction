yields <- read_rds(here("data", "yields_mtdarwin.rds"))

metrics_by_source <- yields %>%
  filter(source != "rain") %>%
  inner_join(
    yields %>%
      filter(source == "rain") %>%
      dplyr::select(station, crop_season, Yield_rain = Yield),
    by = c("station", "crop_season")
  ) %>%
  group_by(station, source) %>%
  summarise(
    MAE = mean(abs(Yield - Yield_rain), na.rm = TRUE),
    MdAE = median(abs(Yield - Yield_rain), na.rm = TRUE),
    cor = cor(Yield, Yield_rain, use = "complete.obs")
  )

metrics_by_source %>%
  mutate(across(starts_with("M"), ~ sprintf("%.1f", .x)),
         cor = sprintf("%.3f", cor)) %>%
  write.csv(here("results", "TableCM.csv"), row.names = FALSE)

yields_plot <- yields %>%
  dplyr::select(station, crop_season, source, Yield) %>%
  group_by(station, source) %>%
  complete(crop_season = full_seq(crop_season, 1)) %>%
  ungroup() %>%
  mutate(source = recode(source,
                         "rain" = "Gauge",
                         "agera5_rain" = "AgERA5",
                         "est_loci" = "LOCI",
                         "est_loci_mk" = "MC LOCI",
                         "est_qm_gamma" = "QM",
                         "est_qm_gamma_mk" = "MC QM"),
         source = factor(source, levels = c("Gauge", "AgERA5", "LOCI", "MC LOCI", "QM", "MC QM")))

col_values_amt <- c(
  AgERA5 = "#E31A1C",
  LOCI = "dodgerblue2",
  `MC LOCI` = "green4",
  QM = "#FF7F00",
  `MC QM` = "gold1"
)

col_scale_amt <- scale_colour_manual(
  values =  col_values_amt
)

ggplot(yields_plot,
       aes(x = crop_season, y = Yield, colour = source, group = source)) +
  geom_line(na.rm = FALSE) +
  geom_point() +
  labs(x = "Crop season", y = "Simulated yield", colour = "Rainfall source") +
  col_scale_amt +
  base_theme()

yields_wide <- yields %>%
  dplyr::select(station, crop_season, source, Yield) %>%
  pivot_wider(names_from = source, values_from = Yield) %>%
  arrange(station, crop_season)

yields_plot <- yields_wide %>%
  pivot_longer(cols = c(agera5_rain, est_loci, est_loci_mk,
                        est_qm_gamma, est_qm_gamma_mk),
               names_to = "source", values_to = "yield") %>%
  mutate(source = recode(source,
                         "agera5_rain" = "AgERA5",
                         "est_loci" = "LOCI",
                         "est_loci_mk" = "MC LOCI",
                         "est_qm_gamma" = "QM",
                         "est_qm_gamma_mk" = "MC QM"),
         source = factor(source, levels = c("AgERA5", "LOCI", "MC LOCI", "QM", "MC QM")))

ggplot(yields_plot, aes(x = rain, y = yield)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(colour = source)) +
  coord_fixed(xlim = c(2700, 6500), ylim = c(2700, 6500)) +
  facet_wrap(vars(source)) +
  labs(x = "Gauge-based yield (kg/ha)", y = "Esimate-based yield (kg/ha)") +
  col_scale_amt +
  base_theme()

ggsave(here("results", "FigS3.png"), dpi = 600,
       bg = "white", width = 12, height = 9)

# Bootstrap code

station_boot <- function(dat) {
  
  boot_stat <- function(data, indices) {
    d <- data[indices, ]
    
    calc_metrics <- function(est, obs) {
      c(
        MAE = mean(abs(est - obs), na.rm = TRUE),
        MdAE = median(abs(est - obs), na.rm = TRUE),
        cor = cor(est, obs, use = "complete.obs")
      )
    }
    loci <- calc_metrics(d$est_loci, d$rain)
    mc_loci <- calc_metrics(d$est_loci_mk, d$rain)
    qm <- calc_metrics(d$est_qm_gamma, d$rain)
    mc_qm <- calc_metrics(d$est_qm_gamma_mk, d$rain)
    c(mc_loci - loci, mc_qm - qm)
  }
  
  b <- boot::boot(
    data = dat,
    statistic = boot_stat,
    R = 10000
  )
  
  metric <- rep(c("MAE", "MdAE", "cor"), 2)
  comparison <- rep(c("MC LOCI - LOCI", "MC QM - QM"), each = 3)
  
  cis <- lapply(seq_along(b$t0), function(i) {
    boot::boot.ci(b, conf = 0.95, type = "bca", index = i)
  })
  
  results <- tibble(
    comparison = comparison,
    metric = metric,
    estimate = b$t0,
    lower = sapply(cis, function(ci) ci$bca[4]),
    upper = sapply(cis, function(ci) ci$bca[5])
  )
  list(
    results = results,
    boot = b
  )
}

yield_data <- yields %>%
  dplyr::select(station, crop_season, source, Yield) %>%
  pivot_wider(names_from = source, values_from = Yield) %>%
  arrange(crop_season)

set.seed(6)

yield_boot <- yield_data %>%
  group_by(station) %>%
  nest() %>%
  mutate(boot_out = map(data, station_boot))

yield_boot_results <- yield_boot %>%
  mutate(results = map(boot_out, "results")) %>%
  dplyr::select(station, results) %>%
  unnest(results) %>%
  arrange(metric, comparison)

yield_boot_results %>%
  mutate(
    Difference = case_when(
      metric %in% c("MAE", "MdAE") ~ sprintf("%.1f", estimate),
      metric == "cor" ~ sprintf("%.3f", estimate)
    ),
    `95% CI` = case_when(
      metric %in% c("MAE", "MdAE") ~ sprintf("(%.1f, %.1f)", lower, upper),
      metric == "cor" ~ sprintf("(%.3f, %.3f)", lower, upper)
    )
  ) %>%
  dplyr::select(station, comparison, metric, Difference, `95% CI`) %>%
  write.csv(here("results", "TableCMCIs.csv"), row.names = FALSE)
