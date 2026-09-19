library(here)
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(purrr)
library(verification)
library(naflex)
library(boot)

# Setup -------------------------------------------------------------------

source(here("src", "helper_funs.R"))
source(here("src", "bootstrap_funs.R"))

zimbabwe_bc <- read_rds(here("data", "BC_data", "zimbabwe_agera5_bc.RDS"))

zimbabwe_bc_stack <- zimbabwe_bc %>%
  dplyr::select(station, date, year, month, day, season, rain, agera5_rain, est_loci:est_qm_gamma_mk) %>%
  pivot_longer(cols = c(rain, agera5_rain, est_loci:est_qm_gamma_mk), names_to = "source", values_to = "rr")

blocks <- as.Date(c("1979-01-01", "1989-01-01", "1999-01-01", 
                    "2009-01-01", "2023-07-01"))

zimbabwe_bc_stack <- zimbabwe_bc_stack %>%
  mutate(block = cut(date, breaks = blocks, include.lowest = TRUE, right = FALSE,
                     labels = paste0("block_", 1:4)))

rm(zimbabwe_bc)

# 1 Aug = 214
s_doy_start <- 214

zimbabwe_bc_stack <- zimbabwe_bc_stack %>%
  group_by(station, date) %>%
  # Only include days in analysis when there was a gauge value
  mutate(rr = if (any(is.na(rr))) NA else rr) %>%
  group_by(station, source) %>%
  mutate(rainday = rr > 0.85,
         lag_rainday = lag(rainday),
         month = factor(month(date), levels = c(8:12, 1:7)),
         month_abb = factor(month(date, label = TRUE, abbr = TRUE), levels = c(month.abb[8:12], month.abb[1:7])),
         s_year = ifelse(month %in% 1:7, year - 1, year),
         doy = yday_366(date),
         s_doy = (doy - s_doy_start + 1) %% 366,
         s_doy = ifelse(s_doy == 0, 366, s_doy)
         )

zimbabwe_bc_stack$source <- factor(zimbabwe_bc_stack$source,
  levels = c("rain", "agera5_rain", "est_loci", "est_loci_mk", "est_qm_gamma",  
             "est_qm_gamma_mk", "est_qm_empirical", "est_qm_empirical_mk"),
  labels = c("Gauge", "AgERA5", "LOCI", "MC LOCI", "QM", "MC QM",
             "QM-Empirical", "MC QM-Empirical"))

zimbabwe_bc_stack$station <- recode(zimbabwe_bc_stack$station, 
                                    Buffalo_Range = "Buffalo Range",
                                    Mt_Darwin = "Mt Darwin")

# Only need to evaluate one existing BC method and one modified method
# since rain days are constructed the same way
occurrence_source <- c("Gauge", "AgERA5", "LOCI", "MC LOCI")

amounts_source <- c(occurrence_source, "QM", "MC QM")

# RAINFALL OCCURRENCE -----------------------------------------------------

zimbabwe_bc_stack_occ <- zimbabwe_bc_stack %>%
  filter(source %in% occurrence_source)

zimbabwe_bc_stack_occ$source <- recode(zimbabwe_bc_stack_occ$source,
                                       LOCI = "LOCI/QM", `MC LOCI` = "MC")

col_values_occ <- c(
  Gauge = "black",
  AgERA5 = "#E31A1C",
  `LOCI/QM` = "dodgerblue2",
  MC = "green4"
)

col_scale_occ <- scale_colour_manual(
  values =  col_values_occ
)

col_fill_occ <- scale_fill_manual(
  values = col_values_occ
)

# Monthly climatology -----------------------------------------------------

zim_monthly_occ <- zimbabwe_bc_stack_occ %>%
  group_by(station, source, month_abb, year) %>%
  summarise(n_rain = sum(rainday %>% na_omit_if(n = 10, consec = 4))) %>%
  summarise(n_rain = mean(n_rain, na.rm = TRUE))

ggplot(zim_monthly_occ,
       aes(x = month_abb, y = n_rain, colour = source, group = source,
           linewidth = source)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(station)) +
  labs(colour = "Source",
       linewidth = "Source",
       x = "Month",
       y = "Mean number of rain days") +
  col_scale_occ +
  scale_linewidth_manual(values = c("Gauge" = 1.2,
                                   "AgERA5" = 0.5,
                                   "LOCI/QM" = 0.5,
                                   "MC" = 0.5)) +
  base_theme(panel.grid.minor = TRUE)

ggsave(here("results", "Fig3.png"), bg = "white",
       dpi = 600, width = 12, height = 6)
  
# Annual summaries --------------------------------------------------------

zim_annual_occ <- zimbabwe_bc_stack_occ %>%
  group_by(station, source, s_year) %>%
  # first remove incomplete years
  filter(!(s_year %in% c(1978, 2022))) %>%
  summarise(n_rain = sum(rainday %>% na_omit_if(n = 27, consec = 20))) %>%
  ungroup()

zim_annual_dryspells <- zimbabwe_bc_stack_occ %>%
  group_by(station, source, s_year) %>%
  # first remove incomplete years (2022 retained since complete up to March 2023)
  filter(!(s_year %in% c(1978))) %>%
  filter(month %in% c(10:12, 1:3)) %>%
  summarise(
    max_dry_spell = {
      r <- rle(!rainday)
      max(r$lengths[r$values], na.rm = TRUE)
      },
    na = sum(is.na(rainday)),
    naconsec = na_consec(rainday)) %>%
  mutate(max_dry_spell = if_else(na > 27 | naconsec > 20, NA, max_dry_spell)) %>%
  dplyr::select(-c(na, naconsec)) %>%
  ungroup()

zim_annual_occ <- full_join(zim_annual_occ, zim_annual_dryspells, 
                            by = c("station", "source", "s_year"))

zim_annual_occ_station <- zim_annual_occ %>% 
  filter(source == "Gauge") %>% 
  rename(n_rain_station = n_rain,
         max_dry_spell_station = max_dry_spell) %>% 
  dplyr::select(-source)

zim_annual_occ_wide <- zim_annual_occ %>% 
  filter(source != "Gauge") %>%
  left_join(zim_annual_occ_station, by = c("station", "s_year"))

zim_annual_occ_wide <- zim_annual_occ_wide %>% 
  group_by(station, source) %>%
  mutate(n_rain_diff = n_rain - n_rain_station,
         max_dry_spell_diff = max_dry_spell - max_dry_spell_station)

zim_annual_occ_metrics <- zim_annual_occ_wide %>% 
  group_by(station, source) %>%
  summarise(n_rain_me = mean(n_rain_diff, na.rm = TRUE),
            max_dry_spell_me = mean(max_dry_spell_diff, na.rm = TRUE),
            n_rain_cor = cor(n_rain, n_rain_station, use = "complete.obs"),
            max_dry_spell_cor = cor(max_dry_spell, max_dry_spell_station, use = "complete.obs"))

# Annual number of rain days
# (2022 is already NA for n_rain but filtering removes it from influencing axis)
ggplot(zim_annual_occ %>% filter(s_year != 2022), 
       aes(x = s_year, y = n_rain, colour = source)) +
  geom_line() +
  geom_point(size = 0.6) +
  facet_wrap(vars(station)) +
  col_scale_occ +
  base_theme() +
  labs(colour = "Source",
       x = "Year",
       y = "Number of rain days")

ggsave(here("results", "Fig4.png"), bg = "white",
       dpi = 600, width = 12, height = 6)

tbl_annual_occ_nrain <- zim_annual_occ_metrics %>%
  dplyr::select(station, source, n_rain_me, n_rain_cor) %>%
  pivot_longer(cols = c(n_rain_me, n_rain_cor), names_to = "metric",
    values_to = "value") %>%
  mutate(metric = recode(metric, n_rain_me = "ME", n_rain_cor  = "cor"),
         value = round(value, 3)) %>%
  pivot_wider(names_from = c(metric, source), values_from = value, 
              names_sep = "_") %>%
  dplyr::select(station, sort(grep("^ME_", names(.), value = TRUE)),
    sort(grep("^cor_",  names(.), value = TRUE)))

tbl_annual_occ_nrain %>%
  mutate(across(starts_with("ME"), ~ sprintf("%.2f", .x)),
         across(starts_with("cor"), ~ sprintf("%.2f", .x))) %>%
  write.csv(here("results", "Table6.csv"), row.names = FALSE)

# Bootstrap CIs for nram ME and corr (Table 6)
zim_annual_occ_nrain_wide <- zim_annual_occ %>%
  # (2022 is already NA for n_rain but filtering removes it from influencing number of observations)
  filter(s_year != 2022) %>%
  dplyr::select(station, s_year, source, n_rain) %>%
  pivot_wider(names_from = source, values_from = n_rain) %>%
  arrange(station, s_year) %>%
  dplyr::select(station, s_year, Gauge, `LOCI/QM`, MC) %>%
  dplyr::rename(LOCI_QM = `LOCI/QM`)

set.seed(6)

nrain_cis <- zim_annual_occ_nrain_wide %>%
  group_by(station) %>%
  nest() %>%
  mutate(
    bootstrap = map(data, nrain_bootstrap_station, 
                    statistic = me_corr_occ_stats,
                    R = 10000, l = 3)
  ) %>%
  dplyr::select(station, bootstrap) %>%
  unnest(bootstrap)

nrain_cis %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x))) %>%
  #TODO Update table name once finalised
  write.csv(here("results", "Table6CIs.csv"), row.names = FALSE)

# Annual length of longest dry spell (October to March)
ggplot(zim_annual_occ, 
       aes(x = s_year, y = max_dry_spell, colour = source)) +
  geom_line() +
  facet_wrap(vars(station)) +
  col_scale_occ +
  base_theme() +
  labs(colour = "Source",
       x = "Year",
       y = "Maximum dry spell length (days)")

ggsave(here("results", "Fig5.png"), bg = "white",
       dpi = 600, width = 12, height = 6)

tbl_annual_occ_maxdry <- zim_annual_occ_metrics %>%
  dplyr::select(station, source, max_dry_spell_me, max_dry_spell_cor) %>%
  pivot_longer(cols = c(max_dry_spell_me, max_dry_spell_cor), names_to = "metric",
               values_to = "value") %>%
  mutate(metric = recode(metric, max_dry_spell_me = "ME", max_dry_spell_cor  = "cor"),
         value = round(value, 3)) %>%
  pivot_wider(names_from = c(metric, source), values_from = value, 
              names_sep = "_") %>%
  dplyr::select(station, sort(grep("^ME_", names(.), value = TRUE)),
                sort(grep("^cor_",  names(.), value = TRUE)))

tbl_annual_occ_maxdry %>%
  mutate(across(starts_with("ME"), ~ sprintf("%.2f", .x)),
         across(starts_with("cor"), ~ sprintf("%.2f", .x))) %>%
  write.csv(here("results", "Table7.csv"), row.names = FALSE)

# Bootstrap CIs for max dry spell ME and corr (Table 7)
zim_annual_occ_maxdryspell_wide <- zim_annual_occ %>%
  dplyr::select(station, s_year, source, max_dry_spell) %>%
  pivot_wider(names_from = source, values_from = max_dry_spell) %>%
  arrange(station, s_year) %>%
  dplyr::select(station, s_year, Gauge, `LOCI/QM`, MC) %>%
  dplyr::rename(LOCI_QM = `LOCI/QM`)

set.seed(6)

maxdryspell_cis <- zim_annual_occ_maxdryspell_wide %>%
  group_by(station) %>%
  nest() %>%
  mutate(
    bootstrap = map(data, nrain_bootstrap_station, 
                    statistic = me_corr_occ_stats,
                    R = 10000, l = 3)
  ) %>%
  dplyr::select(station, bootstrap) %>%
  unnest(bootstrap)

maxdryspell_cis %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x))) %>%
  #TODO Update table name once finalised
  write.csv(here("results", "Table7CIs.csv"), row.names = FALSE)

# Distribution of wet/dry spells ------------------------------------------

dry_spells <- zimbabwe_bc_stack_occ %>%
  # remove incomplete year
  filter(s_year != 1978) %>%
  group_by(station, source, s_year) %>%
  filter(month %in% c(10:12, 1:3)) %>%
  reframe(dry_spell_length = {
    r <- rle(!rainday)
    r$lengths[r$values]
  })

wet_spells <- zimbabwe_bc_stack_occ %>%
  # remove incomplete year
  filter(s_year != 1978) %>%
  group_by(station, source, s_year) %>%
  filter(month %in% c(10:12, 1:3)) %>%
  reframe(wet_spell_length = {
    r <- rle(rainday)
    r$lengths[r$values]
  })

ggplot(wet_spells, aes(x = wet_spell_length, colour = source)) +
  stat_ecdf(linewidth = 1) +
  scale_x_log10(breaks = c(1, 5, 10, 30, 50)) +
  facet_wrap(vars(station)) +
  labs(
    x = "Wet spell length (days)",
    y = "Cumulative frequency of wet spell lengths",
    colour = "Source"
  ) +
  col_scale_occ +
  base_theme(panel.grid.minor = FALSE)

ggsave(here("results", "Fig6.png"), bg = "white",
       dpi = 600, width = 12, height = 6)

ks_results_wet <- wet_spells %>%
  group_by(station) %>%
  reframe(
    map_dfr(c("AgERA5", "LOCI/QM", "MC"), function(src) {
      test <- ks.test(
        wet_spell_length[source == "Gauge"],
        wet_spell_length[source == src]
      )
      tibble(
        source = src,
        `K-S test statistic` = unname(test$statistic),
        `p value` = test$p.value
      )
    })
  ) %>%
  ungroup()

ks_results_wet %>%
  mutate(across(starts_with("K-S"), ~ sprintf("%.3f", .x)),
         `p value` = ifelse(`p value` < 0.001, "<0.001", 
                            sprintf("%.3f", `p value`))) %>%
  write.csv(here("results", "Table8.csv"), row.names = FALSE)

# Bootstrap CIs for wet spells K-S (Table 8)
set.seed(6)

wet_ks_cis <- wet_spells %>%
  group_by(station) %>%
  nest() %>%
  mutate(
    bootstrap = map(
      data,
      ks_bootstrap_station,
      R = 10000,
      l = 3,
      col = "wet_spell_length"
    )
  ) %>%
  dplyr::select(station, bootstrap) %>%
  unnest(bootstrap)

wet_ks_cis  %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x))) %>%
  #TODO Update table name once finalised
  write.csv(here("results", "Table8CIs.csv"), row.names = FALSE)


ggplot(dry_spells, aes(x = dry_spell_length, colour = source)) +
  stat_ecdf(linewidth = 1) +
  scale_x_log10(breaks = c(1, 5, 10, 30, 50)) +
  facet_wrap(vars(station)) +
  labs(
    x = "Dry spell length (days)",
    y = "Cumulative frequency of dry spell lengths",
    colour = "Source"
  ) +
  col_scale_occ +
  base_theme(panel.grid.minor = FALSE)

ggsave(here("results", "Fig7.png"), bg = "white",
       dpi = 600, width = 12, height = 6)

ks_results_dry <- dry_spells %>%
  group_by(station) %>%
  reframe(
    map_dfr(c("AgERA5", "LOCI/QM", "MC"), function(src) {
      test <- ks.test(
        dry_spell_length[source == "Gauge"],
        dry_spell_length[source == src]
      )
      tibble(
        source = src,
        `K-S test statistic` = unname(test$statistic),
        `p value` = test$p.value
      )
    })
  ) %>%
  ungroup()

ks_results_dry %>%
  mutate(across(starts_with("K-S"), ~ sprintf("%.3f", .x)),
         `p value` = ifelse(`p value` < 0.001, "<0.001", 
                            sprintf("%.3f", `p value`))) %>%
  write.csv(here("results", "Table9.csv"), row.names = FALSE)

# Bootstrap CIs for dry spells K-S (Table 9)
set.seed(6)

dry_ks_cis <- dry_spells %>%
  group_by(station) %>%
  nest() %>%
  mutate(
    bootstrap = map(
      data,
      ks_bootstrap_station,
      R = 10000,
      l = 3,
      col = "dry_spell_length"
    )
  ) %>%
  dplyr::select(station, bootstrap) %>%
  unnest(bootstrap)

dry_ks_cis  %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x))) %>%
  #TODO Update table name once finalised
  write.csv(here("results", "Table9CIs.csv"), row.names = FALSE)

# By block

dry_spells_block <- zimbabwe_bc_stack_occ %>%
  group_by(station, source, block, s_year) %>%
  filter(month %in% c(10:12, 1:3)) %>%
  reframe(dry_spell_length = {
    r <- rle(!rainday)
    r$lengths[r$values]
  })

wet_spells_block <- zimbabwe_bc_stack_occ %>%
  group_by(station, source, block, s_year) %>%
  filter(month %in% c(10:12, 1:3)) %>%
  reframe(wet_spell_length = {
    r <- rle(rainday)
    r$lengths[r$values]
  })

ks_results_wet_block <- wet_spells_block %>%
  group_by(station, block) %>%
  reframe(
    map_dfr(c("AgERA5", "LOCI/QM", "MC"), function(src) {
      test <- ks.test(
        wet_spell_length[source == "Gauge"],
        wet_spell_length[source == src]
      )
      tibble(
        source = src,
        `K-S test statistic` = unname(test$statistic),
        `p value` = test$p.value
      )
    })
  ) %>%
  ungroup()

ks_summary_wet_block <- ks_results_wet_block %>%
  group_by(source) %>%
  summarise(
    mean = mean(`K-S test statistic`, na.rm = TRUE),
    sd = sd(`K-S test statistic`, na.rm = TRUE)
  ) %>%
  mutate(statistic = "K-S test statistic - wet spells")

ks_results_wet_block_wide <- ks_results_wet_block %>%
  dplyr::select(station, block, source, `K-S test statistic`) %>%
  tidyr::pivot_wider(
    names_from = source,
    values_from = `K-S test statistic`
  )

sum(ks_results_wet_block_wide$MC < ks_results_wet_block_wide$`LOCI/QM`, na.rm = TRUE)
sum(ks_results_wet_block_wide$MC < ks_results_wet_block_wide$AgERA5, na.rm = TRUE)

ks_results_dry_block <- dry_spells_block %>%
  group_by(station, block) %>%
  reframe(
    map_dfr(c("AgERA5", "LOCI/QM", "MC"), function(src) {
      test <- ks.test(
        dry_spell_length[source == "Gauge"],
        dry_spell_length[source == src]
      )
      tibble(
        source = src,
        `K-S test statistic` = unname(test$statistic),
        `p value` = test$p.value
      )
    })
  ) %>%
  ungroup()

ks_summary_dry_block <- ks_results_dry_block %>%
  group_by(source) %>%
  summarise(
    mean = mean(`K-S test statistic`, na.rm = TRUE),
    sd = sd(`K-S test statistic`, na.rm = TRUE)
  ) %>%
  mutate(statistic = "K-S test statistic - dry spells")

ks_results_dry_block_wide <- ks_results_dry_block %>%
  dplyr::select(station, block, source, `K-S test statistic`) %>%
  tidyr::pivot_wider(
    names_from = source,
    values_from = `K-S test statistic`
  )

sum(ks_results_dry_block_wide$MC < ks_results_dry_block_wide$`LOCI/QM`, na.rm = TRUE)
sum(ks_results_dry_block_wide$MC < ks_results_dry_block_wide$AgERA5, na.rm = TRUE)

# Seasonal ----------------------------------------------------------------

# Markov Chain Zero Order Rainday Models

fit_zero_order_markov <- function(data) {
  # Fit logistic regression: P(rainday) ~ seasonal harmonics
  glm(rainday ~ 
        sin(2 * pi * s_doy / 366) + cos(2 * pi * s_doy / 366) +
        sin(4 * pi * s_doy / 366) + cos(4 * pi * s_doy / 366) +
        sin(6 * pi * s_doy / 366) + cos(6 * pi * s_doy / 366),
      data = data,
      family = binomial)
}

mc_models_0 <- zimbabwe_bc_stack_occ %>%
  group_by(source, station) %>%
  group_modify(~ tibble(m = list(fit_zero_order_markov(.x)))) %>%
  ungroup()

doy_df <- tibble(s_doy = 1:366, 
                 s_doy_date = as.Date(1:366, origin = as.Date("1999/07/31")))
fitted_list <- list()
for (i in seq_len(nrow(mc_models_0))) {
  
  # Extract info for this model
  src <- mc_models_0$source[i]
  stn <- mc_models_0$station[i]
  mod <- mc_models_0$m[[i]]
  
  # Predict probabilities
  preds <- predict(mod, newdata = doy_df, type = "response")
  
  # Combine into a tibble
  fitted_list[[i]] <- tibble(
    source = src,
    station = stn,
    s_doy = doy_df$s_doy,
    s_doy_date = doy_df$s_doy_date,
    fitted = preds
  )
}
fitted_doy_df_0 <- bind_rows(fitted_list)

ggplot(fitted_doy_df_0, aes(x = s_doy_date, y = fitted, color = source)) +
  geom_line(size = 1) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  facet_wrap(vars(station), axes = "all_x") +
  labs(
    x = "Date",
    y = "Rain day probability",
    color = "Source",
    linewidth = "Source"
  ) +
  col_scale_occ +
  base_theme(panel.grid.minor = FALSE)

ggsave(here("results", "Fig8.png"), bg = "white",
       dpi = 600, width = 12, height = 6)

rain_ref <- fitted_doy_df_0 %>%
  filter(source == "Gauge") %>%
  dplyr::select(station, s_doy, fitted_rain = fitted)

rmse_rainday_0 <- fitted_doy_df_0 %>%
  filter(source != "Gauge") %>%
  left_join(rain_ref, by = c("station", "s_doy")) %>%
  group_by(station, source) %>%
  summarise(RMSE = sqrt(mean((fitted - fitted_rain)^2, na.rm = TRUE))) %>%
  pivot_wider(names_from  = source, values_from = RMSE)

rmse_rainday_0 %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x))) %>%
  write.csv(here("results", "Table9.csv"), row.names = FALSE)

# Markov Chain First Order Rainday Models

fit_first_order_markov <- function(data) {
  # Fit logistic regression: P(rainday) ~ rainday lag + seasonal harmonics
  glm(rainday ~ 
        lag_rainday +
        sin(2 * pi * s_doy / 366) + cos(2 * pi * s_doy / 366) +
        sin(4 * pi * s_doy / 366) + cos(4 * pi * s_doy / 366) +
        sin(6 * pi * s_doy / 366) + cos(6 * pi * s_doy / 366),
      data = data,
      family = binomial)
}

mc_models_1 <- zimbabwe_bc_stack_occ %>%
  group_by(source, station) %>%
  group_modify(~ tibble(m = list(fit_first_order_markov(.x)))) %>%
  ungroup()

doy_df <- expand.grid(lag_rainday = c(TRUE, FALSE), s_doy = 1:366)
doy_df$s_doy_date <- as.Date(doy_df$s_doy, origin = as.Date("1999/07/31"))
fitted_list <- list()
for (i in seq_len(nrow(mc_models_1))) {
  fitted_data <- doy_df
  preds <- predict(mc_models_1$m[[i]], newdata = fitted_data, type = "response")
  fitted_data$fitted <- preds
  fitted_data$source <- mc_models_1$source[i]
  fitted_data$station <- mc_models_1$station[i]
  
  fitted_list[[i]] <- fitted_data
}

fitted_doy_df_1 <- bind_rows(fitted_list)
fitted_doy_df_1$lag_rainday_fct <- 
  factor(ifelse(fitted_doy_df_1$lag_rainday, "Rain", "No Rain"),
         levels = c("Rain", "No Rain"))

ggplot(fitted_doy_df_1, aes(x = s_doy_date, y = fitted, color = source)) +
  geom_line(size = 0.8) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  labs(
    x = "Date",
    y = "Rain day probability",
    color = "Source",
    linetype = "Previous Day State"
  ) +
  col_scale_occ +
  base_theme(panel.grid.minor = FALSE) +
  facet_grid(rows = vars(lag_rainday_fct), cols = vars(station),
             axes = "all_x")

ggsave(here("results", "Fig9.png"), bg = "white",
       dpi = 600, width = 12, height = 6)

rain_ref <- fitted_doy_df_1 %>%
  filter(source == "Gauge") %>%
  dplyr::select(station, s_doy, lag_rainday_fct, fitted_rain = fitted)

rmse_rainday_1 <- fitted_doy_df_1 %>%
  filter(source != "Gauge") %>%
  left_join(rain_ref, by = c("station", "s_doy", "lag_rainday_fct")) %>%
  group_by(station, source, lag_rainday_fct) %>%
  summarise(RMSE = sqrt(mean((fitted - fitted_rain)^2, na.rm = TRUE)))

rmse_rainday_1_wide <- rmse_rainday_1 %>%
  pivot_wider(names_from = source, values_from = RMSE) %>%
  arrange(lag_rainday_fct, station) %>%
  rename(`Previous state` = lag_rainday_fct)

# Include in supplementary material
rmse_rainday_1_wide %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x))) %>%
  write.csv(here("results", "TableS2.csv"), row.names = FALSE)

ggplot(rmse_rainday_1,
       aes(x = lag_rainday_fct, y = RMSE, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(vars(station)) +
  labs(
    x = "Previous Day State",
    y = "RMSE",
    fill = "Source"
    ) +
  col_fill_occ +
  base_theme()

ggsave(here("results", "Fig10.png"), bg = "white",
       dpi = 600, width = 12, height = 6)

# First order RMSE CIs

source(here("src", "bootstrap_funs.R"))
set.seed(6)

mc_first_occ_rmse_cis <- zimbabwe_bc_stack_occ %>%
  group_by(station) %>%
  nest() %>%
  mutate(
    bootstrap = map(
      data,
      mc_first_rmse_bootstrap_station,
      R = 10000,
      parallel = "snow",
      ncpus = 8
    )
  ) %>%
  dplyr::select(station, bootstrap) %>%
  unnest(bootstrap)

mc_first_occ_rmse_cis %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x))) %>%
  write.csv(here("results", "TableS3CIs.csv"), row.names = FALSE)

# By block

# zero order

mc_models_0_block <- zimbabwe_bc_stack_occ %>%
  group_by(source, station, block) %>%
  group_modify(~ tibble(m = list(fit_zero_order_markov(.x)))) %>%
  ungroup()

doy_df <- tibble(s_doy = 1:366, 
                 s_doy_date = as.Date(1:366, origin = as.Date("1999/07/31")))
fitted_list <- list()
for (i in seq_len(nrow(mc_models_0_block))) {
  
  # Extract info for this model
  src <- mc_models_0_block$source[i]
  stn <- mc_models_0_block$station[i]
  blk <- mc_models_0_block$block[i]
  mod <- mc_models_0_block$m[[i]]
  
  # Predict probabilities
  preds <- predict(mod, newdata = doy_df, type = "response")
  
  # Combine into a tibble
  fitted_list[[i]] <- tibble(
    source = src,
    station = stn,
    block = blk,
    s_doy = doy_df$s_doy,
    s_doy_date = doy_df$s_doy_date,
    fitted = preds
  )
}
fitted_doy_df_0_block <- bind_rows(fitted_list)

rain_ref_block <- fitted_doy_df_0_block %>%
  filter(source == "Gauge") %>%
  dplyr::select(station, block, s_doy, fitted_rain = fitted)

rmse_rainday_0_block <- fitted_doy_df_0_block %>%
  filter(source != "Gauge") %>%
  left_join(rain_ref_block, by = c("station", "block", "s_doy")) %>%
  group_by(station, block, source) %>%
  summarise(RMSE = sqrt(mean((fitted - fitted_rain)^2, na.rm = TRUE))) %>%
  group_by(source) %>%
  summarise(mean = mean(RMSE, na.rm = TRUE),
            sd = sd(RMSE, na.rm = TRUE)) 

rmse_rainday_0_block <- rmse_rainday_0_block %>%
  mutate(statistic = "RMSE_0")

# First order

mc_models_1_block <- zimbabwe_bc_stack_occ %>%
  group_by(source, station, block) %>%
  group_modify(~ tibble(m = list(fit_first_order_markov(.x)))) %>%
  ungroup()

doy_df <- expand.grid(lag_rainday = c(TRUE, FALSE), s_doy = 1:366)
doy_df$s_doy_date <- as.Date(doy_df$s_doy, origin = as.Date("1999/07/31"))
fitted_list <- list()
for (i in seq_len(nrow(mc_models_1_block))) {
  fitted_data <- doy_df
  preds <- predict(mc_models_1_block$m[[i]], newdata = fitted_data, type = "response")
  fitted_data$fitted <- preds
  fitted_data$source <- mc_models_1_block$source[i]
  fitted_data$station <- mc_models_1_block$station[i]
  fitted_data$block <- mc_models_1_block$block[i]
  
  fitted_list[[i]] <- fitted_data
}

fitted_doy_df_1_block <- bind_rows(fitted_list)
fitted_doy_df_1_block$lag_rainday_fct <- 
  factor(ifelse(fitted_doy_df_1_block$lag_rainday, "Rain", "No Rain"),
         levels = c("Rain", "No Rain"))

rain_ref_block <- fitted_doy_df_1_block %>%
  filter(source == "Gauge") %>%
  dplyr::select(station, block, s_doy, lag_rainday_fct, fitted_rain = fitted)

rmse_rainday_1_block <- fitted_doy_df_1_block %>%
  filter(source != "Gauge") %>%
  left_join(rain_ref_block, by = c("station", "block", "s_doy", "lag_rainday_fct")) %>%
  group_by(station, block, source, lag_rainday_fct) %>%
  summarise(RMSE = sqrt(mean((fitted - fitted_rain)^2, na.rm = TRUE))) %>%
  group_by(source, lag_rainday_fct) %>%
  summarise(mean = mean(RMSE, na.rm = TRUE),
            sd = sd(RMSE, na.rm = TRUE))

rmse_rainday_1_block <- rmse_rainday_1_block %>%
  mutate(statistic = ifelse(lag_rainday_fct == "Rain", "RMSE_1(W)", "RMSE_1(D)")) %>%
  dplyr::select(-lag_rainday_fct) %>%
  arrange(statistic, source)

# Rainfall occurrence detection -------------------------------------------

zimbabwe_bc_stack_station_occ <- zimbabwe_bc_stack_occ %>%
  ungroup() %>%
  filter(source == "Gauge") %>%
  rename(rr_station = rr,
         rainday_station = rainday) %>%
  dplyr::select(station, date, rainday_station, rr_station)

zimbabwe_bc_comp_occ <- zimbabwe_bc_stack_occ %>%
  filter(source != "Gauge" & source %in% c("AgERA5", "LOCI/QM", "MC")) %>%
  left_join(zimbabwe_bc_stack_station_occ, by = c("station", "date"))

zimbabwe_pod_hss_occ <- zimbabwe_bc_comp_occ %>%
  group_by(station, source) %>%
  summarise(ver = list(verify(rainday_station, rainday, frcst.type = "binary",
                              obs.type = "binary")),
            POD = map_dbl(ver, ~ .x$POD),
            FAR = map_dbl(ver, ~ .x$FAR),
            HSS = map_dbl(ver, ~ .x$HSS)) %>%
  dplyr::select(-ver) %>%
  pivot_longer(cols = c(POD, FAR, HSS), names_to = "metric", values_to = "value",
               names_ptypes = factor(levels = c("POD", "FAR", "HSS")))

ggplot(zimbabwe_pod_hss_occ,
       aes(x = metric, y = value, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(vars(station), axes = "all_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Metric",
    y = "Metric value",
    fill = "Source"
  ) +
  col_fill_occ +
  base_theme(panel.grid.minor = FALSE)

ggsave(here("results", "Fig11.png"), bg = "white",
       dpi = 600, width = 12, height = 6)

# Include table in supplementary material
zimbabwe_pod_hss_occ_wide <- zimbabwe_pod_hss_occ %>%
  pivot_wider(names_from = metric, values_from = value)

zimbabwe_pod_hss_occ_wide %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x))) %>%
  write.csv(here("results", "TableS3.csv"), row.names = FALSE)


# Station-block summaries -------------------------------------------------

station_block_summaries <- bind_rows(ks_summary_wet_block,
                                     ks_summary_dry_block,
                                     rmse_rainday_0_block,
                                     rmse_rainday_1_block)

station_block_summaries <- station_block_summaries %>%
  relocate(statistic)

station_block_summaries <- station_block_summaries %>% 
  mutate(value = sprintf("%.3f (%.3f)", mean, sd)) %>%
  dplyr::select(statistic, source, value) %>%
  pivot_wider(names_from = source, values_from = value)

station_block_summaries %>%
  write.csv(here("results", "Table5.csv"), row.names = FALSE)

# RAINFALL AMOUNTS --------------------------------------------------------

zimbabwe_bc_stack_amt <- zimbabwe_bc_stack %>%
  filter(source %in% amounts_source)

col_values_amt <- c(
  Gauge = "black",
  AgERA5 = "#E31A1C",
  LOCI = "dodgerblue2",
  `MC LOCI` = "green4",
  QM = "#FF7F00",
  `MC QM` = "gold1"
)

col_scale_amt <- scale_colour_manual(
  values =  col_values_amt
)

col_fill_amt <- scale_fill_manual(
  values = col_values_amt
)

# Monthly climatology -----------------------------------------------------

zim_monthly_amt <- zimbabwe_bc_stack_amt %>%
  group_by(station, source, month_abb, year) %>%
  summarise(n_rain = sum(rainday %>% na_omit_if(n = 10, consec = 4)),
            t_rain = sum(rr %>% na_omit_if(n = 10, consec = 4)),
            mean_rain = t_rain / n_rain,
            mean_rain = ifelse(is.infinite(mean_rain), NA, mean_rain),
            max_rain = max(rr %>% na_omit_if(n = 10, consec = 4)),
            max_rain = ifelse(is.infinite(max_rain), NA, max_rain)) %>%
  summarise(t_rain = mean(t_rain, na.rm = TRUE),
            mean_rain = mean(mean_rain, na.rm = TRUE),
            max_rain = mean(max_rain, na.rm = TRUE))

ggplot(zim_monthly_amt, 
       aes(x = month_abb, y = t_rain, colour = source, group = source)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(station)) +
  labs(colour = "Source",
       x = "Month",
       y = "Mean monthly rainfall (mm)") +
  col_scale_amt +
  base_theme()

ggsave(here("results", "Fig12.jpeg"),
       width = 12, height = 6)

ggplot(zim_monthly_amt, 
       aes(x = month_abb, y = mean_rain, colour = source, group = source)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(station)) +
  labs(colour = "Source",
       x = "Month",
       y = "Mean rainfall per rain day (mm/rain day)") +
  col_scale_amt +
  base_theme()

ggsave(here("results", "Fig13.jpeg"),
       width = 12, height = 6)

zim_monthly_amt_metrics <- zim_monthly_amt %>%
  dplyr::select(station, source, month_abb, mean_rain) %>%
  pivot_wider(names_from = source, values_from = mean_rain) %>%
  pivot_longer(cols = -c(station, month_abb, Gauge),
               names_to = "source", values_to = "mean_rain_src") %>%
  filter(!is.na(Gauge), !is.na(mean_rain_src)) %>%
  group_by(station, source) %>%
  summarise(RMSE = sqrt(mean((mean_rain_src - Gauge)^2))) %>%
  pivot_wider(names_from = source, values_from = RMSE)

zim_monthly_amt_metrics %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x))) %>%
  write.csv(here("results", "TableS5.csv"), row.names = FALSE)

ggplot(zim_monthly_amt, 
       aes(x = month_abb, y = max_rain, colour = source, group = source)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(station)) +
  labs(colour = "Source",
       x = "Month",
       y = "Maximum daily rainfall (mm)") +
  col_scale_amt +
  base_theme()

ggsave(here("results", "Fig14.jpeg"),
       width = 12, height = 6)

# Annual summaries --------------------------------------------------------

zim_annual_amt <- zimbabwe_bc_stack_amt %>%
  group_by(station, source, s_year) %>%
  # first remove incomplete years
  filter(!(s_year %in% c(1978, 2022))) %>%
  summarise(n_rain = sum(rainday %>% na_omit_if(n = 27, consec = 20)),
            t_rain = sum(rr %>% na_omit_if(n = 27, consec = 20)),
            mean_rain = t_rain / n_rain,
            max_rain = max(rr %>% na_omit_if(n = 27, consec = 20))) %>%
  ungroup()

zim_annual_amt_station <- zim_annual_amt %>% 
  filter(source == "Gauge") %>% 
  rename(t_rain_station = t_rain, max_rain_station = max_rain,
         mean_rain_station = mean_rain) %>% 
  dplyr::select(-source)

zim_annual_amt_wide <- zim_annual_amt %>% 
  filter(source != "Gauge") %>%
  left_join(zim_annual_amt_station, by = c("station", "s_year"))

zim_annual_amt_wide <- zim_annual_amt_wide %>% 
  group_by(station, source) %>%
  mutate(t_rain_diff = t_rain - t_rain_station,
         max_rain_diff = max_rain - max_rain_station,
         mean_rain_diff = mean_rain - mean_rain_station)

zim_annual_amt_metrics <- zim_annual_amt_wide %>% 
  group_by(station, source) %>%
  summarise(t_rain_me = mean(t_rain_diff, na.rm = TRUE),
            t_rain_cor = cor(t_rain, t_rain_station, use = "complete.obs"),
            t_rain_rsd = hydroGOF::rSD(t_rain, t_rain_station),
            mean_rain_me = mean(mean_rain_diff, na.rm = TRUE),
            mean_rain_cor = cor(mean_rain, mean_rain_station, use = "complete.obs"),
            mean_rain_rsd = hydroGOF::rSD(mean_rain, mean_rain_station),
            max_rain_me = mean(max_rain_diff, na.rm = TRUE),
            max_rain_cor = cor(max_rain, max_rain_station, use = "complete.obs"),
            max_rain_rsd = hydroGOF::rSD(max_rain, max_rain_station))

zim_annual_amt_metrics %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x))) %>%
  write.csv(here("results", "TableS6.csv"), row.names = FALSE)

# Annual total rainfall
ggplot(zim_annual_amt, 
       aes(x = s_year, y = t_rain, colour = source)) +
  geom_line() +
  geom_point(size = 0.6) +
  facet_wrap(vars(station)) +
  labs(colour = "Source",
       x = "Year",
       y = "Total rainfall (mm)") +
  col_scale_amt + 
  base_theme()

ggsave(here("results", "FigS1.png"), dpi = 600,
       width = 12, height = 6)

# Annual mean rainfall
ggplot(zim_annual_amt, 
       aes(x = s_year, y = mean_rain, colour = source)) +
  geom_line() +
  geom_point(size = 0.6) +
  facet_wrap(vars(station)) +
  labs(colour = "Source",
       x = "Year",
       y = "Mean rainfall per rain day (mm/rain day)") +
  col_scale_amt + 
  base_theme()

ggsave(here("results", "Fig15.png"), dpi = 600,
       width = 12, height = 6)

# Annual max rainfall
ggplot(zim_annual_amt, 
       aes(x = s_year, y = max_rain, colour = source)) +
  geom_line() +
  geom_point(size = 0.6) +
  facet_wrap(vars(station)) +
  labs(colour = "Source",
       x = "Year",
       y = "Maximum daily rainfall (mm)") +
  col_scale_amt + 
  base_theme()

ggsave(here("results", "FigS2.png"), dpi = 600,
       width = 12, height = 6)

# Seasonal ----------------------------------------------------------------

# Markov Chain Zero Order Rainfall Models

fit_zero_order_markov_amounts <- function(data) {
  data_rain <- data %>% filter(rainday)
  glm(rr ~ 
        sin(2 * pi * s_doy / 366) + cos(2 * pi * s_doy / 366) +
        sin(4 * pi * s_doy / 366) + cos(4 * pi * s_doy / 366),
      data = data_rain,
      family = Gamma(link = "log"))
}

mc_models_0_amounts <- zimbabwe_bc_stack_amt %>%
  group_by(source, station) %>%
  group_modify(~ tibble(m = list(fit_zero_order_markov_amounts(.x)))) %>%
  ungroup()

doy_df <- tibble(s_doy = 1:366,
                 s_doy_date = as.Date(1:366, origin = as.Date("1999/07/31")))
fitted_list <- list()
for (i in seq_len(nrow(mc_models_0_amounts))) {
  
  src <- mc_models_0_amounts$source[i]
  stn <- mc_models_0_amounts$station[i]
  mod <- mc_models_0_amounts$m[[i]]
  
  preds <- predict(mod, newdata = doy_df, type = "response")
  
  fitted_list[[i]] <- tibble(
    source = src,
    station = stn,
    s_doy = doy_df$s_doy,
    s_doy_date = doy_df$s_doy_date,
    fitted = preds
  )
}
fitted_doy_df_0_amounts <- bind_rows(fitted_list)

ggplot(fitted_doy_df_0_amounts, 
       aes(x = s_doy_date, y = fitted, color = source)) +
  geom_line(size = 0.8) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  facet_wrap(vars(station), axes = "all_x") +
  labs(
    x = "Date",
    y = "Mean rainfall per rain day (mm/rain day)",
    color = "Source",
    linetype = "Source"
  ) +
  col_scale_amt + 
  base_theme()

ggsave(here("results", "Fig16.jpeg"),
       width = 12, height = 6)

rain_ref <- fitted_doy_df_0_amounts %>%
  filter(source == "Gauge") %>%
  dplyr::select(station, s_doy, fitted_rain = fitted)

rmse_rain_amounts_0 <- fitted_doy_df_0_amounts %>%
  filter(source != "Gauge") %>%
  left_join(rain_ref, by = c("station", "s_doy")) %>%
  group_by(station, source) %>%
  summarise(RMSE = sqrt(mean((fitted - fitted_rain)^2, na.rm = TRUE))) %>%
  pivot_wider(names_from  = source, values_from = RMSE)

rmse_rain_amounts_0 %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x))) %>%
  write.csv(here("results", "Table12.csv"), row.names = FALSE)

# By block

mc_models_0_amounts_block <- zimbabwe_bc_stack_amt %>%
  group_by(source, station, block) %>%
  group_modify(~ tibble(m = list(fit_zero_order_markov_amounts(.x)))) %>%
  ungroup()

doy_df <- tibble(s_doy = 1:366,
                 s_doy_date = as.Date(1:366, origin = as.Date("1999/07/31")))
fitted_list <- list()
for (i in seq_len(nrow(mc_models_0_amounts_block))) {
  
  src <- mc_models_0_amounts_block$source[i]
  stn <- mc_models_0_amounts_block$station[i]
  blk <- mc_models_0_amounts_block$block[i]
  mod <- mc_models_0_amounts_block$m[[i]]
  
  preds <- predict(mod, newdata = doy_df, type = "response")
  
  fitted_list[[i]] <- tibble(
    source = src,
    station = stn,
    block = blk,
    s_doy = doy_df$s_doy,
    s_doy_date = doy_df$s_doy_date,
    fitted = preds
  )
}
fitted_doy_df_0_amounts_block <- bind_rows(fitted_list)

ggplot(fitted_doy_df_0_amounts_block, 
       aes(x = s_doy_date, y = fitted, color = source)) +
  geom_line(size = 0.8) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  facet_grid(vars(station), vars(block), axes = "all_x") +
  labs(
    x = "Date",
    y = "Mean rainfall per rain day (mm/rain day)",
    color = "Source",
    linetype = "Source"
  ) +
  col_scale_amt + 
  base_theme()

rain_ref <- fitted_doy_df_0_amounts_block %>%
  filter(source == "Gauge") %>%
  dplyr::select(station, block, s_doy, fitted_rain = fitted)

rmse_rain_amounts_0_block <- fitted_doy_df_0_amounts_block %>%
  filter(source != "Gauge") %>%
  left_join(rain_ref, by = c("station", "block", "s_doy")) %>%
  group_by(station, block, source) %>%
  summarise(RMSE = sqrt(mean((fitted - fitted_rain)^2, na.rm = TRUE)))

rmse_rain_amounts_0_block_summary <- rmse_rain_amounts_0_block %>%
  group_by(source) %>%
  summarise(mean = mean(RMSE, na.rm = TRUE),
            sd = sd(RMSE, na.rm = TRUE)) %>%
  mutate(statistic = "RMSE_0")

# Markov Chain First Order Rainfall Models

fit_first_order_markov_amounts <- function(data) {
  data_rain <- data %>% filter(rainday)
  glm(rr ~ 
        lag_rainday +
        sin(2 * pi * s_doy / 366) + cos(2 * pi * s_doy / 366) +
        sin(4 * pi * s_doy / 366) + cos(4 * pi * s_doy / 366),
      data = data_rain,
      family = Gamma(link = "log"))
}

mc_models_1_amounts <- zimbabwe_bc_stack_amt %>%
  group_by(source, station) %>%
  group_modify(~ tibble(m = list(fit_first_order_markov_amounts(.x)))) %>%
  ungroup()

doy_df <- expand.grid(lag_rainday = c(TRUE, FALSE), s_doy = 1:366)
doy_df$s_doy_date <- as.Date(doy_df$s_doy, origin = as.Date("1999/07/31"))
fitted_list <- list()
for (i in seq_len(nrow(mc_models_1_amounts))) {
  fitted_data <- doy_df
  preds <- predict(mc_models_1_amounts$m[[i]], newdata = fitted_data, type = "response")
  fitted_data$fitted <- preds
  fitted_data$source <- mc_models_1_amounts$source[i]
  fitted_data$station <- mc_models_1_amounts$station[i]
  
  fitted_list[[i]] <- fitted_data
}
fitted_doy_df_1_amounts <- bind_rows(fitted_list)
fitted_doy_df_1_amounts$lag_rainday_fct <- 
  factor(ifelse(fitted_doy_df_1_amounts$lag_rainday, "Rain", "No Rain"),
         levels = c("Rain", "No Rain"))

ggplot(fitted_doy_df_1_amounts, 
       aes(x = s_doy_date, y = fitted, color = source)) +
  geom_line(size = 0.8) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  facet_grid(vars(lag_rainday_fct), vars(station), axes = "all_x") +
  labs(
    x = "Date",
    y = "Mean rainfall per rain day (mm/rain day)",
    color = "Source"
  ) +
  base_theme() +
  col_scale_amt

ggsave(here("results", "Fig17.jpeg"),
       width = 12, height = 6)

rain_ref <- fitted_doy_df_1_amounts %>%
  filter(source == "Gauge") %>%
  dplyr::select(station, s_doy, lag_rainday_fct, fitted_rain = fitted)

rmse_rain_amounts_1 <- fitted_doy_df_1_amounts %>%
  filter(source != "Gauge") %>%
  left_join(rain_ref, by = c("station", "s_doy", "lag_rainday_fct")) %>%
  group_by(station, source, lag_rainday_fct) %>%
  summarise(RMSE = sqrt(mean((fitted - fitted_rain)^2, na.rm = TRUE)))

# Include table in supplementary material
rmse_rain_amounts_1 %>%
  pivot_wider(names_from = lag_rainday_fct, values_from = RMSE) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x))) %>%
  write.csv(here("results", "TableS4.csv"), row.names = FALSE)

ggplot(rmse_rain_amounts_1,
       aes(x = lag_rainday_fct, y = RMSE, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(vars(station), axes = "all_x") +
  labs(
    x = "Lagged rainday",
    y = "RMSE",
    fill = "Source",
  ) +
  base_theme() +
  col_fill_amt

ggsave(here("results", "Fig18.jpeg"),
       width = 12, height = 6)

# By block

mc_models_1_amounts_block <- zimbabwe_bc_stack_amt %>%
  group_by(source, station, block) %>%
  group_modify(~ tibble(m = list(fit_first_order_markov_amounts(.x)))) %>%
  ungroup()

doy_df <- expand.grid(lag_rainday = c(TRUE, FALSE), s_doy = 1:366)
doy_df$s_doy_date <- as.Date(doy_df$s_doy, origin = as.Date("1999/07/31"))
fitted_list <- list()
for (i in seq_len(nrow(mc_models_1_amounts_block))) {
  fitted_data <- doy_df
  preds <- predict(mc_models_1_amounts_block$m[[i]], newdata = fitted_data, type = "response")
  fitted_data$fitted <- preds
  fitted_data$source <- mc_models_1_amounts_block$source[i]
  fitted_data$station <- mc_models_1_amounts_block$station[i]
  fitted_data$block <- mc_models_1_amounts_block$block[i]
  
  fitted_list[[i]] <- fitted_data
}
fitted_doy_df_1_amounts_block <- bind_rows(fitted_list)
fitted_doy_df_1_amounts_block$lag_rainday_fct <- 
  factor(ifelse(fitted_doy_df_1_amounts_block$lag_rainday, "Rain", "No Rain"),
         levels = c("Rain", "No Rain"))

ggplot(fitted_doy_df_1_amounts_block, 
       aes(x = s_doy_date, y = fitted, color = source)) +
  geom_line(size = 0.8) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  facet_grid(vars(lag_rainday_fct), vars(station, block), axes = "all_x") +
  labs(
    x = "Date",
    y = "Mean rainfall per rain day (mm/rain day)",
    color = "Source"
  ) +
  base_theme() +
  col_scale_amt

rain_ref <- fitted_doy_df_1_amounts_block %>%
  filter(source == "Gauge") %>%
  dplyr::select(station, block, s_doy, lag_rainday_fct, fitted_rain = fitted)

rmse_rain_amounts_1_block <- fitted_doy_df_1_amounts_block %>%
  filter(source != "Gauge") %>%
  left_join(rain_ref, by = c("station", "block", "s_doy", "lag_rainday_fct")) %>%
  group_by(station, block, source, lag_rainday_fct) %>%
  summarise(RMSE = sqrt(mean((fitted - fitted_rain)^2, na.rm = TRUE)))

rmse_rain_amounts_1_block_summary <- rmse_rain_amounts_1_block %>%
  group_by(source, lag_rainday_fct) %>%
  summarise(mean = mean(RMSE, na.rm = TRUE),
            sd = sd(RMSE, na.rm = TRUE))

rmse_rain_amounts_1_block_summary <- rmse_rain_amounts_1_block_summary %>%
  mutate(statistic = ifelse(lag_rainday_fct == "Rain", "RMSE_1(W)", "RMSE_1(D)")) %>%
  dplyr::select(-lag_rainday_fct) %>%
  arrange(statistic, source)

# Station-block summaries -------------------------------------------------

station_block_summaries_amounts <- bind_rows(rmse_rain_amounts_0_block_summary,
                                             rmse_rain_amounts_1_block_summary)

station_block_summaries_amounts <- station_block_summaries_amounts %>%
  relocate(statistic)

station_block_summaries_amounts <- station_block_summaries_amounts %>% 
  mutate(value = sprintf("%.3f (%.3f)", mean, sd)) %>%
  dplyr::select(statistic, source, value) %>%
  pivot_wider(names_from = source, values_from = value)

station_block_summaries_amounts %>%
  write.csv(here("results", "Table11.csv"), row.names = FALSE)

# POD and HSS for rainfall categories -------------------------------------

cat_labs <- c("No Rain", "Light Rain", 
              "Moderate Rain", "Heavy Rain", 
              "Violent Rain")

zimbabwe_bc_stack_amt <- zimbabwe_bc_stack_amt %>%
  mutate(rain_cat = cut(rr, c(0, 0.85, 5, 20, 40, Inf), include.lowest = TRUE,
                         right = FALSE, labels = cat_labs))

zimbabwe_bc_stack_amt_station <- zimbabwe_bc_stack_amt %>% 
  ungroup() %>%
  filter(source == "Gauge") %>% 
  rename(rain_cat_station = rain_cat, rr_station = rr) %>% 
  dplyr::select(station, date, rain_cat_station, rr_station)

zimbabwe_bc_stack_amt_wide <- zimbabwe_bc_stack_amt %>% 
  filter(source != "Gauge") %>%
  left_join(zimbabwe_bc_stack_amt_station, by = c("station", "date"))

zimbabwe_pod_hss_amt <- zimbabwe_bc_stack_amt_wide %>%
  group_by(station, source) %>%
  filter(!is.na(rain_cat) & !is.na(rain_cat_station)) %>%
  filter(month %in% c(10:12, 1:3)) %>%
  summarise(
    pod_no = sum(rain_cat == "No Rain" & rain_cat_station == "No Rain", na.rm = TRUE) / sum(rain_cat_station == "No Rain", na.rm = TRUE),
    pod_light = sum(rain_cat == "Light Rain" & rain_cat_station == "Light Rain", na.rm = TRUE) / sum(rain_cat_station == "Light Rain", na.rm = TRUE),
    pod_moderate = sum(rain_cat == "Moderate Rain" & rain_cat_station == "Moderate Rain", na.rm = TRUE) / sum(rain_cat_station == "Moderate Rain", na.rm = TRUE),
    pod_heavy = sum(rain_cat == "Heavy Rain" & rain_cat_station == "Heavy Rain", na.rm = TRUE) / sum(rain_cat_station == "Heavy Rain", na.rm = TRUE),
    pod_violent = sum(rain_cat == "Violent Rain" & rain_cat_station == "Violent Rain", na.rm = TRUE) / sum(rain_cat_station == "Violent Rain", na.rm = TRUE),
    ver = list(verify(rain_cat_station, rain_cat, frcst.type = "cat",
                      obs.type = "cat")),
    hss = map_dbl(ver, ~ .x$hss))

zimbabwe_pod_hss_amt_format <- zimbabwe_pod_hss_amt %>%
  dplyr::select(-ver) %>%
  pivot_wider(names_from = source, values_from = hss)

zimbabwe_pod_hss_amt_long <- zimbabwe_pod_hss_amt %>%
  dplyr::select(station, source, starts_with("pod_"), hss) %>%
  pivot_longer(
    cols = c(starts_with("pod_"), "hss"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = gsub("pod_", "", metric),
    metric = factor(metric,
                      levels = c("no", "light", "moderate", "heavy", "violent", "hss"),
                      labels = c("No Rain", "Light Rain", "Moderate Rain", "Heavy Rain", "Violent Rain", "HSS"))
  )

ggplot(zimbabwe_pod_hss_amt_long, aes(x = metric, y = value, fill = source)) +
  geom_col(position = position_dodge()) +
  facet_wrap(vars(station)) +
  labs(
    x = "Rainfall Category / Metric",
    y = "Metric Value",
    fill = "Source",
  ) +
  col_fill_amt +
  base_theme() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
  )

ggsave(here("results", "Fig19.jpeg"),
       width = 12, height = 6)
