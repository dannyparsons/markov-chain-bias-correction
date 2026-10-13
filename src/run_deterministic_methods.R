library(here)
library(readr)
library(dplyr)
library(fitdistrplus)
library(lubridate)
library(ggplot2)
library(purrr)
library(tidyr)
library(tibble)

# Setup -------------------------------------------------------------------

zimbabwe <- read_rds(here("data", "zimbabwe_with_sre.RDS"))

zimbabwe <- zimbabwe %>%
  filter(date >= as.Date("1979-01-01") & date <= as.Date("2023-06-30"))

zimbabwe_stations <- readr::read_csv(here("data", "zimbabwe_stations.csv"))

source(here("src", "methods.R"))
source(here("src", "helper_funs.R"))

blocks <- as.Date(c("1979-01-01", "1989-01-01", "1999-01-01", 
                    "2009-01-01", "2023-07-01"))

# Testing method correctness ----------------------------------------------

m_thresh <- markov_thresholds(zimbabwe, obs_col = "rain", est_col = "agera5_rain",
                              season_col = "season", station_col = "station",
                              tol = 1e-2, max_it = 20)

m_thresh$season <- factor(m_thresh$season, levels = c("dry", 10:12, 1:4))

zimbabwe_t <- zimbabwe %>% 
  left_join(m_thresh %>% dplyr::select(station, season, t0, t_w, t_d), 
            by = c("station", "season"))

zimbabwe_t <- zimbabwe_t %>%
  mutate(agera5_bc_rainday = conditional_wd(agera5_rain, t_w, t_d, t0),
         agera5_bc_rainday_lag = lag(agera5_bc_rainday))

zimbabwe_t_by_season <- zimbabwe_t %>% 
  group_by(station, season) %>%
  summarise(p = mean(rainday, na.rm = TRUE),
            p_bc = mean(agera5_bc_rainday, na.rm = TRUE),
            p_w = mean(rainday[rainday_lag], na.rm = TRUE),
            p_w_bc = mean(agera5_bc_rainday[agera5_bc_rainday_lag], na.rm = TRUE),
            p_d = mean(rainday[!rainday_lag], na.rm = TRUE),
            p_d_bc = mean(agera5_bc_rainday[!agera5_bc_rainday_lag], na.rm = TRUE))

# Generate BC data --------------------------------------------------------

source(here("src", "methods.R"))

zimbabwe_bc <- markov_loci(zimbabwe, obs_col = "rain", est_col = "agera5_rain", 
                 season_col = "season", station_col = "station",
                 blocks = blocks)

saveRDS(zimbabwe_bc, here("data", "BC_data", "zimbabwe_agera5_bc.RDS"))

# Calibration Results -----------------------------------------------------

m_thresh <- list()
stations <- unique(zimbabwe$station)
cal_sets <- c("B1", "B2", "B3", "B4")
for (st in stations) {
  print(st)
  data_st <- zimbabwe %>% filter(station == st)
  for (b in 1:(length(blocks) - 1)) {
    data_cal <- data_st %>%
      filter(!(date >= blocks[b] & date < blocks[b + 1]))
    # Calculate thresholds on calibration data
    m_thresh_i <- markov_thresholds(data_cal, obs_col = "rain", est_col = "agera5_rain", 
                                    season_col = "season", station_col = "station")
    m_thresh_i$calibration_set <- paste0("{", paste0(cal_sets[-i], collapse = ", "), "}")
    m_thresh[[length(m_thresh) + 1]] <- m_thresh_i
  }
}
m_thresh <- bind_rows(m_thresh)
m_thresh$station <- recode(m_thresh$station,
                           Buffalo_Range = "Buffalo Range",
                           Mt_Darwin = "Mt Darwin")

m_thresh_res <- m_thresh %>%
  dplyr::select(calibration_set, station:iterations, 
                p0_obs:p_d_est,
                s_all, s_wet, s_dry,
                gamma_est_all, gamma_est_wet, gamma_est_dry) %>%
  mutate(gamma_all_shape = map_dbl(gamma_est_all, ~ .x$estimate["shape"]),
         gamma_all_rate = map_dbl(gamma_est_all, ~ .x$estimate["rate"]),
         gamma_wet_shape = map_dbl(gamma_est_wet, ~ if(!is.null(.x)) .x$estimate["shape"] else NA_real_),
         gamma_wet_rate = map_dbl(gamma_est_wet, ~ if(!is.null(.x)) .x$estimate["rate"] else NA_real_),
         gamma_dry_shape = map_dbl(gamma_est_dry, ~ .x$estimate["shape"]),
         gamma_dry_rate = map_dbl(gamma_est_dry, ~ .x$estimate["rate"])) %>%
  dplyr::select(-c(gamma_est_all:gamma_est_dry))

m_thresh_res <- m_thresh_res %>% 
  mutate(t_diff = t_w - t_d) %>%
  relocate(t_diff, .after = t_d)

# Table in supplementary material 1
m_thresh_res %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x))) %>%
  write.csv(here("results", "TableS1.csv"), row.names = FALSE)

m_thresh_p <- m_thresh %>%
  dplyr::select(calibration_set, station:iterations, p0_obs:p_d_est) %>%
  pivot_longer(
    cols = c(p0_obs:p_d_est),
    names_to = c("ptype", "source"),
    names_pattern = "(p0|p_w|p_d)_(obs|est)",
    values_to = "p"
  ) %>%
  pivot_wider(
    names_from  = source,
    values_from = p
  ) %>%
  mutate(
    ptype = factor(ptype, levels = c("p0", "p_w", "p_d")),
    ptype_mc = factor(ifelse(ptype == "p0", "Unconditional", "Conditional"),
                      levels = c("Unconditional", "Conditional")),
    diff = est - obs
  )

ggplot(m_thresh_p, 
       aes(x = obs, y = est, colour = ptype)) +
  geom_abline() +
  geom_point() +
  labs(x = "Target Probability", y = "Calibrated Probabilty", colour = "Probability") +
  scale_colour_manual(
    labels = c(p0 = expression(p[0]), p_w = expression(p[w]), p_d = expression(p[d])),
    values = c(p0 = "#E31A1C", p_d = "dodgerblue2", p_w = "green4")) +
  coord_fixed(ratio = 1) +
  facet_grid(rows = vars(ptype_mc), cols = vars(station), axes = "all_x") +
  base_theme()

ggsave(here("results", "Fig2.png"), bg = "white",
       dpi = 600, width = 12, height = 6)
