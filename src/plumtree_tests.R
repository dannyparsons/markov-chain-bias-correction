plumtree_n_rain <- zim_annual_amt %>%
  filter(station == "Plumtree", source == "Gauge") %>%
  arrange(s_year)

pettitt_n_rain <- BreakPoints::pettit(plumtree_n_rain$n_rain, n_period = 5)
snht_n_rain <- BreakPoints::SNHT(plumtree_n_rain$n_rain, n_period = 5, simulations = 10000)

pettitt_n_rain
snht_n_rain

plumtree_t_rain <- zim_annual_amt %>%
  filter(station == "Plumtree", source == "Gauge") %>%
  arrange(s_year)

pettitt_t_rain <- BreakPoints::pettit(plumtree_t_rain$t_rain, n_period = 5)
snht_t_rain <- BreakPoints::SNHT(plumtree_t_rain$t_rain, n_period = 5, simulations = 10000)

pettitt_t_rain
snht_t_rain
