yday_366 <- function(date) {
  temp_doy <- lubridate::yday(date)
  temp_leap <- lubridate::leap_year(date)
  temp_doy[(!is.na(temp_doy)) & temp_doy > 59 & (!temp_leap)] <- 1 + temp_doy[(!is.na(temp_doy)) & temp_doy > 59 & (!temp_leap)]
  return(temp_doy)
}

base_theme <- function(panel.grid.minor = TRUE) {
  th <- theme_minimal() +
    theme(strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          legend.position = "bottom")
  if (!panel.grid.minor) th <- th + theme(panel.grid.minor = element_blank())
  th
}

c25 <- c(
  "dodgerblue2", 
  "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  # "black", remove black so can be used as reference
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

closest_point <- function(points, target) {
  if(is.vector(target)) {
    ind <- which.min(sp::spDistsN1(pts = as.matrix(points), pt = target))
    return(c(points[ind, 1], points[ind, 2]))
  }
  else if(is.data.frame(target)) {
    inds <- apply(target, 1, function(r) which.min(sp::spDistsN1(pts = as.matrix(points),
                                                                 pt = c(r[[1]], r[[2]]))))
    return(points[inds, ])
  }
}