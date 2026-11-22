spain_annual_covariates <- function(years = 2003:2023) {
  # World Bank NY.GDP.PCAP.CD and EN.ATM.CO2E.PC for Spain, in year order.
  data.frame(
    year = years,
    gdppc = c(
      21510.83625, 24907.00085, 26429.15094, 28389.07858, 32591.35031,
      35510.72223, 32169.50285, 30532.48051, 31677.90031, 28322.94659,
      29077.18206, 29513.65118, 25754.36103, 26537.15949, 28185.32137,
      30379.72111, 29581.51855, 26984.29628, 30488.82095, 29674.54429,
      32676.98371
    ) / 10000,
    co2_ems = c(
      7.553341743, 7.817928356, 8.029197431, 7.69812925, 7.842320847,
      7.056479349, 6.214687048, 5.885763493, 5.893329176, 5.77854921,
      5.229400787, 5.206009215, 5.538855216, 5.340585339, 5.681590779,
      5.521059835, 5.131798801, 4.279594571, NA, NA, NA
    )
  )
}

# Each station row stores years from most recent to earliest. Prepending the
# 365-day blocks puts 2003 in the first rows, in line with rep(years, each = n).
stack_years <- function(values, n_year = 21, days_per_year = 365) {
  stacked <- NULL
  for (i in seq_len(n_year)) {
    block <- values[, (days_per_year * (i - 1) + 1):(days_per_year * i)]
    stacked <- rbind(block, stacked)
  }
  if (anyNA(stacked)) zoo::na.approx(stacked, na.rm = FALSE, rule = 2) else stacked
}

prepare_latitude_summer <- function(dat, summer_days = 153:273, years = 2003:2023) {
  stations <- dat[dat$altitude <= 650 & dat$altitude >= 30 & dat$longitude >= -13, ]
  stations <- stations[order(stations$latitude), ]
  n_station <- nrow(stations)
  n_year <- length(years)
  n <- n_year * n_station
  temperature <- stack_years(stations[, "temperature"], n_year)[, summer_days]
  precipitation <- stack_years(stations[, "precipitation"], n_year)
  windspeed <- stack_years(stations[, "windspeed"], n_year)
  annual <- spain_annual_covariates(years)
  # Length-21 covariates are repeated in data.frame order across the
  # n_year * n_station rows, matching the original data.frame() call.
  data <- data.frame(
    ID = rep(seq_len(n_year), each = n_station),
    year = rep(years, each = n_station),
    gdppc = rep(annual$gdppc, length.out = n),
    co2_ems = rep(annual$co2_ems, length.out = n),
    whether_late = factor(ifelse(rep(years, each = n_station) <= 2008, 0, 1)),
    max_prec = apply(precipitation[, summer_days], 1, function(x) max(x, na.rm = TRUE)),
    max_wind = apply(windspeed[, summer_days], 1, function(x) max(x, na.rm = TRUE)),
    prec = I(precipitation[, summer_days]),
    wind = I(windspeed[, summer_days]),
    temperature = I(temperature)
  )
  list(data = data, S = n_station)
}

fit_temperature <- function(formula, data, S, knots,
                            bootstrap.opt = list(B = 50, M = 10000, method = "nonparametric")) {
  fit <- gfmm2d_gaussian(
    formula = formula,
    data = data,
    S = S,
    robust = FALSE,
    smoother = "te",
    knots = knots,
    bootstrap.opt = bootstrap.opt,
    parallel = FALSE,
    pcb = TRUE,
    scb = TRUE
  )
  print(fit$qn)
  print(pvalues.gfmm2d(fit, const = 0, tol = 0.001, delta = 0))
  fit
}
