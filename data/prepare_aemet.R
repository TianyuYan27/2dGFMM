#!/usr/bin/env Rscript
## Build data/aemet_dat.RData, the station table used by demo_aemet.R.
##
## The raw data come from the Spanish Meteorological Agency (AEMET) OpenData
## service through the climaemet package.  A free API key is required for the
## download; see data/README.md for the registration link and the setup.
##
## Run from the repository root (or anywhere; the script finds its own folder):
##
##   Rscript data/prepare_aemet.R
##
## Responses are cached one RData file per station and period under
## data/cache/, so an interrupted run resumes where it stopped.  Set the
## environment variable AEMET_CACHE_DIR to point the cache somewhere else.
##
## The script writes data/aemet_dat.RData and prints the number of stations
## kept.  Stations for which AEMET returns no records are skipped.

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    getwd()
  }
})

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------

years <- 2003:2023
days_per_year <- 365

## AEMET is queried in two blocks rather than one 21-year request.
periods <- list(
  "data2003-2012" = c(start = "2003-01-01", end = "2012-12-31"),
  "data2013-2023" = c(start = "2013-01-01", end = "2023-12-31")
)

cache_root <- Sys.getenv("AEMET_CACHE_DIR",
                         unset = file.path(script_dir, "cache"))
period_dirs <- file.path(cache_root, names(periods))
out_file <- file.path(script_dir, "aemet_dat.RData")

needed <- c("fda.usc", "climaemet", "lubridate")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install the missing packages first:\n  install.packages(c(",
       paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE)
}
suppressPackageStartupMessages({
  library(fda.usc)
  library(climaemet)
  library(lubridate)
})

## The key is only needed when something still has to be downloaded, so a
## run that resumes from, or is served by, a complete cache works offline.
require_api_key <- function(stations, dirs) {
  needed <- unlist(lapply(dirs, function(dir) {
    file.path(dir, paste0(stations, ".RData"))
  }))
  if (length(needed) == 0 || all(file.exists(needed))) {
    return(invisible(TRUE))
  }
  if (!isTRUE(climaemet::aemet_detect_api_key())) {
    stop("No AEMET API key found.\n",
         "Register at https://opendata.aemet.es/centrodedescargas/altaUsuario\n",
         "and store the key once with\n",
         "  climaemet::aemet_api_key(\"YOUR_KEY\", install = TRUE)\n",
         "See data/README.md for details.", call. = FALSE)
  }
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## Station list and download
## ---------------------------------------------------------------------------

## fda.usc::aemet ships the station code, name, longitude and latitude of the
## Spanish climate stations.  Two stations share the indicativo 1387; AEMET
## uses 1387E for A Coruna/Alvedro, so the code is renamed here.
aemet_station_lookup <- function() {
  env <- new.env()
  utils::data("aemet", package = "fda.usc", envir = env)
  lookup <- env$aemet$df
  lookup$ind[lookup$name == "A CORUÑA/ALVEDRO"] <- "1387E"
  lookup
}

download_period <- function(station, period, dir) {
  path <- file.path(dir, paste0(station, ".RData"))
  if (file.exists(path)) {
    return(invisible(path))
  }
  message("Downloading station ", station, " (", basename(dir), ")")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (attempt in 1:3) {
    ok <- tryCatch({
      data_daily <- climaemet::aemet_daily_clim(
        station = station,
        start = period[["start"]],
        end = period[["end"]],
        progress = FALSE
      )
      save(data_daily, file = path)
      TRUE
    }, error = function(e) {
      message("    ", conditionMessage(e))
      FALSE
    })
    if (ok) {
      return(invisible(path))
    }
    Sys.sleep(5)
  }
  warning("No data for station ", station, " in ", basename(dir),
          "; the station will be skipped.", call. = FALSE)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## Tidy one station-year
## ---------------------------------------------------------------------------

## Put the days of one year on a fixed 365-day calendar (29 February is
## removed in leap years), turn trace precipitation "Ip" into 0.05, and use
## 0.01 for a reported zero so that the value can be logged.
align_station_year <- function(prec, temp, wind, dates, year) {
  full_year_dates <- seq.Date(as.Date(paste0(year, "-01-01")),
                              as.Date(paste0(year, "-12-31")),
                              by = "day")
  missing_dates <- as.Date(setdiff(full_year_dates, as.Date(dates)))
  if (length(missing_dates) != 0) {
    names(prec) <- names(temp) <- dates
    extended_prec <- setNames(rep(NA, length(full_year_dates)), full_year_dates)
    extended_prec[names(prec)] <- prec
    extended_temp <- setNames(rep(NA, length(full_year_dates)), full_year_dates)
    extended_temp[names(temp)] <- temp
    extended_wind <- setNames(rep(NA, length(full_year_dates)), full_year_dates)
    extended_wind[names(wind)] <- wind
    prec <- extended_prec
    temp <- extended_temp
    wind <- extended_wind
  }
  if (year %in% c(2004, 2008, 2012, 2016, 2020)) {
    prec <- prec[-60]
    temp <- temp[-60]
    wind <- wind[-60]
  }
  prec[which(prec == "Ip")] <- "0,05"
  # "Acum" marks a value accumulated over several days; it is left missing and
  # interpolated later by prepare_latitude_summer().
  prec[which(prec == "Acum")] <- NA
  prec <- as.numeric(gsub(",", ".", prec))
  prec <- ifelse(prec == 0, 0.01, prec)
  list(prec = prec, temp = temp, wind = wind)
}

## AEMET does not report the same variables for every station and period
## (for instance, wind speed is absent for some stations after 2012).  Keep
## the variables the model uses and fill a variable that is absent in one
## period with NA before stacking the two periods.
merge_periods <- function(pieces) {
  keep <- c("fecha", "indicativo", "nombre", "provincia", "altitud",
            "tmed", "prec", "velmedia")
  pieces <- lapply(pieces, function(d) {
    for (nm in setdiff(keep, names(d))) {
      d[[nm]] <- NA
    }
    d[, keep, drop = FALSE]
  })
  do.call(rbind, pieces)
}

## Stack the 21 fixed-length years into a single 7665-day vector with the
## earliest year first.
station_record <- function(data_daily, station, lookup, years = 2003:2023) {
  if (is.null(data_daily) || nrow(data_daily) == 0) {
    return(NULL)
  }
  precipitation <- temperature <- windspeed <- NULL
  for (y in years) {
    data_sub <- data_daily[lubridate::year(data_daily$fecha) == y, ]
    ## A period that could not be downloaded leaves whole years empty; keep
    ## the 21-year layout with NA so that the vectors stay aligned.
    if (nrow(data_sub) == 0) {
      precipitation <- c(rep(NA_real_, days_per_year), precipitation)
      temperature <- c(rep(NA_real_, days_per_year), temperature)
      windspeed <- c(rep(NA_real_, days_per_year), windspeed)
      next
    }
    wind <- if (is.null(data_sub$velmedia)) rep(NA, days_per_year) else data_sub$velmedia
    aligned <- align_station_year(data_sub$prec, data_sub$tmed, wind,
                                  data_sub$fecha, y)
    precipitation <- c(aligned$prec, precipitation)
    temperature <- c(aligned$temp, temperature)
    windspeed <- c(aligned$wind, windspeed)
  }
  meta <- lookup[lookup$ind == station, ]
  data.frame(
    index = unique(data_daily$indicativo),
    station_name = unique(data_daily$nombre),
    province = unique(data_daily$provincia),
    altitude = unique(data_daily$altitud),
    longitude = meta$longitude,
    latitude = meta$latitude,
    precipitation = I(t(precipitation)),
    temperature = I(t(temperature)),
    windspeed = I(t(windspeed)),
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## Main
## ---------------------------------------------------------------------------

main <- function() {
  lookup <- aemet_station_lookup()
  stations <- lookup$ind

  ## Download.  Existing cache files are reused, so this is cheap to restart.
  for (dir in period_dirs) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  require_api_key(stations, period_dirs)
  for (station in stations) {
    for (i in seq_along(periods)) {
      download_period(station, periods[[i]], period_dirs[i])
    }
  }

  ## Combine the two periods per station and build the table.
  dat <- data.frame()
  for (station in stations) {
    pieces <- lapply(period_dirs, function(dir) {
      path <- file.path(dir, paste0(station, ".RData"))
      if (!file.exists(path)) return(NULL)
      env <- new.env()
      load(path, envir = env)
      env$data_daily
    })
    pieces <- Filter(function(x) !is.null(x) && nrow(x) > 0, pieces)
    if (length(pieces) == 0) {
      next
    }
    data_daily <- merge_periods(pieces)
    record <- station_record(data_daily, station, lookup, years)
    if (!is.null(record)) {
      dat <- rbind(record, dat)
    }
  }
  dat <- dat[order(dat$latitude), ]
  rownames(dat) <- NULL

  save(dat, file = out_file)
  message("\nSaved ", nrow(dat), " stations x ", ncol(dat), " fields to ",
          out_file)
  invisible(dat)
}

if (sys.nframe() == 0) {
  main()
}
