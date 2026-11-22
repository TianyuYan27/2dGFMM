# Spanish daily climate data

`demo_aemet.R` reads `data/aemet_dat.RData`, a station table built from daily
records of the Spanish Meteorological Agency (AEMET). The file is **not
included in this repository**: AEMET asks users to download its data
themselves, and the table can be rebuilt with `prepare_aemet.R`. AEMET answers
each multi-year request in a few minutes, so a full run usually takes several
hours; the cache makes it safe to interrupt and resume. This folder contains the script and documents exactly what
it does, so the demo remains fully reproducible.

## Requirements

* R packages:

  ```r
  install.packages(c("fda.usc", "climaemet", "lubridate", "zoo"))
  ```

  `fda.usc` supplies the station catalogue (code, name, longitude, latitude),
  `climaemet` queries AEMET OpenData, and `lubridate` handles the dates.
  `zoo` is used by the demo when it interpolates missing days.

* A free AEMET OpenData API key. Register at
  <https://opendata.aemet.es/centrodedescargas/altaUsuario>, then store the key
  once:

  ```r
  climaemet::aemet_api_key("YOUR_KEY", install = TRUE, overwrite = TRUE)
  ```

  This writes `AEMET_API_KEY` to your `.Renviron`; restart R afterwards. To
  use a key without storing it, run `Sys.setenv(AEMET_API_KEY = "YOUR_KEY")`
  in the session that runs the script.

## Build the file

From the repository root:

```bash
Rscript data/prepare_aemet.R
```

The script prints `Downloading <station>` for the 73 catalogue stations and
finishes with `Saved 66 stations x 9 fields to .../data/aemet_dat.RData`.

Downloads are cached one RData file per station and period in `data/cache/`.
If the run is interrupted, start it again: cached stations are skipped and
only the missing ones are requested. To use another cache location, set
`AEMET_CACHE_DIR`, for example

```bash
AEMET_CACHE_DIR=/path/to/cache Rscript data/prepare_aemet.R
```

## What the script does

| Step | Operation |
|------|-----------|
| 1. Station list | Read `fda.usc::aemet`, which provides 73 indicators (`ind`), names, longitudes and latitudes. Two stations share the indicator `1387`; the Alvedro airport is stored as `1387E`, so its code is renamed. |
| 2. Download | Call `climaemet::aemet_daily_clim()` for 2003-01-01 to 2012-12-31 and for 2013-01-01 to 2023-12-31 for every station. Each response is saved as `data/cache/<period>/<station>.RData`. |
| 3. Align the calendar | Put every station-year on a fixed 365-day calendar: missing dates are inserted as `NA`, and 29 February is removed in 2004, 2008, 2012, 2016 and 2020. Trace precipitation `"Ip"` becomes `0.05`; a reported zero becomes `0.01` so that it can be logged; `"Acum"` (a multi-day accumulated value) is left missing. |
| 4. Stack | Stack the 21 years of each station into 7665-day vectors, 2003 first. If AEMET reports different variables in the two periods (for instance, wind speed is missing for some stations after 2012), the two blocks are aligned by name and absent variables are filled with `NA`. |
| 5. Save | Keep the variables used by the model (`fecha`, station metadata, `tmed`, `prec`, `velmedia`), drop stations whose downloads are empty, sort by latitude, and save `dat` to `data/aemet_dat.RData`. |

## The resulting object

`aemet_dat.RData` contains a single data frame `dat` with one row per station:

| Column | Content |
|--------|---------|
| `index` | AEMET indicator |
| `station_name`, `province` | Station name and province |
| `altitude` | Station altitude in metres |
| `longitude`, `latitude` | Coordinates in decimal degrees |
| `precipitation` | 66 x 7665 matrix, daily precipitation (mm), one row per station, columns from 2003-01-01 to 2023-12-31 |
| `temperature` | 66 x 7665 matrix, daily mean temperature (degrees Celsius) |
| `windspeed` | 66 x 7665 matrix, daily mean wind speed (m/s) |

`nrow` is 66 rather than 73 because AEMET returns no records for seven
catalogue stations (`1212`, `200`, `76`, `367`, `5530`, `1024`, `3191`); they
are skipped.

## How `demo_aemet.R` uses it

`prepare_latitude_summer()` in `R/helpers/aemet_models.R`:

1. keeps stations with altitude between 30 and 650 m and longitude at or east
   of `-13`, which leaves 31 stations;
2. takes days 153 to 273 (June to September) as the functional index, so the
   grid is 31 latitudes x 121 days;
3. stacks the years, giving one subject per year (`ID`) and `S = 31` visits
   (latitudes);
4. attaches the annual CO<sub>2</sub> emissions and GDP values hard-coded in
   `spain_annual_covariates()`, the post-2008 indicator, the summer maximum
   precipitation and the summer maximum wind speed;
5. interpolates remaining missing days with `zoo::na.approx()`.

The demo then fits `temperature ~ whether_late + co2_ems + wind + max_prec`
with the tensor-product smoother and `knots = c(10, 8)`.

## Notes and troubleshooting

* **API errors.** AEMET occasionally returns a temporary error or throttles
  repeated requests. The script retries each download three times; rerun it to
  continue from the cache. If a long date range is rejected, shorten the
  periods in `periods` (for example to five years each) and keep the same
  folder names.
* **Empty stations.** A station with no records is skipped with the message
  `No data for station ...`. This is expected for the seven indicators listed
  above.
* **Reproducibility.** AEMET revises its archives, so a newly built table may
  differ slightly from the one used in the paper. The estimation and testing
  results reported in the README were obtained from an earlier download; small
  numerical differences do not change the conclusions.
* **Licence.** The data are distributed by AEMET through its OpenData portal.
  Check the reuse conditions at
  <https://opendata.aemet.es/centrodedescargas/inicio> before redistributing
  them.
* **No API key is needed to run the demo** once `aemet_dat.RData` exists; the
  key is only needed to rebuild the table.
