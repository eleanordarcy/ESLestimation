#' Sea level data from Lowestoft (UK)
#'
#' Peak tide and skew surge observations (in meters, against chart datum) at Lowestoft (UK). 
#' This is a processed and quality controlled version of that from the British Oceanographic Data Centre, as part of the UK National Tide Gauge Network.
#' 
#' @name Lowestoft
#' 
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(Lowestoft)
#'
#' @rdname Lowestoft-data
#'
#' @format A dataframe with 8 columns and 26848 rows.
#' The first four columns correspond to the observation time, they are year, month, day (day in year) and hour. 
#' Observations are recorded for every tidal cycle, i.e., every 12.5 hours approximately. Although there is some missingness. 
#' The fifth and sixth column are skew surge (\code{skews}) and peak tide (\code{maxTide}) observations, respectively.
#' The seventh column gives a date object for the date of observation.
#' The eighth column is a standardised version of the peak tide column (\code{stTide}), where the standardisation is done by subtracting the mean and dividing by the standard deviation of all peak tides.
#' The original data is documented in \insertCite{DArcy2023;textual}{ESLestimation}
#' 
#' @examples
#' library(ESLestimation)
#' data(Lowestoft)
#'
#' # Look at skew surges over time
#' plot(Lowestoft$datetime, Lowestoft$skews, xlab='Date', ylab='Skew surge (m)')
#' plot(Lowestoft$day, Lowestoft$skews, xlab='Day in year', ylab='Skew surge (m)')
#'
#' Look at peak tides over time
#' plot(Lowestoft$datetime, Lowestoft$maxTide, xlab='Date', ylab='Peak tide (m)')
#' 
#' @references \insertRef{DArcy2023}{ESLestimation}
#'
NULL