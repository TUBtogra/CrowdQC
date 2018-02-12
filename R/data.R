#' Example netatmo Data
#'
#' A dataset containing example netatmo Data. The data was obtained by hourly
#' issueing the getpublicdata function of the netatmo php api. The requested
#' rectangle was lon ne 13.647242, lat ne 52.727861 and lon sw 12.960249, lat sw
#' 52.340471, split into 64 evenly spaced tiles.
#' Missing values are set to NaN
#'
#' @format A data frame with 1885807 rows and 5 variables: \describe{
#'   \item{time}{hour were the value belongs to}
#'   \item{p_id}{an unique identifier for each netatmo station, changes on station replacement}
#'   \item{ta}{air temperature in degree C}
#'   \item{lon}{longitude}
#'   \item{lat}{latitude}
#' }
"netatmoBer"
