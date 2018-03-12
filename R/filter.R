#' Main filter step m1
#'
#' Flag values with FALSE that are missing or belong to a p_id having lon lat
#' values which occur more than the 'cutOff' value.
#'
#' @param data data set formated as the sample data (netatmoBer)
#' @param cutOff how much stations are allowed to have the same coordinates,
#'   default is 1. This means that if two p_ids share the same lat lon values,
#'   data at these stations are set to FALSE.
#'
#' @return data.table
#' @export
#'
#' @examples
#' data(netatmoBer)
#' x <- m1(netatmoBer)
m1 <- function(data, cutOff = 1){
  val  <- data[!is.na(ta),.(a = 1), by = .(p_id,lon,lat)]
  bad_s  <- val[,.(anz = sum(lon == val$lon & lat == val$lat)), by = p_id]
  bad_s  <- bad_s[anz > cutOff,]$p_id
  data[,m1:= T]
  data[p_id %in% bad_s,"m1"] <- F
  data[is.na(ta),"m1"] <- F
  return(data)
}

#' Version of the Qn function respecting NaN values
#'
#' @param x
#'
#' @return NaN for no valid values, otherwise robust scale estimator 'Qn' without NaN
Qnr <- function(x){
  x <- x[!is.na(x)]
  return(robustbase::Qn(x))
}

#' Calculate robust z-score
#'
#' @param x vector to calculate the robust z-score from
#'
#' @return vector with z-score for elements in x or NaN
getZ <- function(x){
  q <- Qnr(x)
  if(is.na(q)){
    return(NaN)
  }
  res <- (x / q) - (median(x, na.rm = T)/q)
  return(res)
}

#' Main filter step m2
#'
#' Flags all values as FALSE if robust z-score is not within the critical values
#' obtained from low and high. This approach is based on the distribution
#' function of all values. !Careful: if number of stations is < 200, it would be
#' better to use the Student t distribution (to be implemented)!
#'
#' @param data data.table object obtained from m1
#' @param low 0 < low < high < 1
#' @param high 0 < low < high < 1
#' @param heightCorrection if set to true (default) and the column "z" exists in
#'   the input data, the temperatures used in calculating the the z-score are
#'   corrected. The applied formular is t_cor = t - ((0.0065 * (z - mz)) where
#'   mz is the spatial mean at each time step
#' @param debug set to true to keep intermediate results
#'
#' @return data.table
#' @export
#'
#' @examples
#' data(netatmoBer)
#' x <- m1(netatmoBer)
#' y <- m2(x)
#'
m2 <-function(data, low = 0.01, high = 0.95, heightCorrection = T, debug = F){
  data[,rem_ta := ta]
  # ensures that all what is wrong in m1 is wrong in m2 too
  data[!m1, "rem_ta"] <- NaN
  if(heightCorrection & "z" %in% colnames(data)){
    agg <- data[,.(mz = mean(z, na.rm = T), by=time)]
    data <- merge(data,agg)
    data[, rem_ta := rem_ta - (0.0065 * (z - mz))]
  }
  data[, z_ta := getZ(rem_ta), by = time]
  data[, m2 := T]
  data[z_ta < qnorm(low) | z_ta > qnorm(high) | is.nan(z_ta), "m2"] <- F
  if(!debug){
    data$rem_ta <- NULL
    data$z_ta <- NULL
    data$mz <- NULL
  }
  return(data)
}

#' Monthly correlation of individual CWS with median
#'
#' Calculates the correlation of each CWS vs. a data.set containing an
#' aggregated time series in the column 'med' per month. Both series need to
#' have the same length and values at the same position are expected to belong
#' to the same position in time. Function is called internally by filter m4.
#'
#' @param x values of unaggregated time series
#' @param y data.table containing column 'med' holding aggregated values and
#'   month holding the month the time series belongs to
#' @param m month to base the calculation on
#' @param cutOff value below which FALSE is returned.
#'
#' @return TRUE if correlation for the given month is higher than cutOff, FALSE
#'   otherwise
#'
#' @examples
#' see m4
cor_month <- function(x, y, m, cutOff){
  if(length(x) != length(y[month == m,]$med)){
    stop("Dimensions are off, are you sure your data set contain an NaN value for each p_id at each missing time step?")
  }
  c <- suppressWarnings(cor(x, y[month == m,]$med, use="pairwise.complete.obs")) #supress warning if no pairwise complete obs exists just return false
  #maybe move out of loop ?
  if(is.na(c)){
    return(F) #treat na as no corrleation
  }
  if(c < cutOff ){
    return(F)
  }
  return(T)
}

#' Main filter step m3
#'
#' Flag values with FALSE if more than cutOff percent values are removed during
#' m2 per month. This is done since it is assumed that if too many individual
#' values are flagged FALSE in m2, the station is too suspicious to be kept.
#'
#' @param data data.table object obained from m2
#' @param cutOff value above which data are flagged with FALSE, 0 < cutOff < 1.
#'   Default is 0.2, i.e., 20 percent of data.
#'
#' @return data.table
#' @export
#'
#' @examples
#' y <- m2(x)
#' z <- m3(y)
m3 <- function(data, cutOff = 0.2){
  has_m <- has_month(data)
  if(!has_m){
    data[,month := lubridate::floor_date(time,"month")]
  }
  data[, m3 := m2 & sum(m2)/.N > cutOff, by = .(month, p_id)]
  if(!has_m){
    data$month <- NULL
  }
  return(data)
}

#' Main filter m4
#'
#' Flag values with FALSE if they belong to a month in which the correlation
#' with the median of all stations is lower than cutOff.
#'
#' @param data data.table as returned by m3
#' @param cutOff value of correlation coefficient below which data are flagged
#'   with FALSE, 0 < cutOff < 1. Default is 0.9.
#'
#' @return data.table
#' @export
#'
#' @examples
#' y <- m3(x)
#' z <- m4(y)
m4 <- function(data, cutOff = 0.9){
  has_m <- has_month(data)
  if(!has_m){
    data[,month := lubridate::floor_date(time,"month")]
  }
  data[,rem_ta := ta]
  data[!m3, "rem_ta"] <- NaN
  data_agg <- data[,.(med = median(rem_ta, na.rm = T)), by=.(month, time)]
  cor_station <- data[,.(c = cor_month(rem_ta, data_agg, unique(month), cutOff)), by = .(month, p_id)]
  data <- merge(data, cor_station, all.x = T, by = c("month", "p_id"))
  data[, m4 := c & m3]
  data$c <- NULL
  data$rem_ta <- NULL
  if(!has_m){
    data$month <- NULL
  }
  return(data)
}

#' Interpolation
#'
#' This function takes a numerical vector x and fills NaNs with linearly
#' interpolated values. The allowed length of the gap, i.e., the number of
#' consecutive NaNs to be interpolated and replaces is smaller or equal
#' maxLength. Internally called by o1.
#'
#' @param x a numeric vector
#' @param maxLength allowed length of the gap to interpolate, default is 1.
#'
#' @return vector
#' @export
#'
#' @examples
#' x <- x(1, NaN, 3, NaN NaN, 6, NaN, 8)
#' interpol(x)
#' interpol(x, 2)
interpol <- function(x, maxLength = 1){
  rl <- rle(is.na(x))
  re <- rl$lengths <= maxLength & rl$values
  idx <- rep(re,rl$lengths)
  if(sum(!is.na(x)) < 2){
    return(x)
  }
  x[idx] <- approx(1:length(x), x, (1:length(x))[idx])$y
  return(x)
}

#' Optional filter step o1
#'
#' In this step missing data is interpolated, default is to perform linear
#' interpolation on gaps of maximal length = 1.
#'
#' @param data data.table as returned from m4
#' @param fun  function to use for interpolation, default is interpol
#' @param ...  additional parameters for interpolation function
#'
#' @return data.table
#' @export
#'
#' @examples
#' #default
#' o_1 <- o1(m_4)
#' #interpolate gaps up to 5 hours
#' o_1 <- o1(m_4, maxLength = 5)
o1 <- function(data, fun = interpol, ...){
  data[,ta_int := ta]
  data[!m4, "ta_int"] <- NA
  data[,ta_int := fun(ta_int, ...), by = .(p_id)]
  data[, o1:= (is.na(ta) & !is.na(ta_int)) | m4] #has been interpolated or not
  data[, ta := ta_int]
  data$ta_int <- NULL
  return(data)
}

#' Optional filter step o2
#'
#' Optional filter for temporal data availability. Flags all values in a
#' calendar day as FALSE if less than 'cutOff' percent of valid values are
#' available for that day.
#'
#' @param data data.table as returned from o1
#' @param cutOff percentage of values that must be present for each day before
#'   all values of that day are flagged with FALSE, expressed in fraction: 0 <
#'   cutOff < 1. Default is 0.8, i.e, 80 percent of data.
#'
#' @return data.table
#' @export
#'
#' @examples
#' o_2 <- o2(o_1)
o2 <- function(data, cutOff = 0.8){
  has_d <- "day" %in% colnames(data)
  if(!has_d){
    data[, day := lubridate::floor_date(time,"day")]
  }
  data[, o2 := o1 & sum(o1)/.N < cutOff, by = .(day, p_id)]
  if(!has_d){
    data$day <- NULL
  }
  return(data)
}

#' Optional filter step o3
#'
#' Optional filter for temporal data availability. Flags all values in a month
#' as FALSE if less than 'cutOff' percent of valid values are available for that
#' month.
#'
#' @param data data.table as returned from o2
#' @param cutOff percentage of values that must be present for each month before
#'   all values of that month are flagged with FALSE, expressed in fraction: 0 <
#'   cutOff < 1. Default is 0.8, i.e, 80 percent of data.
#'
#' @return data.table
#' @export
#'
#' @examples
#' o_3 <- o3(o_2)
o3 <- function(data, cutOff = 0.8){
  has_m <- has_month(data)
  if(!has_m){
    data[, month := lubridate::floor_date(time,"month")]
  }
  data[, o3 := o2 & sum(o2)/.N < cutOff, by = .(month, p_id)]
  if(!has_m){
    data$month <- NULL
  }
  return(data)
}

#' has_month
#'
#' @param data
#'
#' @return true if data contains a column month
has_month <- function(data){
  return( "month" %in% colnames(data))
}


#' Complete quality control (QC) of CWS data
#'
#' Performs all QC/filter steps in consecutive order. All settings are according
#' to Napoly et al. (2018). This is the default function to carry out the
#' complete QC procedure. Each filter step takes the result of the previous
#' filter step as input. Thus, e.g., when applying filter step m2, the column
#' 'm1' must be in the input data, for filter step m3 the column 'm2' must be
#' present, etc. Each individual step could be skipped by renaming the columns
#' in the input data. After each filter step a new column of type BYTE is
#' included in the output, containing TRUE or FALSE flags. Flags of the previous
#' levels are carried along, i.e., if a value failed in step m2, this FALSE is
#' kept throughout the remianing filter steps. In the end, only those data
#' values containing TRUE after the all filter steps are valid according to this
#' QC.
#'
#' @param data input data in the format as the example data (netatmoBer)
#' @param m1_cutOff see cutOff in ?m1
#' @param m2_low see low in ?m2
#' @param m2_high see high in ?m2
#' @param m3_cutOff see cutOff in ?m3
#' @param m4_cutOff see cutOff in ?m4
#' @param o1_fun see fun in ?o1
#' @param o2_cutOff see cutOff in ?o2
#' @param o3_cutOff see cutOff in ?o3
#' @param includeOptional set to TRUE if the filter steps o1 to o3 shall be
#'   performed, default: TRUE
#' @param ... additional parameters used in o1. For details see ?o1
#'
#' @return data.table
#' @export
#'
#' @examples
#' data(netatmoBer)
#' y <- filterCWS(netatmoBer)
filterCWS <- function(data,
                   m1_cutOff = 1,
                   m2_low = 0.1, m2_high = 0.95,
                   m3_cutOff = 0.2,
                   m4_cutOff = 0.9,
                   o1_fun = interpol,
                   o2_cutOff = 0.8,
                   o3_cutOff = 0.8,
                   includeOptional = T,
                   ...){
    data[, month := lubridate::floor_date(time,"month")]
    data <- m1(data, cutOff = m1_cutOff)
    data <- m2(data, low = m2_low, high = m2_high)
    data <- m3(data, cutOff = m3_cutOff)
    data <- m4(data, cutOff = m4_cutOff)
    if(includeOptional){
      data <- o1(data, fun = o1_fun, ...)
      data <- o2(data, cutOff = o2_cutOff)
      data <- o3(data, cutOff = o3_cutOff)
    }
    data$month <- NULL
    return(data)
}
