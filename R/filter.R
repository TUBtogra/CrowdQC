#' perform filter step m1
#'
#' flag values with false that are missing or belong to a p_id having lon lat values
#' which appear more than cutOff
#'
#' @param data data set formated as netatmoBer
#' @param cutOff how much stations are allowed to have the same coordinates,
#'   default is 1 meaning no p_id's with the same lat lon values are allowed
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
  bad_s  <- bad_s[anz > 1,]$p_id
  data[,m1:= T]
  data[p_id %in% bad_s,"m1"] <- F
  data[is.na(ta),"m1"] <- F
  return(data)
}

#' Version of Qn respecting NaN values
#'
#' @param x
#'
#' @return NaN for no valid values otherwise Qn without NaN
Qnr <- function(x){
  x <- x[!is.na(x)]
  return(robustbase::Qn(x))
}

#' calculate robust z-score
#'
#' @param x vector to caluclate the robust score from
#'
#' @return vector with z-Score for elements in x or NaN
getZ <- function(x){
  q <- Qnr(x)
  if(is.na(q)){
    return(NaN)
  }
  res <- (x / q) - (median(x, na.rm = T)/q)
  return(res)
}

#' perform filter step m2
#'
#' flags all values with robust z-score is not within the critical values
#' obtained from low and high
#' Steps can be skipped by renaming columns in the input data
#'
#' @param data data.table object obtained from m1
#' @param low value > 0
#' @param high value < 1
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
m2 <-function(data, low = 0.01, high = 0.95, debug = F){
  data[,rem_ta := ta]
  # ensures that all what is wrong in m1 is wrong in m2 too
  data[!m1, "rem_ta"] <- NaN
  data[, z_ta := getZ(rem_ta), by = time]
  data[, m2 := T]
  data[z_ta < qnorm(low) | z_ta > qnorm(high) | is.nan(z_ta), "m2"] <- F
  if(!debug){
    data$rem_ta <- NULL
    data$z_ta <- NULL
  }
  return(data)
}

#' cor_month
#'
#' calculates the correlation on monthly basis vs. an data.set containing an
#' aggregated time series in the column "med". Both series need too have the
#' same length and values at the same position are expected to belong to the
#' same position in time.
#'
#' @param x values of unaggregated time series
#' @param y data.table containing column med holding aggregated value and month
#'   holding the month the time series belongs to
#' @param m month to base the calculation on
#' @param cutOff value below which False is returned
#'
#' @return True if correlation for the given month is higher than cutOff false otherwise
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

#' m3
#'
#' flag values with false if more than cutOff percent values are removed during
#' m1. Steps can be skipped by renaming columns in the input data.
#'
#' @param data data.table object obained from m2
#' @param cutOff
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

#' m4
#'
#' flag values with false if they belong to a month in which the correlation with
#' median of all stations is lower than cutOff. Steps can be skipped by renaming
#' columns in the input data.
#'
#' @param data data.table as returned by m3
#' @param cutOff
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

#' interpol
#'
#' This function takes a numerical vector x and fills NaN's with linear
#' interpolated values. But only if the NaN gap is smaller or equal maxLength
#'
#' @param x a numeric vector
#' @param maxLength How long a gap can be before it remains
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

#' step o1 interpolate missing data
#'
#' In this step missing data is interpolated, default is to perform linear
#' interpolation on gaps of maximal 1h. Steps can be skipped by renaming
#' columns in the input data.
#'
#' @param data data.table as returned from m4
#' @param fun  function to use for interpolation default is interpol
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

#' filter step o2
#'
#' For consitency, filter if lesser than 80 percent of values for a day are
#' present. Steps can be skipped by renaming columns in the input data.
#'
#' @param data data.table as returned from o1
#' @param cutOff percentage of values that must be present before all values of
#'   a day are flagged with false
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

#' step o3
#'
#' For consitency, filter if lesser than 80 percent of values for a day are
#' present. Steps can be skipped by renaming columns in the input data.
#'
#' @param data data.table as returned from o2
#' @param cutOff percentage of values that must be present before the month is flagged false
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
#' @return true if data  contains a column month
has_month <- function(data){
  return( "month" %in% colnames(data))
}


#' filter cws data
#'
#' performs all filter steps.
#'
#' @param data input data in a format like data(netatmoBer)
#' @param m1_cutOff see cutoff in ?m1
#' @param m2_low see low in ?m2
#' @param m2_high see high in ?m2
#' @param m3_cutOff see cutoff in ?m3
#' @param m4_cutOff see cutoff in ?m4
#' @param o1_fun see fun in ?o1
#' @param o2_cutOff see cutoff in ?o2
#' @param o3_cutOff see cutoff in ?o3
#' @param includeOptional if the steps o1 to o3 should be performed
#' @param ... additional parameters used in o1 for details see ?o1
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
