% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{m3}
\alias{m3}
\title{m3}
\usage{
m3(data, cutOff = 0.2)
}
\arguments{
\item{data}{data.table object obained from m2}

\item{cutOff}{}
}
\value{
data.table
}
\description{
flag values with fales if more than cutOff percent values are removed during
m1. Steps can be skipped by renaming columns in the input data.
}
\examples{
y <- m2(x)
z <- m3(y)
}