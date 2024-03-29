% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{wind}
\alias{wind}
\title{Saturna Island wind directions}
\format{
A data frame with 239 rows and 2 columns; the column "angle" provides the angular direction (in radian)
and the column day provides the days on which the data points were collected (ranges between 1-10, corresponding to
October 1-10, 2016).
}
\source{
Environment Canada: \url{https://climate.weather.gc.ca/climate_data/data_quality_e.html}.

CBC news on the storm: \url{https://www.cbc.ca/news/canada/british-columbia/storm-bc-1.3795204}.
}
\usage{
data(wind)
}
\description{
A dataset consisting of 239 observations on wind direction in radians (original measurements were
in 10s of degrees), measured at Saturna Island, British Columbia,
Canada during October 1-10, 2016 (obtained from Environment Canada website). There was a severe storm
during October 4-7, which caused significant fluctuations among the wind directions. As a result the
angular data show a clear multimodality.
}
\keyword{datasets}
