\name{series_nonlinearity}
\alias{series_nonlinearity}
\title{Non linearity test of a Data Frame}
\usage{
series_nonlinearity(Y)
}
\arguments{
\item{Y}{Univariate time series}

}
\value{
\itemize{
  \item nonlinearity_list: A list with a single element (data frame) . Element is named as the name of the series provided. The element is such that first column mentions different statistics and other columns are the Statistics values of the particular dimension.
  }
}
\description{
Performs non linearity test result for a series. Provides output as a single element (data frame) list. First column mentions different statistics (eps). Other columns are the Statistics values of the particular dimension.
}
\examples{
my_series <- rnorm(100, 100, 10)
nonlinearity <- series_nonlinearity(my_series)
nonlinearity$my_series
}
\references{
\itemize{
  \item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.
  \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian svrual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
  }
}
