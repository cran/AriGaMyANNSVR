\name{series_stationarity}
\alias{series_stationarity}
\title{Stationarity Tests Of A Series}
\usage{
series_stationarity(Y)
}
\arguments{
\item{Y}{Univariate time series}

}
\value{
\itemize{
  \item stationarity_table - List of three data frames: 'ADF', 'PP', 'KPSS'
  }
}
\description{
Provides a list of three data frames: 'ADF', 'PP', 'KPSS'. Also indicates whether the data is stationary or not according to the null hypothesis of the corresponding tests.
}
\examples{
Y <- rnorm(100, 100, 10)
result <- series_stationarity(Y)
}
\references{
\itemize{
  \item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.
  \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian svrual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
  }
}

