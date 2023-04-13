\name{series_desstat}
\alias{series_descstat}
\title{Descriptive Statistics Of A Series}
\usage{
series_descstat(Y)
}
\arguments{
\item{Y}{Univariate time series}

}
\value{
\itemize{
  \item desc_table: A table contains 10 descriptive statistics row-wise
  }
}
\description{
Provides descriptive statistics of a particular series. First column in the output result mentions 10 different statistics and second column contains the Statistics values of the particular series.
}
\examples{
Y <- rnorm(100, 100, 10)
result <- series_descstat(Y)
}
\references{
\itemize{
  \item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.
  \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian svrual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
  }
}
