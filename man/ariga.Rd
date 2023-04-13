\name{ariga}
\alias{ariga}
\title{ARIMA-GARCH Hybrid Modeling}
\usage{
ariga(Y, ratio = 0.9, n_lag = 4)
}
\arguments{
\item{Y}{Univariate time series}

\item{ratio}{Ratio of number of observations in training and testing sets}

\item{n_lag}{Lag of the provided time series data}
}
\value{
\itemize{
  \item Output_ariga: List of three data frames containing three data frames namely predict_compare, forecast_compare, and metrics
  }
}
\description{
First fits the time series data by using ARIMA model. If the residuals are having "arch" effect, then GARCH is fitted. Based on the previously mentioned condition final prediction is obtained.
}
\examples{
Y <- rnorm(100, 100, 10)
result <- ariga(Y, ratio = 0.8, n_lag = 4)
}
\references{
\itemize{
  \item Paul, R. K., & Garai, S. (2021). Performance comparison of wavelets-based machine learning technique for forecasting agricultural commodity prices. Soft Computing, 25(20), 12857-12873.
  \item Paul, R. K., & Garai, S. (2022). Wavelets based artificial neural network technique for forecasting agricultural prices. Journal of the Indian Society for Probability and Statistics, 23(1), 47-61.
  \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
  }
}

