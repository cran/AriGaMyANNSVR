#' @title Descriptive Statistics Of A Series
#' @description Provides descriptive statistics of a particular series. First column in the output result mentions 10 different statistics and second column contains the Statistics values of the particular series.
#' @param Y Univariate series
#' @import dplyr, psych, e1071, stats
#'
#' @return
#' \itemize{
#'   \item desc_table - A table contains 10 descriptive statistics row-wise
#' }
#' @export
#'
#' @examples
#' Y <- rnorm(100, 100, 10)
#' result <- series_descstat(Y)
#' @references
#' \itemize{
#'\item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.

#' \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.

#' }
series_descstat <- function(Y){
  sw_test <- shapiro.test(Y)
  p_value <- sw_test$p.value
  sw_stat <- round(sw_test$statistic, 3)

  if (p_value <= 0.01) {
    shapiro_wilk <- paste0("***")
  } else if (p_value <= 0.05) {
    shapiro_wilk <- paste0("**")
  } else if (p_value <= 0.1) {
    shapiro_wilk <- paste0("*")
  } else {
    shapiro_wilk <- paste0("")
  }
  # embedding
  diff.1 <- embed(Y, 2)
  Yt <- diff.1[,1]
  Yt_1 <- diff.1[,2]
  y <- log(Yt/Yt_1)
  stats <- rbind(length(Y),
                 round(min(Y), 3),
                 round(max(Y), 3),
                 round(mean(Y), 3),
                 round(sd(Y), 3),
                 round(sd(y), 3),
                 round((sd(Y) / mean(Y) * 100), 3),
                 round((psych::skew(Y)), 3),
                 round(e1071::kurtosis(Y), 3),
                 paste0(sw_stat, shapiro_wilk))

  desc_table <- data.frame(cbind(c('N', 'Minimum', 'Maximum', 'Mean', 'SD', 'Cond_SD', 'CV(%)',
                                   'Skewness', 'Kurtosis', 'Shapiro-Wilk'), stats))
  colnames(desc_table) <- c('Statistics', 'Values')
  return(desc_table)
}

#'@title Non linearity test of a Data Frame
#' @description Performs non linearity test result for a series. Provides output as a single element (data frame) list. First column mentions different statistics (eps). Other columns are the Statistics values of the particular dimension.
#' @param df Data Frame with first column as serial number or date
#' @import DescribeDF
#'
#' @return
#' \itemize{
#'   \item nonlinearity_list - A list with a single element (data frame) . Element is named as the name of the series provided. The element is such that first column mentions different statistics and other columns are the Statistics values of the particular dimension.
#' }
#' @export
#'
#' @examples
#' my_series <- rnorm(100, 100, 10)
#' nonlinearity <- series_nonlinearity(my_series)
#' nonlinearity$my_series
#' @references
#' \itemize{
#'\item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.

#' \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., ... & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.

#' }
series_nonlinearity <- function(Y){
  sl_no <- 1:length(Y)
  df <- data.frame(sl_no = sl_no, X = Y)
  colnames(df)[2] <- paste(deparse(substitute(Y)))
  nonlinearity_list <- DescribeDF::df_nonlinearity(df)
  return(nonlinearity_list)
}

#'@title Stationarity Tests Of A Series
#' @description Provides a list of three data frames: 'ADF', 'PP', 'KPSS'. Also indicates whether the data is stationary or not according to the null hypothesis of the corresponding tests.
#' @param Y Univariate time series
#' @import DescribeDF
#'
#' @return
#' \itemize{
#'   \item stationarity_table - List of three data frames: 'ADF', 'PP', 'KPSS'
#' }
#' @export
#'
#' @examples
#' my_series <- rnorm(100, 100, 10)
#' series_stationarity(my_series)
#' @references
#' \itemize{
#'\item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.

#' \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., ... & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.

#' }
series_stationarity <- function(Y){
  sl_no <- 1:length(Y)
  df <- data.frame(sl_no = sl_no, X = Y)
  colnames(df)[2] <- paste(deparse(substitute(Y)))
  stationarity_table <- DescribeDF::df_stationarity(df)
  return(stationarity_table)
}

#' @title ARIMA-GARCH Hybrid Modeling
#' @description First fits the time series data by using ARIMA model. If the residuals are having "arch" effect, then GARCH is fitted. Based on the previously mentioned condition final prediction is obtained. More can be found from Paul and Garai (2021)  <doi:10.1007/s00500-021-06087-4>.
#' @param Y Univariate time series
#' @param ratio Ratio of number of observations in training and testing sets
#' @param n_lag  Lag of the provided time series data
#' @import DescribeDf
#' @return
#' \itemize{
#'   \item Output_ariga: List of three data frames: predict_compare, forecast_compare, and metrics
#'   }
#' @export
#'
#' @examples
#' Y <- rnorm(100, 100, 10)
#' result <- ariga(Y, ratio = 0.8, n_lag = 4)
#' @references
#' \itemize{
#'   \item Paul, R. K., & Garai, S. (2021). Performance comparison of wavelets-based machine learning technique for forecasting agricultural commodity prices. Soft Computing, 25(20), 12857-12873.
#'   \item Paul, R. K., & Garai, S. (2022). Wavelets based artificial neural network technique for forecasting agricultural prices. Journal of the Indian Society for Probability and Statistics, 23(1), 47-61.
#'   \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
#'   }

ariga <- function(Y, ratio = 0.9, n_lag = 4){
  # embedding for finding log return

  diff.1 <- embed(Y, 2)
  Yt <- diff.1[,1]
  Yt_1 <- diff.1[,2]
  y <- log(Yt/Yt_1)
  # Compute the average of non-zero contents in the data
  nonzero_data <- y[y != 0 & !is.na(y)]
  average_nonzero <- mean(nonzero_data)

  # Replace NaNs with the average of non-zero contents in the data
  y[is.nan(y)] <- average_nonzero

  # Check the result
  y
  # embedding for finding lag series of actual series
  embed_size_y <- n_lag+1 # same lag (lags_y-1) for every sub series
  diff.2 <- embed(Yt,embed_size_y)
  dim(diff.2)
  Y_actual <- diff.2[,1]
  Y_actual_1 <- diff.2[,2]
  # train-test split
  n <- length(Y_actual) # this is already (original length-embed_y)
  Y_train <- Y_actual[1:(n*ratio)]
  Y_test <- Y_actual[(n*ratio+1):n]
  Y_train_1 <- Y_actual_1[1:(n*ratio)]
  Y_test_1 <- Y_actual_1[(n*ratio+1):n]
  # embedding for finding lag series of log return series
  diff.3 <- embed(y, embed_size_y)
  y_actual <- diff.3[,1]
  y_train <- y_actual[1:(n*ratio)]
  y_test <- y_actual[(n*ratio+1):n]

  # ARIMA
  arima_ori <- forecast::auto.arima(y_train)
  order_ori <- forecast::arimaorder(arima_ori)
  model_arima_ori<-arima(y_train,order=c(order_ori[1], order_ori[2], order_ori[3]))
  pred_arima_ori <-arima_ori$fitted
  forecast_arima_ori <- data.frame(predict(arima_ori,n.ahead=(n-(n*ratio))))
  forecast_arima_ori <-forecast_arima_ori$pred
  #GARCH Model ####
  ARCH_pvalue_ori <- as.numeric(FinTS::ArchTest(model_arima_ori$residuals)$p.value)
  #ARCH_pvalue<-1
  if(ARCH_pvalue_ori<=0.05){
    garch.fit_ori <- fGarch::garchFit(~garch(1, 1), data = y_train, trace = FALSE)
    pred_V_ori <- garch.fit_ori@fitted
    forecast_V_ori <- predict(garch.fit_ori, n.ahead=(n-(n*ratio)))
    forecast_V_ori  <- forecast_V_ori$meanForecast

    Resid_V_ori <- garch.fit_ori@residuals
    for_resid_ori<-as.ts(y_test-forecast_V_ori)
  }else {
    pred_V_ori <- pred_arima_ori
    forecast_V_ori <- forecast_arima_ori

    Resid_V_ori <- as.ts(model_arima_ori$residuals)
    for_resid_ori<-as.vector(y_test-as.vector(forecast_arima_ori))
  }
  ARIGA_train <- exp(pred_V_ori)*Y_train_1
  ARIGA_test <- exp(forecast_V_ori)*Y_test_1
  # validation_ARIGA
  metrics_ariga_train <- data.frame(AllMetrics::all_metrics(Y_train, ARIGA_train))
  metrics_ariga_test <- data.frame(AllMetrics::all_metrics(Y_test, ARIGA_test))

  metrics <- cbind(metrics_ariga_train$Metrics,
                   as.numeric(metrics_ariga_train$Values),
                   as.numeric(metrics_ariga_test$Values))
  colnames(metrics) <- c('Metrics', 'ARIGA_Train',
                         'ARIGA_Test')
  predict_compare <- data.frame(cbind(train_actual = Y_train,
                                      predicted = ARIGA_train))
  colnames(predict_compare) <- c('train_actual', 'train_predicted')
  forecast_compare <- data.frame(cbind(test_actual = Y_test,
                                       forecast = ARIGA_test))
  colnames(forecast_compare) <- c('test_actual', 'test_predicted')
  Output_ariga <- list(predict_compare=predict_compare,
                       forecast_compare=forecast_compare,
                       metrics=metrics)
  return(Output_ariga)
}

#' @title Specially Designed SVR-Based Modeling
#' @description Fits a ANN model to the uni-variate time series data. The contribution is related to the PhD work of the maintainer. More can be found from Paul and Garai (2021)  <doi:10.1007/s00500-021-06087-4>.
#' @param Y Univariate time series
#' @param ratio Ratio of number of observations in training and testing sets
#' @param n_lag  Lag of the provided time series data
#' @import forecast, AllMetrics, neuralnet
#' @return
#' \itemize{
#'   \item Output_ann: List of three data frames: predict_compare, forecast_compare, and metrics
#'   }
#' @export
#'
#' @examples
#' Y <- rnorm(100, 100, 10)
#' result <- my_ann(Y, ratio = 0.8, n_lag = 4)
#' @references
#' \itemize{
#'   \item Paul, R. K., & Garai, S. (2021). Performance comparison of wavelets-based machine learning technique for forecasting agricultural commodity prices. Soft Computing, 25(20), 12857-12873.
#'   \item Paul, R. K., & Garai, S. (2022). Wavelets based artificial neural network technique for forecasting agricultural prices. Journal of the Indian Society for Probability and Statistics, 23(1), 47-61.
#'   \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
#'   }
my_ann <- function(Y, ratio = 0.9, n_lag = 4){
  # embedding for finding log return

  diff.1 <- embed(Y, 2)
  Yt <- diff.1[,1]
  Yt_1 <- diff.1[,2]
  y <- log(Yt/Yt_1)
  # Compute the average of non-zero contents in the data
  nonzero_data <- y[y != 0 & !is.na(y)]
  average_nonzero <- mean(nonzero_data)

  # Replace NaNs with the average of non-zero contents in the data
  y[is.nan(y)] <- average_nonzero

  # Check the result
  y
  # embedding for finding lag series of actual series
  embed_size_y <- n_lag+1 # same lag (lags_y-1) for every sub series
  diff.2 <- embed(Yt,embed_size_y)
  dim(diff.2)
  Y_actual <- diff.2[,1]
  Y_actual_1 <- diff.2[,2]
  # train-test split

  n <- length(Y_actual) # this is already (original length-embed_y)
  Y_train <- Y_actual[1:(n*ratio)]
  Y_test <- Y_actual[(n*ratio+1):n]
  Y_train_1 <- Y_actual_1[1:(n*ratio)]
  Y_test_1 <- Y_actual_1[(n*ratio+1):n]
  # embedding for finding lag series of log return series
  diff.3 <- embed(y, embed_size_y)
  diff.3_train <- diff.3[1:(n*ratio),]
  diff.3_test <- diff.3[(n*ratio+1):n,]
  # Set the column names
  colnames(diff.3_train) <- c("y_train", paste0("x", 1:(embed_size_y-1)))
  colnames(diff.3_test) <- c("y_test", paste0("x", 1:(embed_size_y-1)))
  # extract the predictors and response from the current matrix
  predictors_train <- diff.3_train[, -1]
  response_train <- diff.3_train[, 1]
  allvars <- colnames(predictors_train)
  predictorvars <- paste(allvars, collapse="+")
  form <- as.formula(paste("y_train ~", predictorvars))

  y_actual <- diff.3[,1]
  y_train <- y_actual[1:(n*ratio)]
  y_test <- y_actual[(n*ratio+1):n]

  # ANN ####
  nn <- neuralnet::neuralnet(form, diff.3_train, hidden = c(4))
  train_predicted <- predict(nn, predictors_train)

  predictors_test <- diff.3_test[, -1]
  response_test <- diff.3_test[, 1]
  test_predicted <- predict(nn, predictors_test)

  ann_train <- exp(train_predicted)*Y_train_1
  ann_test <- exp(test_predicted)*Y_test_1
  # validation_ann
  metrics_ann_train <- data.frame(AllMetrics:: all_metrics(Y_train, ann_train))
  metrics_ann_test <- data.frame(AllMetrics:: all_metrics(Y_test, ann_test))

  metrics <- cbind(metrics_ann_train$Metrics,
                   as.numeric(metrics_ann_train$Values),
                   as.numeric(metrics_ann_test$Values))
  colnames(metrics) <- c('Metrics', 'ANN_Train',
                         'ANN_Test')
  predict_compare <- data.frame(cbind(train_actual = Y_train,
                                      predicted = ann_train))
  colnames(predict_compare) <- c('train_actual', 'train_predicted')
  forecast_compare <- data.frame(cbind(test_actual = Y_test,
                                       forecast = ann_test))
  colnames(forecast_compare) <- c('test_actual', 'test_predicted')
  Output_ann <- list(predict_compare=predict_compare,
                     forecast_compare=forecast_compare,
                     metrics=metrics)
  return(Output_ann)
}

#' @title Specially Designed SVR-Based Modeling
#' @description Fits a SVR model to the uni-variate time series data. The contribution is related to the PhD work of the maintainer. More can be found from Paul and Garai (2021)  <doi:10.1007/s00500-021-06087-4>.
#' @param Y Univariate time series
#' @param ratio Ratio of number of observations in training and testing sets
#' @param n_lag  Lag of the provided time series data
#' @import forecast, AllMetrics, e1071
#' @return
#' \itemize{
#'   \item Output_ann: List of three data frames: predict_compare, forecast_compare, and metrics
#'   }
#' @export
#'
#' @examples
#' Y <- rnorm(100, 100, 10)
#' result <- my_svr(Y, ratio = 0.8, n_lag = 4)
#' @references
#' \itemize{
#'   \item Paul, R. K., & Garai, S. (2021). Performance comparison of wavelets-based machine learning technique for forecasting agricultural commodity prices. Soft Computing, 25(20), 12857-12873.
#'   \item Paul, R. K., & Garai, S. (2022). Wavelets based artificial neural network technique for forecasting agricultural prices. Journal of the Indian Society for Probability and Statistics, 23(1), 47-61.
#'   \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
#'   }
my_svr <- function(Y, ratio = 0.9, n_lag = 4){
  # embedding for finding log return

  diff.1 <- embed(Y, 2)
  Yt <- diff.1[,1]
  Yt_1 <- diff.1[,2]
  y <- log(Yt/Yt_1)
  # Compute the average of non-zero contents in the data
  nonzero_data <- y[y != 0 & !is.na(y)]
  average_nonzero <- mean(nonzero_data)

  # Replace NaNs with the average of non-zero contents in the data
  y[is.nan(y)] <- average_nonzero

  # Check the result
  y
  # embedding for finding lag series of actual series
  embed_size_y <- n_lag+1 # same lag (lags_y-1) for every sub series
  diff.2 <- embed(Yt,embed_size_y)
  dim(diff.2)
  Y_actual <- diff.2[,1]
  Y_actual_1 <- diff.2[,2]
  # train-test split

  n <- length(Y_actual) # this is already (original length-embed_y)
  Y_train <- Y_actual[1:(n*ratio)]
  Y_test <- Y_actual[(n*ratio+1):n]
  Y_train_1 <- Y_actual_1[1:(n*ratio)]
  Y_test_1 <- Y_actual_1[(n*ratio+1):n]
  # embedding for finding lag series of log return series
  diff.3 <- embed(y, embed_size_y)
  diff.3_train <- diff.3[1:(n*ratio),]
  diff.3_test <- diff.3[(n*ratio+1):n,]
  # Set the column names
  colnames(diff.3_train) <- c("y_train", paste0("x", 1:(embed_size_y-1)))
  colnames(diff.3_test) <- c("y_test", paste0("x", 1:(embed_size_y-1)))
  # extract the predictors and response from the current matrix
  predictors_train <- diff.3_train[, -1]
  response_train <- diff.3_train[, 1]
  allvars <- colnames(predictors_train)
  predictorvars <- paste(allvars, collapse="+")
  form <- as.formula(paste("y_train ~", predictorvars))

  y_actual <- diff.3[,1]
  y_train <- y_actual[1:(n*ratio)]
  y_test <- y_actual[(n*ratio+1):n]

  # SVR ####
  # train a svm model with radial basis function as kernel
  svmmodel <- e1071::svm(form, diff.3_train, type="eps-regression", kernel = 'radial')

  # predict the response for the current matrix in the training set
  train_predicted <- as.numeric(predict(svmmodel, predictors_train))

  predictors_test <- diff.3_test[, -1]
  response_test <- diff.3_test[, 1]

  # predict the response for the current matrix in the testing set
  test_predicted <- as.numeric(predict(svmmodel, predictors_test))

  svr_train <- exp(train_predicted)*Y_train_1
  svr_test <- exp(test_predicted)*Y_test_1

  # validation_svr
  metrics_svr_train <- data.frame(AllMetrics:: all_metrics(Y_train, svr_train))
  metrics_svr_test <- data.frame(AllMetrics:: all_metrics(Y_test, svr_test))

  metrics <- cbind(metrics_svr_train$Metrics,
                   as.numeric(metrics_svr_train$Values),
                   as.numeric(metrics_svr_test$Values))
  colnames(metrics) <- c('Metrics', 'SVR_Train',
                         'SVR_Test')
  predict_compare <- data.frame(cbind(train_actual = Y_train,
                                      predicted = svr_train))
  colnames(predict_compare) <- c('train_actual', 'train_predicted')
  forecast_compare <- data.frame(cbind(test_actual = Y_test,
                                       forecast = svr_test))
  colnames(forecast_compare) <- c('test_actual', 'test_predicted')
  Output_svr <- list(predict_compare=predict_compare,
                     forecast_compare=forecast_compare,
                     metrics=metrics)
  return(Output_svr)
}
