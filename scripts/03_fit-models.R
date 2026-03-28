library(forecast)

# fit models
persist_fc <- function(train, holdout, iprint = FALSE) {
  ntrain <- length(train)
  nh <- length(holdout)
  sse <- 0
  fc_vec <- numeric(nh)
  fc <- train[ntrain] 
  yt <- holdout[1]
  fc_vec[1] <- fc
  fcerror <- yt - fc
  sse <- sse + (fcerror^2)
  for (i in 2:nh) {
    yt <- holdout[i]
    fc <- holdout[i-1]
    fc_vec[i] <- fc
    fcerror <- yt - fc
    sse <- sse + (fcerror^2)
  }
  return(list(rmse = sqrt(sse/nh), fc = fc_vec))
}


iid_fc <- function(train, holdout, iprint = FALSE) {
  n_train <- length(train)
  n_holdout <- length(holdout)
  avg_val <- mean(train) 
  sse <- 0
  fc_vec <- rep(avg_val, n_holdout) # Every forecast is the same constant
  for (i in 1:n_holdout) {
    yt <- holdout[i]
    fc <- avg_val
    fcerror <- yt - fc
    sse <- sse + (fcerror^2)
  }
  return(list(rmse = sqrt(sse / n_holdout), fc = fc_vec))
}


train = ts(as.numeric(sub("%","",merged_data$sp500_ret)), start = c(1990,2), end = c(2017,1), frequency = 1)
acf(train)
pacf(train)

auto.arima(train,stationary = T, seasonal = F, stepwise = F)

arima(train,order = c(0,0,3))

