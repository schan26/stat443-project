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

esm_fc <- function(train, holdout, alpha, level, iprint = FALSE) {
  nh <- length(holdout)
  sse <- 0
  fc_vec <- numeric(nh)
  L <- level 
  for (i in 1:nh) {
    fc <- L
    fc_vec[i] <- fc
    yt <- holdout[i]
    fcerror <- yt - fc
    sse <- sse + (fcerror^2)
    L <- alpha * yt + (1 - alpha) * L
  }
  return(list(rmse = sqrt(sse/nh), fc = fc_vec))
}

lholt_fc <- function(train, holdout, alpha, beta, level, slope, iprint = FALSE) {
  nh <- length(holdout)
  sse <- 0
  fc_vec <- numeric(nh)
  L <- level
  b <- slope
  for (i in 1:nh) {
    fc <- L + b
    fc_vec[i] <- fc
    yt <- holdout[i]
    fcerror <- yt - fc
    sse <- sse + (fcerror^2)
    L_old <- L
    L <- alpha * yt + (1 - alpha) * (L_old + b)
    b <- beta * (L - L_old) + (1 - beta) * b
  }
  return(list(rmse = sqrt(sse / nh), fc = fc_vec))
}

arima_fc = function(tsdata,ntrain,order,seasonal,method,traincoef,include.mean,
                    iprint=F)
{obj <<- arima(tsdata,order=order,seasonal=seasonal,init=traincoef,fixed=traincoef, method=method, include.mean=include.mean,
                optim.control=list(maxit=0))
  # yt = y_{t|t-1} + innov ; fc = yt-innov
  fc = tsdata - obj$residuals
  ntotal = length(tsdata)
  holdout_fc = fc[(ntrain+1):ntotal]
  holdout = tsdata[(ntrain+1):ntotal] 
  if(iprint) print(cbind(holdout,holdout_fc)) 
  rmse = sqrt(mean((holdout-holdout_fc)^2)) 
  return(list(rmse=rmse, fc=holdout_fc))
}


data <- read.csv("data/processed/cleanedData.csv")
data_ts <- ts(data$sp500_ret, start = c(1990,2), frequency = 12)

train = window(data_ts, start = c(1990,2), end = c(2015,4)) #train data from 1990-02 to 2015-04
holdout = window(data_ts, start = c(2015,5)) #holdout data from 2015-05 onwards

#Baseline methods - Persistence and IID
persist_results = persist_fc(train, holdout)
iid_results = iid_fc(train, holdout)

#Exponential smoothing methods - Holt linear and simple exponential smoothing
esm_model = HoltWinters(train, beta = F, gamma = F)
esm_results = esm_fc(train, holdout, esm_model$alpha, esm_model$coefficients[1])

lholt_model = HoltWinters(train, gamma = F)
lholt_results = lholt_fc(train, holdout, lholt_model$alpha, lholt_model$beta, lholt_model$coefficients[1], lholt_model$coefficients[2])

#Simple exponential smoothing performs better than linear Holt due to the lack of a consistent trend in the data.

#ARMA
arma_model = auto.arima(train, stationary = T, seasonal = F)
arma_results = arima_fc(data_ts,length(train),c(0,0,0),c(0,0,0),"ML",arma_model$coef,include.mean = T)

#Vary P and Q a little and compare AIC and holdout set RMSE
alt_arma1 = arima(train,order = c(1,0,1))
alt_arma2 = arima(train,order = c(0,0,1))
alt_arma3 = arima(train,order = c(0,0,2))
alt_arma4 = arima(train,order = c(1,0,0))
alt_arma5 = arima(train,order = c(2,0,0))

arma_alt1_results = arima_fc(data_ts,length(train),c(1,0,1),c(0,0,0),"ML",alt_arma1$coef,include.mean = T)
arma_alt2_results = arima_fc(data_ts,length(train),c(0,0,1),c(0,0,0),"ML",alt_arma2$coef,include.mean = T)
arma_alt3_results = arima_fc(data_ts,length(train),c(0,0,2),c(0,0,0),"ML",alt_arma3$coef,include.mean = T)
arma_alt4_results = arima_fc(data_ts,length(train),c(1,0,0),c(0,0,0),"ML",alt_arma4$coef,include.mean = T)
arma_alt5_results = arima_fc(data_ts,length(train),c(2,0,0),c(0,0,0),"ML",alt_arma5$coef,include.mean = T)

#Model fitted by auto.arima, which is white noise, outperforms alternate models on both AIC and holdout-set RMSE.

#ARMAX

#VIX forecasting - ESM
vix_ts = ts(data$VIXCLS, start = c(1990,2), frequency = 12)
vix_pct = as.numeric(diff(log(vix_ts)))
vix_pct_ts = ts(vix_pct, start = c(1990,2), frequency = 12)

vix_train = window(vix_pct_ts, start = c(1990,2), end = c(2015,4))
vix_holdout = window(vix_pct_ts, start = c(2015,5))

vix_model =  HoltWinters(vix_train, beta = F, gamma = F)
vix_results = esm_fc(vix_train, vix_holdout, vix_model$alpha, vix_model$coefficients[1])

vix_fc = c(vix_model$fitted[,"xhat"],vix_results$fc)

#Copper forecasting - ARMA

copper_ts = ts(data$Copper, start = c(1990,2), frequency = 12)
copper_pct = as.numeric(diff(log(copper_ts)))
copper_pct_ts = ts(copper_pct, start = c(1990,2), frequency = 12)

copper_train = window(copper_pct_ts, start = c(1990,2), end = c(2015,4))
copper_holdout = window(copper_pct_ts, start = c(2015,5))

copper_model =  auto.arima(copper_train, stationary = T, seasonal = F)
copper_results = arima_fc(copper_pct_ts,length(copper_train),c(1,0,0),c(0,0,0),"ML",copper_model$coef,include.mean = F)

copper_fc = c(copper_model[["fitted"]], copper_results$fc)

#ARMAX set up
ntotal = nrow(data)

df = data.frame(vix_fc = vix_fc[1:(ntotal-2)],
                copper_fc = copper_fc[1:(ntotal-2)],
                sp500_ret = data$sp500_ret[2:(ntotal-1)]
                )

ntrain = floor(nrow(df) * 0.7)

reg = lm(sp500_ret ~ vix_fc + copper_fc, data = df[1:ntrain,])
print(summary(reg))

acf(reg$residuals)
pacf(reg$residuals)

#Residuals look like white noise so no further model needs to be fitted for the residuals.

pred_reg = predict(reg, df[(ntrain+1):nrow(df),])

rmse_reg = sqrt(mean((df$sp500_ret[(ntrain+1):nrow(df)] - pred_reg)^2))

