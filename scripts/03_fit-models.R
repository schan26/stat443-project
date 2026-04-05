library(forecast)

# ── Fit models ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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

data_start = c(1990,2)

data <- read.csv("data/processed/cleanedData.csv")
data_ts <- ts(data$sp500_ret, start = data_start, frequency = 12)

ntrain = ceiling(nrow(data)*0.7)
nholdout = nrow(data) - ntrain
train_end = c(2015,4)
holdout_start = c(2015,5)

train = window(data_ts, start = data_start, end = train_end) #train data from 1990-02 to 2015-04
holdout = window(data_ts, start = holdout_start) #holdout data from 2015-05 onwards

#────── Baseline methods - Persistence and IID ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
persist_results = persist_fc(train, holdout)
iid_results = iid_fc(train, holdout)

#────── Exponential smoothing methods - Holt linear and simple exponential smoothing ──────────────────────────────────────────────────────────────────────────────
esm_model = HoltWinters(train, beta = F, gamma = F)
esm_results = esm_fc(train, holdout, esm_model$alpha, esm_model$coefficients[1])

lholt_model = HoltWinters(train, gamma = F)
lholt_results = lholt_fc(train, holdout, lholt_model$alpha, lholt_model$beta, lholt_model$coefficients[1], lholt_model$coefficients[2])

#Simple exponential smoothing performs better than linear Holt due to the lack of a consistent trend in the data.

#────── ARMA ──────────────────────────────────────────────────────────────────────────────
arma_model = auto.arima(train, stationary = T, seasonal = F)
arma_results = arima_fc(data_ts,length(train),c(0,0,0),c(0,0,0),"ML",arma_model$coef,include.mean = T)

#Vary P and Q a little and compare AIC and holdout set RMSE
alt_arma101 = arima(train,order = c(1,0,1))
alt_arma001 = arima(train,order = c(0,0,1))
alt_arma002 = arima(train,order = c(0,0,2))
alt_arma100 = arima(train,order = c(1,0,0))
alt_arma200 = arima(train,order = c(2,0,0))

arma_alt101_results = arima_fc(data_ts,length(train),c(1,0,1),c(0,0,0),"ML",alt_arma101$coef,include.mean = T)
arma_alt001_results = arima_fc(data_ts,length(train),c(0,0,1),c(0,0,0),"ML",alt_arma001$coef,include.mean = T)
arma_alt002_results = arima_fc(data_ts,length(train),c(0,0,2),c(0,0,0),"ML",alt_arma002$coef,include.mean = T)
arma_alt100_results = arima_fc(data_ts,length(train),c(1,0,0),c(0,0,0),"ML",alt_arma100$coef,include.mean = T)
arma_alt200_results = arima_fc(data_ts,length(train),c(2,0,0),c(0,0,0),"ML",alt_arma200$coef,include.mean = T)

#Model fitted by auto.arima, which is white noise, outperforms alternate models on both AIC and holdout-set RMSE.

#────── ARMAX ──────────────────────────────────────────────────────────────────────────────

#BAA forecasting - ARMA
BAA_ts = ts(diff(data$BAA10YM), start = data_start, frequency = 12)
BAA_train = window(BAA_ts, start = data_start, end = train_end)
BAA_holdout = window(BAA_ts, start = holdout_start)

ccf(BAA_train,train, main = "CCF of differenced BAA spread & SP500 returns", ylab = "CCF")

acf(BAA_train, main = "ACF of differenced BAA spread")
pacf(BAA_train, main = "PACF of differenced BAA spread")

BAA_model = auto.arima(BAA_train, stationary = T, seasonal = F)
BAA_results = arima_fc(BAA_ts,length(BAA_train),c(2,0,1),c(0,0,0),"ML",BAA_model$coef,include.mean = F)

BAA_fc <- c(BAA_train,BAA_results$fc)

#VIX forecasting - ARMA
vix_ts = ts(data$VIXCLS, start = data_start, frequency = 12)
vix_diff = as.numeric(diff(vix_ts))
vix_diff_ts = ts(vix_diff, start = data_start, frequency = 12)

vix_train = window(vix_diff_ts, start = data_start, end = train_end)
vix_holdout = window(vix_diff_ts, start = holdout_start)

acf(vix_train, main = "ACF of differenced VIX" )
pacf(vix_train, main = "PACF of differenced VIX")

ccf(vix_train,train, main = "CCF of differenced VIX & SP500 returns")

vix_model = auto.arima(vix_train, stationary = T, seasonal = F)
vix_results = arima_fc(vix_diff_ts,length(vix_train),c(0,0,3),c(0,0,0),"ML",vix_model$coef,include.mean = F)

vix_fc = c(vix_train,vix_results$fc)

#ARMAX set up
ntotal = nrow(data) #This equals 432

df = data.frame(vix_fc = vix_fc[1:(ntotal-1)],
                sp500_ret = data$sp500_ret[1:(ntotal-1)],
                BAA_fc = BAA_fc[1:(ntotal-1)])

reg = lm(sp500_ret ~ vix_fc + BAA_fc, data = df[1:ntrain,])
print(summary(reg))

acf(reg$residuals, main = "ACF of regression residuals")
pacf(reg$residuals, main = "PACF of regression residuals")

#Residuals look like white noise so no further model needs to be fitted for the residuals.

pred_reg = predict(reg, df[(ntrain+1):nrow(df),])

rmse_reg = sqrt(mean((df$sp500_ret[(ntrain+1):nrow(df)] - pred_reg)^2))

#Alternate ARMAX
indpro_ts = ts(diff(data$INDPRO), start = data_start, frequency = 12)

indpro_train = window(indpro_ts, start = data_start, end = train_end)
indpro_holdout = window(indpro_ts, start = holdout_start)

ccf(indpro_train, train, main = "CCF of differenced INDPRO & SP500 returns")

df2 = data.frame(indpro_l2 = c(indpro_ts)[1:(ntotal-3)],
                 vix_fc = vix_fc[3:(ntotal-1)],
                 sp500_ret = data$sp500_ret[3:(ntotal-1)],
                 BAA_fc = BAA_fc[3:(ntotal-1)])

reg_alt = lm(sp500_ret ~ vix_fc + BAA_fc + indpro_l2, data = df2[1:(ntrain-2),])
print(summary(reg_alt))

acf(reg_alt$residuals, main = "ACF of regression residuals")
pacf(reg_alt$residuals, main = "PACF of regression residuals")

pred_reg_alt = predict(reg_alt, df2[(ntrain-1):nrow(df2),])
rmse_reg_alt = sqrt(mean((df2$sp500_ret[(ntrain-1):nrow(df2)] - pred_reg_alt)^2))

#Alternate ARMAX 2
copper_ts = ts(diff(log(data$Copper)), start = data_start, frequency = 12)

copper_train = window(copper_ts, start = data_start, end = train_end)
copper_holdout = window(copper_ts, start = holdout_start)

ccf(copper_train, train, main = "CCF of log differenced Copper & SP500 returns")

df3 = data.frame(copper_l2 = c(copper_ts)[1:(ntotal-3)],
                 vix_fc = vix_fc[3:(ntotal-1)],
                 sp500_ret = data$sp500_ret[3:(ntotal-1)],
                 BAA_fc = BAA_fc[3:(ntotal-1)])

reg_alt2 = lm(sp500_ret ~ vix_fc + BAA_fc + copper_l2, data = df3[1:(ntrain-2),])
print(summary(reg_alt2))

acf(reg_alt2$residuals, main = "ACF of regression residuals")
pacf(reg_alt2$residuals, main = "PACF of regression residuals")

pred_reg_alt2 = predict(reg_alt2, df3[(ntrain-1):nrow(df3),])
rmse_reg_alt2 = sqrt(mean((df3$sp500_ret[(ntrain-1):nrow(df3)] - pred_reg_alt2)^2))

#Initially, the focus was on finding leading predictors but most predictors had either a weak leading effect or no leading effect at all over the S&P 500, which led to a worser model. 
#Hence, the focus shifted to finding variables that had a strong lag 0 correlation with the SP500 returns, and use the 1-step forecasts of these variables in the holdout set to build a good ARMAX model. 
#Two such variables were VIX and Baa Corporate bond yield relative to yield on 10-Year Treasury Constant Maturity, which had moderate lag 0 correlations with the SP500 returns over the training set.

rmse_table = cbind(iid_results$rmse, persist_results$rmse, esm_results$rmse, lholt_results$rmse, arma_results$rmse, rmse_reg)
colnames(rmse_table) = c("iid_rmse","persist_rmse","esm_rmse","lholt_rmse","arma_rmse","armax_rmse")
rownames(rmse_table) = NA

fc_table = cbind(data$sp500_ret[(ntrain+1):(ntotal-1)],iid_results$fc[1:(nholdout-1)],persist_results$fc[1:(nholdout-1)],esm_results$fc[1:(nholdout-1)],lholt_results$fc[1:(nholdout-1)],arma_results$fc[1:(nholdout-1)],c(pred_reg))
colnames(fc_table) = c("holdout","iid_fc","persist_fc","esm_fc","lholt_fc","arma_fc","armax_fc")
rownames(fc_table) = data$observation_date[(ntrain+1):(ntotal-1)]

#pairs(cbind(iid_results$fc,esm_results$fc))

