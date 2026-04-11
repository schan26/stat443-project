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

regression_fc = function(data,ntrain) {
  
  reg = lm(sp500_ret ~ ., data = data[1:ntrain,])
  print(summary(reg))
  
  acf(reg$residuals, main = "ACF of regression residuals")
  pacf(reg$residuals, main = "PACF of regression residuals")
  
  pred_reg = predict(reg, data[(ntrain+1):nrow(data),])
  
  rmse_reg = sqrt(mean((data$sp500_ret[(ntrain+1):nrow(data)] - pred_reg)^2))
  
  return(list(rmse = rmse_reg, fc = pred_reg, resid = reg$residuals, model = reg))
}

data_start = c(1990,2)

data <- read.csv("data/processed/cleanedData.csv")
data_ts <- ts(data$sp500_ret, start = data_start, frequency = 12)

ntotal = nrow(data) # this equals 432
ntrain = ceiling(nrow(data)*0.7)
nholdout = nrow(data) - ntrain

train_end = c(2015,4)
holdout_start = c(2015,5)

train = window(data_ts, start = data_start, end = train_end) # train data from 1990-02 to 2015-04
holdout = window(data_ts, start = holdout_start) # holdout data from 2015-05 onwards

#────── Baseline methods - Persistence and IID ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
persist_results = persist_fc(train, holdout)
iid_results = iid_fc(train, holdout)

#────── Exponential smoothing methods - Holt linear and simple exponential smoothing ──────────────────────────────────────────────────────────────────────────────
esm_model = HoltWinters(train, beta = F, gamma = F)
esm_results = esm_fc(train, holdout, esm_model$alpha, esm_model$coefficients[1])

lholt_model = HoltWinters(train, gamma = F)
lholt_results = lholt_fc(train, holdout, lholt_model$alpha, lholt_model$beta, lholt_model$coefficients[1], lholt_model$coefficients[2])

# Simple exponential smoothing performs better than linear Holt due to the lack of a consistent trend in the data.

#────── ARMA ──────────────────────────────────────────────────────────────────────────────
arma_model = auto.arima(train, stationary = T, seasonal = F)
arma_results = arima_fc(data_ts,length(train),c(0,0,0),c(0,0,0),"ML",arma_model$coef,include.mean = T)

# Vary P and Q a little and compare AIC and holdout set RMSE
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

# Model fitted by auto.arima, which is white noise, outperforms alternate models on both AIC and holdout-set RMSE.

#────── ARMAX ──────────────────────────────────────────────────────────────────────────────

# BAA forecasting - ARMA
BAA_ts = ts(diff(data$BAA10YM), start = data_start, frequency = 12)
BAA_train = window(BAA_ts, start = data_start, end = train_end)
BAA_holdout = window(BAA_ts, start = holdout_start)

plot(dates, BAA_train,
     type = "l", lwd = 1.5,
     xlab = "Date", ylab = "Differenced Baa Spread (%)",
     main = "Time Series — Differenced Baa Spread")

abline(h = mean(BAA_train), col = "red", lty = 2, lwd = 1)

BAA_ccf <- ccf(as.numeric(BAA_train),as.numeric(train), main = "CCF of differenced Baa spread & SP500 returns", ylab = "CCF", xlab = "Lags (months)", cex.main = 1.2, lag.max = 10)

BAA_model = auto.arima(BAA_train, stationary = T, seasonal = F)
BAA_results = arima_fc(BAA_ts,length(BAA_train),c(2,0,1),c(0,0,0),"ML",BAA_model$coef,include.mean = F)

BAA_fc <- c(BAA_train,BAA_results$fc)

# VIX forecasting - ARMA
vix_ts = ts(diff(data$VIXCLS), start = data_start, frequency = 12)

vix_train = window(vix_ts, start = data_start, end = train_end)
vix_holdout = window(vix_ts, start = holdout_start)

plot(dates, vix_train,
     type = "l", lwd = 1.5,
     xlab = "Date", ylab = "Differenced VIX (%)",
     main = "Time Series — Differenced VIX")

abline(h = mean(vix_train), col = "red", lty = 2, lwd = 1)

VIX_ccf <- ccf(as.numeric(vix_train),as.numeric(train), main = "CCF of differenced VIX & SP500 returns", ylab = "CCF", xlab = "Lags (months)", cex.main = 1.2, lag.max = 10)

vix_model = auto.arima(vix_train, stationary = T, seasonal = F)
vix_results = arima_fc(vix_ts,length(vix_train),c(0,0,3),c(0,0,0),"ML",vix_model$coef,include.mean = F)

vix_fc = c(vix_train,vix_results$fc)

# ARMAX set up
indpro_diff = diff(data$INDPRO)
copper_ldiff = diff(log(data$Copper))
silver_ldiff = diff(log(data$Silver))
gold_ldiff = diff(log(data$Gold))

indpro_ccf <- ccf(indpro_diff[1:ntrain], c(train), main = "CCF of diff(INDPRO) & SP500 returns",  ylab = "CCF", lag.max = 10)
CPI_ccf <- ccf(data$CPI[1:ntrain], c(train), main = "CCF of CPI % change & SP500 returns",  ylab = "CCF", lag.max = 10)
copper_ccf <- ccf(copper_ldiff[1:ntrain], c(train), main = "CCF of log.diff(Copper) & SP500 returns",  ylab = "CCF", lag.max = 10)
silver_ccf <- ccf(silver_ldiff[1:ntrain], c(train), main = "CCF of log.diff(Silver) & SP500 returns",  ylab = "CCF", lag.max = 10)
gold_ccf <- ccf(gold_ldiff[1:ntrain], c(train), main = "CCF of log.diff(Gold) & SP500 returns",  ylab = "CCF", lag.max = 10)


df = data.frame(gold_l10 = gold_ldiff[1:(ntotal-11)],
                cpi_l4 = data$CPI[7:(ntotal-5)],
                indpro_l2 = indpro_diff[9:(ntotal-3)],
                copper_l2 = copper_ldiff[9:(ntotal-3)],
                silver_l1 = silver_ldiff[10:(ntotal-2)],
                vix_fc = vix_fc[11:(ntotal-1)],
                sp500_ret = data$sp500_ret[11:(ntotal-1)],
                BAA_fc = BAA_fc[11:(ntotal-1)])

# Variables ordered in descending order according to absolute correlation are vix_fc, BAA_fc, indpro_l2, CPI_l4, copper_l2, silver_l1 and gold_l10

# ntrain has to change so that data after 2015-04 cannot be used for training

# Forward selection of variables with holdout RMSE threshold of 0.01 - variables have to improve holdout RMSE by 1% to be added to the model

# 1 variable
reg_results = regression_fc(subset(df, select = c(sp500_ret, vix_fc)), ntrain = ntrain - 10)

# 2 variables
reg_2_results = regression_fc(subset(df, select = c(sp500_ret, vix_fc, BAA_fc)), ntrain = ntrain - 10) 

# 3 variables
reg_3_results = regression_fc(subset(df, select = c(sp500_ret, vix_fc, BAA_fc, indpro_l2)), ntrain = ntrain - 10)
reg_4_results = regression_fc(subset(df, select = c(sp500_ret, vix_fc, BAA_fc, cpi_l4)), ntrain = ntrain - 10)
reg_5_results = regression_fc(subset(df, select = c(sp500_ret, vix_fc, BAA_fc, copper_l2)), ntrain = ntrain - 10)
reg_6_results = regression_fc(subset(df, select = c(sp500_ret, vix_fc, BAA_fc, silver_l1)), ntrain = ntrain - 10)

# 4 variables
reg_7_results = regression_fc(subset(df, select = c(sp500_ret, vix_fc, BAA_fc, silver_l1, gold_l10)), ntrain = ntrain - 10)

# Residuals look like white noise so no further model needs to be fitted for the residuals.

# Initially, the focus was on finding leading predictors but most predictors had either a weak leading effect or no leading effect at all over the S&P 500, which led to a worser model. 
# Hence, the focus shifted to finding variables that had a strong lag 0 correlation with the SP500 returns, and use the 1-step forecasts of these variables in the holdout set to build a good ARMAX model. 
# Two such variables were VIX and Baa Corporate bond yield relative to yield on 10-Year Treasury Constant Maturity, which had moderate lag 0 correlations with the SP500 returns over the training set.

rmse_table = cbind(iid_results$rmse, persist_results$rmse, esm_results$rmse, lholt_results$rmse, arma_results$rmse, reg_6_results$rmse)
colnames(rmse_table) = c("iid_rmse","persist_rmse","esm_rmse","lholt_rmse","arma_rmse","armax_rmse")
rownames(rmse_table) = NA

fc_table = cbind(data$sp500_ret[(ntrain+1):(ntotal-1)],iid_results$fc[1:(nholdout-1)],persist_results$fc[1:(nholdout-1)],esm_results$fc[1:(nholdout-1)],lholt_results$fc[1:(nholdout-1)],arma_results$fc[1:(nholdout-1)],reg_6_results$fc)
colnames(fc_table) = c("holdout","iid_fc","persist_fc","esm_fc","lholt_fc","arma_fc","armax_fc")
rownames(fc_table) = data$observation_date[(ntrain+1):(ntotal-1)]
