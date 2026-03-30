library(forecast)

data <- read.csv("data/processed/cleanedData.csv")

train = window(data_ts, start = c(1990,2), end = c(2015,4)) #train data from 1990-02 to 2015-04
holdout = window(data_ts, start = c(2015,5)) #holdout data from 2015-05 onwards
data$observation_date <- as.Date(data$observation_date)

vars <- c("CPI", "VIXCLS", "T10Y2Y", "INDPRO", "sp500_ret", "Gold", "Silver")

plot_acf_pacf <- function(var_name, data, lag.max = 36) {
  x <- na.omit(data[[var_name]])
  
  par(mfrow = c(1, 2))
  
  acf(x,
      lag.max = lag.max,
      main    = paste("ACF —", var_name),
      col     = "steelblue",
      lwd     = 2,
      ci.col  = "firebrick")
  
  pacf(x,
       lag.max = lag.max,
       main    = paste("PACF —", var_name),
       col     = "steelblue",
       lwd     = 2,
       ci.col  = "firebrick")
}

for (v in vars) {
  plot_acf_pacf(v, data)
}

par(mfrow = c(1, 1))

