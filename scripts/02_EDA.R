library(forecast)

data <- read.csv("data/processed/cleanedData.csv")

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

