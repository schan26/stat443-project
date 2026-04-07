train_end   <- c(2015, 4)
holdout_start <- c(2015, 5)

data <- read.csv("data/processed/cleanedData.csv")

# FIX: Parse "YYYY-MM" dates by appending "-01" so as.Date() can read them
data$observation_date <- as.Date(paste0(data$observation_date, "-01"))

data_ts <- ts(data$sp500_ret, start = c(1990,2), frequency = 12)

train   <- window(data_ts, start = c(1990,2), end = train_end)
holdout <- window(data_ts, start = holdout_start)

# Subset data to training period only
train_data <- data[1:length(train),]

vars <- c("CPI", "VIXCLS", "INDPRO", "Gold", "Silver")

# ── times series plot analysis ────────────────────────────────────────────────
plot(dates, na.omit(train_data[["sp500_ret"]]) ,
     type = "l", lwd = 1.5,
     xlab = "Date", ylab = "sp500_ret",
     main = "Time Series — sp500_ret", cex.main = 2.5)

abline(h = mean(na.omit(train_data[["sp500_ret"]])), col = "red", lty = 2, lwd = 1)


par(mfrow = c(2, 2))

for (v in vars) {
  # FIX: Use train_data instead of train (since train is a ts object)
  x     <- na.omit(train_data[[v]]) 
  dates <- train_data$observation_date[!is.na(train_data[[v]])]
  
  plot(dates, x,
       type = "l", lwd = 1.5,
       xlab = "Date", ylab = v,
       main = paste("Time Series —", v), cex.main = 2)
  
  abline(h = mean(x), col = "red", lty = 2, lwd = 1)
}

par(mfrow = c(1, 1))


# ── acf/pacf analysis ────────────────────────────────────────────────
vars <- c("sp500_ret", vars)
plot_acf_pacf <- function(var_name, df, lag.max = 36) {
  x <- na.omit(df[[var_name]])
  par(mfrow = c(1, 2))
  
  acf(x,  lag.max = lag.max,
      main = paste("ACF —", var_name), cex.main = 3)
  
  pacf(x, lag.max = lag.max,
       main = paste("PACF —", var_name), cex.main = 3)
}

for (v in vars) {
  plot_acf_pacf(v, train_data)
}


# ── Variogram function ────────────────────────────────────────────────
variogram <- function(y, lagmax = 10, iprint = FALSE) {
  G  <- rep(1, lagmax)
  n  <- length(y)
  if (lagmax > n) lagmax <- n - 2
  y1 <- y[-1]; y2 <- y[-n]
  d1 <- y1 - y2; denom <- var(d1)
  for (k in 2:lagmax) {
    y1 <- y[(k+1):n]; y2 <- y[1:(n-k)]
    dk <- y1 - y2
    G[k] <- var(dk) / denom
  }
  ac <- c(acf(y, plot = FALSE, lag.max = lagmax)$acf)
  H  <- (1 - ac[-1]) / (1 - ac[2])
  if (iprint) print(cbind(G, H))
  list(G = G, H = H)
}



# ── variogram analysis ────────────────────────────────────────────────
par(mfrow = c(2, 2))

for (v in vars) {
  x    <- na.omit(train_data[[v]])
  vg   <- variogram(x, lagmax = 12)
  lags <- 1:12
  
  matplot(lags, cbind(vg$G, vg$H),
          type = "b", pch = c(16, 17), 
          col = c("steelblue", "darkorange"), lwd = 2.5,
          xlab = "Lag (months)", ylab = "variogram",
          main = paste("Variogram —", v), cex.main = 2.5)
  
  abline(h = 1, col = "black", lty = 2, lwd = 1.5)
  
  legend("bottomright",
         legend = c("G", "H"),
         col = c("steelblue", "darkorange"),
         pch = c(16, 17), lwd = 2, bty = "n", 
         cex = 3)
}

par(mfrow = c(1, 1))


# ── ccf analysis ────────────────────────────────────────────────
indpro_diff = diff(data$INDPRO)
copper_ldiff = diff(log(data$Copper))
silver_ldiff = diff(log(data$Silver))
gold_ldiff = diff(log(data$Gold))

ccf(indpro_diff[1:ntrain], c(train), main = "CCF of diff(INDPRO) & SP500 returns",  ylab = "CCF", lag.max = 10)
ccf(data$CPI[1:ntrain], c(train), main = "CCF of CPI % change & SP500 returns",  ylab = "CCF", lag.max = 10)
ccf(copper_ldiff[1:ntrain], c(train), main = "CCF of log.diff(Copper) & SP500 returns",  ylab = "CCF", lag.max = 10)
ccf(silver_ldiff[1:ntrain], c(train), main = "CCF of log.diff(Silver) & SP500 returns",  ylab = "CCF", lag.max = 10)
ccf(gold_ldiff[1:ntrain], c(train), main = "CCF of log.diff(Gold) & SP500 returns",  ylab = "CCF", lag.max = 10)
