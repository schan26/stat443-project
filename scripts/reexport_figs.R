### -------------------
### Run ALL other scripts before running this script
### -------------------


# ── Common style settings ──────────────────────────────────────────────────────
main_cex <- 2.0
lab_cex  <- 1.5
axis_cex <- 1.2

# ── Variogram function ─────────────────────────────────────────────────────────
variogram <- function(y, lagmax = 12, iprint = FALSE) {
  G <- rep(1, lagmax)
  n <- length(y)
  if (lagmax > n) lagmax <- n - 2
  
  y1 <- y[-1]
  y2 <- y[-n]
  d1 <- y1 - y2
  denom <- var(d1)
  
  for (k in 2:lagmax) {
    y1 <- y[(k + 1):n]
    y2 <- y[1:(n - k)]
    dk <- y1 - y2
    G[k] <- var(dk) / denom
  }
  
  ac <- c(acf(y, plot = FALSE, lag.max = lagmax)$acf)
  H <- (1 - ac[-1]) / (1 - ac[2])
  
  if (iprint) print(cbind(G, H))
  list(G = G, H = H)
}



# ── Load and prepare data ──────────────────────────────────────────────────────
data <- read.csv("data/processed/cleanedData.csv")
data$observation_date <- as.Date(paste0(data$observation_date, "-01"))

train_end <- c(2015, 4)

data_ts <- ts(data$sp500_ret, start = c(1990, 2), frequency = 12)
train <- window(data_ts, start = c(1990, 2), end = train_end)
train_data <- data[1:length(train), ]

sp500_x <- na.omit(train_data$sp500_ret)
sp500_dates <- train_data$observation_date[!is.na(train_data$sp500_ret)]



# ── Time series plot ───────────────────────────────────────────────────────────
if (!dir.exists("figs/ts")) {
  dir.create("figs/ts")
}

png("figs/ts/sp500_tsplot.png", width = 2000, height = 1300, res = 220)


par(
  mar = c(5.5, 5.5, 4.5, 2) + 0.1,
  cex.main = main_cex,
  cex.lab = lab_cex,
  cex.axis = axis_cex
)

plot(
  sp500_dates, sp500_x,
  type = "l",
  lwd = 1.5,
  xlab = "Date",
  ylab = "Monthly return (%)",
  main = "Time Series — sp500 returns"
)

abline(h = mean(sp500_x), col = "red", lty = 2, lwd = 1.2)

dev.off()

# ── ACF / PACF plot ────────────────────────────────────────────────────────────

if (!dir.exists("figs/acf_pacf")) {
  dir.create("figs/acf_pacf")
}

png("figs/acf_pacf/sp500_ret_acfpacf.png", width = 2400, height = 1200, res = 220)
par(
  mfrow = c(1, 2),
  mar = c(5.5, 5.5, 4.5, 2) + 0.1,
  cex.main = 1.5,
  cex.lab = lab_cex,
  cex.axis = axis_cex
)

acf(
  sp500_x,
  lag.max = 20,
  lwd = 2,
  xlab = "Lag (Months)",
  ylab = "ACF",
  main = "ACF — sp500 returns"
)

pacf(
  sp500_x,
  lag.max = 20,
  lwd = 2,
  xlab = "Lag (Months)",
  ylab = "Partial ACF",
  main = "PACF — sp500 returns"
)

dev.off()

png("figs/acf_pacf/BAA_acfpacf.png", width = 2400, height = 1200, res = 220)
par(
  mfrow = c(1, 2),
  mar = c(5.5, 5.5, 4.5, 2) + 0.1,
  cex.main = 1.5,
  cex.lab = lab_cex,
  cex.axis = axis_cex
)

acf(BAA_train, 
    lag.max = 20,
    lwd = 2,
    xlab = "Lag (Months)",
    ylab = "ACF",
    main = "ACF of differenced Baa spread")
pacf(BAA_train, 
     lag.max = 20,
     lwd = 2,
     xlab = "Lag (Months)",
     ylab = " Partial ACF",
     main = "PACF of differenced Baa spread")
dev.off()

png("figs/acf_pacf/VIX_acfpacf.png", width = 2400, height = 1200, res = 220)
par(
  mfrow = c(1, 2),
  mar = c(5.5, 5.5, 4.5, 2) + 0.1,
  cex.main = 1.5,
  cex.lab = lab_cex,
  cex.axis = axis_cex
)

acf(vix_train, 
    main = "ACF of differenced VIX",
    lag.max = 20,
    lwd = 2,
    xlab = "Lag (Months)",
    ylab = "ACF")
pacf(vix_train, 
     main = "PACF of differenced VIX",
     lag.max = 20,
     lwd = 2,
     xlab = "Lag (Months)",
     ylab = "Partial ACF")
dev.off()


# ── Variogram plot of S&P500 returns───────────────────────────────────────────
if (!dir.exists("figs/vario")) {
  dir.create("figs/vario")
}

vg <- variogram(sp500_x, lagmax = 12)
lags <- 1:12

png("figs/vario/sp500_ret_vario.png", width = 2200, height = 1300, res = 220)

par(
  mar = c(5.5, 5.5, 4.5, 2) + 0.1,
  cex.main = main_cex,
  cex.lab = lab_cex,
  cex.axis = axis_cex
)

matplot(
  lags, cbind(vg$G, vg$H),
  type = "b",
  pch = c(16, 17),
  col = c("steelblue", "darkorange"),
  lwd = 2.5,
  xlab = "Lag (months)",
  ylab = "Variogram",
  main = "Variogram — sp500 returns"
)

abline(h = 1, col = "black", lty = 2, lwd = 1.2)

legend(
  "bottomright",
  legend = c("G", "H"),
  col = c("steelblue", "darkorange"),
  pch = c(16, 17),
  lwd = 2,
  bty = "n",
  cex = 1.3
)

dev.off()



# ── ccf plots ─────────────────────────────────────────────────────────────
if (!dir.exists("figs/ccf")) {
  dir.create("figs/ccf")
}

png("figs/ccf/Copper_ccf.png", width = 1497, height = 1044, res = 220)
ccf(copper_ldiff[1:ntrain], c(train), main = "CCF of log differenced Copper & SP500 returns",  ylab = "CCF", xlab="Lag (Months)", lag.max = 10)
dev.off()

png("figs/ccf/INDPRO_ccf.png", width = 1497, height = 1044, res = 220)
ccf(indpro_diff[1:ntrain], c(train), main = "CCF of differenced INDPRO & SP500 returns",  ylab = "CCF", xlab="Lag (Months)", lag.max = 10)
dev.off()

png("figs/ccf/CPI_ccf.png", width = 1497, height = 1044, res = 220)
ccf(data$CPI[1:ntrain], c(train), main = "CCF of CPI % change & SP500 returns",  ylab = "CCF", xlab="Lag (Months)", lag.max = 10)
dev.off()

png("figs/ccf/Silver_ccf.png", width = 1497, height = 1044, res = 220)
ccf(silver_ldiff[1:ntrain], c(train), main = "CCF of log differenced Silver & SP500 returns",  ylab = "CCF", xlab="Lag (Months)", lag.max = 10)
dev.off()

png("figs/ccf/Gold_ccf.png", width = 1497, height = 1044, res = 220)
ccf(gold_ldiff[1:ntrain], c(train), main = "CCF of log differenced Gold & SP500 returns",  ylab = "CCF", xlab="Lag (Months)", lag.max = 10)
dev.off()

png("figs/ccf/BAA_ccf.png", width = 1497, height = 1044, res = 220)
ccf(BAA_train,train, main = "CCF of differenced Baa spread & SP500 returns", ylab = "CCF", xlab="Lag (Months)", lag.max = 10)
dev.off()

png("figs/ccf/VIX_ccf.png", width = 1497, height = 1044, res = 220)
ccf(vix_train,train, main = "CCF of differenced VIX & SP500 returns", ylab = "CCF", xlab="Lag (Months)", lag.max = 10) 
dev.off()


# ───────────────────────────────────────────────────────────────
cat("Done. New figures saved in figs/\n")