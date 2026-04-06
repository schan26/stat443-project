# Re-export cleaner S&P 500 report figures
# Save this as: scripts/reexport_figs.R

data <- read.csv("data/processed/cleanedData.csv")
data$observation_date <- as.Date(paste0(data$observation_date, "-01"))

train_end <- c(2015, 4)

data_ts <- ts(data$sp500_ret, start = c(1990, 2), frequency = 12)
train <- window(data_ts, start = c(1990, 2), end = train_end)
train_data <- data[1:length(train), ]

if (!dir.exists("figs/report")) {
  dir.create("figs/report", recursive = TRUE)
}

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

# -----------------------------
# Common style settings
# -----------------------------
main_cex <- 2.0
lab_cex  <- 1.5
axis_cex <- 1.2

sp500_x <- na.omit(train_data$sp500_ret)
sp500_dates <- train_data$observation_date[!is.na(train_data$sp500_ret)]

# -----------------------------
# 1) Time series plot
# -----------------------------
png("figs/report/sp500_tsplot.png", width = 2000, height = 1300, res = 220)

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
  main = "Time Series — sp500_ret"
)

abline(h = mean(sp500_x), col = "red", lty = 2, lwd = 1.2)

dev.off()

# -----------------------------
# 2) ACF / PACF plot
# -----------------------------
png("figs/report/sp500_ret_acfpacf.png", width = 2400, height = 1200, res = 220)

par(
  mfrow = c(1, 2),
  mar = c(5.5, 5.5, 4.5, 2) + 0.1,
  cex.main = 1.8,
  cex.lab = lab_cex,
  cex.axis = axis_cex
)

acf(
  sp500_x,
  lag.max = 36,
  lwd = 2,
  xlab = "Lag",
  ylab = "ACF",
  main = "ACF — sp500_ret"
)

pacf(
  sp500_x,
  lag.max = 36,
  lwd = 2,
  xlab = "Lag",
  ylab = "Partial ACF",
  main = "PACF — sp500_ret"
)

dev.off()

# -----------------------------
# 3) Variogram plot
# -----------------------------
vg <- variogram(sp500_x, lagmax = 12)
lags <- 1:12

png("figs/report/sp500_ret_vario.png", width = 2200, height = 1300, res = 220)

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
  main = "Variogram — sp500_ret"
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

cat("Done. New figures saved in figs/report/\n")