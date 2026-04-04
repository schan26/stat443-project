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

vars <- c("CPI", "VIXCLS", "T10Y2Y", "INDPRO", "Gold", "Silver")

# ── Time series plots ──────────────────────────────────────────────────────────
if (!dir.exists("figs/ts")) {
  dir.create("figs/ts")
}

png("figs/ts/sp500_tsplot.png", width = 800, height = 600)
plot(dates, na.omit(train_data[["sp500_ret"]]) ,
     type = "l", lwd = 1.5,
     xlab = "Date", ylab = "sp500_ret",
     main = "Time Series — sp500_ret", cex.main = 2.5)

abline(h = mean(na.omit(train_data[["sp500_ret"]])), col = "red", lty = 2, lwd = 1)
dev.off()


png("figs/ts/tsplots_%03d.png", width = 1000, height = 800)
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
dev.off()

par(mfrow = c(1, 1))

# ── ACF / PACF ────────────────────────────────────────────────────────────────
if (!dir.exists("figs/acf_pacf")) {
  dir.create("figs/acf_pacf")
}


vars <- c("sp500_ret", vars)
plot_acf_pacf <- function(var_name, df, lag.max = 36) {
  png(paste0("figs/acf_pacf/", v, "_acfpacf.png"), width = 1200, height = 650)
  x <- na.omit(df[[var_name]])
  par(mfrow = c(1, 2))
  
  acf(x,  lag.max = lag.max,
      main = paste("ACF —", var_name), cex.main = 3)
  
  pacf(x, lag.max = lag.max,
       main = paste("PACF —", var_name), cex.main = 3)
  dev.off()
}

for (v in vars) {
  plot_acf_pacf(v, train_data)
}

# ── Variogram function ────────────────────────────────────────────────────────
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

# ── Variograms (levels) ───────────────────────────────────────────────────────
if (!dir.exists("figs/vario")) {
  dir.create("figs/vario")
}

png("figs/vario/sp500_ret_vario.png", width = 1400, height = 800)
x    <- na.omit(train_data[["sp500_ret"]])
vg   <- variogram(x, lagmax = 12)
lags <- 1:12

matplot(lags, cbind(vg$G, vg$H),
        type = "b", pch = c(16, 17), 
        col = c("steelblue", "darkorange"), lwd = 2.5,
        xlab = "Lag (months)", ylab = "variogram",
        main = "Variogram — sp500_ret", cex.main = 2.5)

abline(h = 1, col = "black", lty = 2, lwd = 1.5)

legend("bottomright",
       legend = c("G", "H"),
       col = c("steelblue", "darkorange"),
       pch = c(16, 17), lwd = 2, bty = "n",
       cex = 3)
dev.off()


png(paste0("figs/vario/vario_%03d.png"), width = 1400, height = 800)
par(mfrow = c(2, 2))

vars <- vars[vars != "sp500_ret"]
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
dev.off()

par(mfrow = c(1, 1))

# ── Stationarity classification ───────────────────────────────────────────────
#stationary_vars     <- c("CPI", "VIXCLS", "T10Y2Y")
#non_stationary_vars <- c("Gold", "Silver", "INDPRO")
#predictors          <- c(stationary_vars, non_stationary_vars)

# ── Differenced series plots ──────────────────────────────────────────────────
#par(mfrow = c(2, 2))

#for (v in non_stationary_vars) {
#  sub <- train_data[!is.na(train_data[[v]]), c("observation_date", v)]
#  sub <- sub[order(sub$observation_date), ]
  
#  x          <- sub[[v]]
#  dx         <- diff(x)
#  dates_diff <- sub$observation_date[-1]
  
#   plot(dates_diff, dx,
#        type = "l", col = "steelblue", lwd = 1.5,
#        xlab = "Date", ylab = paste("diff(", v, ")"),
#        main = paste("Differenced —", v))
#   
#   abline(h = 0, col = "firebrick", lty = 2, lwd = 1)
# }
# 
# par(mfrow = c(1, 1))

# ── Variograms (differenced) ──────────────────────────────────────────────────
# par(mfrow = c(2, 2))
# 
# for (v in non_stationary_vars) {
#   x    <- na.omit(train_data[[v]])
#   dx   <- diff(x)
#   vg   <- variogram(dx, lagmax = 36)
#   lags <- 1:36
#   
#   matplot(lags, cbind(vg$G, vg$H),
#           type = "b", pch = c(16, 17),
#           col = c("steelblue", "darkorange"), lwd = 2,
#           xlab = "Lag (months)", ylab = expression(gamma(h)),
#           main = paste("Variogram diff —", v))
#   
#   abline(h = 1, col = "firebrick", lty = 2, lwd = 1.5)
#   
#   legend("bottomright",
#          legend = c("G (empirical)", "H (stationary)"),
#          col = c("steelblue", "darkorange"),
#          pch = c(16, 17), lwd = 2, bty = "n")
# }
# 
# par(mfrow = c(1, 1))

# ── CCF vs sp500_ret ──────────────────────────────────────────────────────────
sp500_train <- na.omit(train_data$sp500_ret)

if (!dir.exists("figs/ccf")) {
  dir.create("figs/ccf")
}

par(mfrow = c(1, 1),cex.main = 5)
predictors <- c("CPI", "INDPRO")
for (v in predictors) {
  png(paste0("figs/ccf/", v, "_ccf.png"), width = 800, height = 600)
  x <- na.omit(train_data[[v]])
  
    if (v == "INDPRO") {
      x <- diff(x)
      # FIX: Drop first element to align lengths instead of double differencing returns
      sp500_aligned <- sp500_train[-1]
      label <- paste0("diff(", v, ")")
    } else {
      sp500_aligned <- sp500_train
      label <- v
    }
  
    min_len <- min(length(x), length(sp500_aligned))
    x <- tail(x, min_len)
    sp500_aligned <- tail(sp500_aligned, min_len)
  
    ccf(x, sp500_aligned,
        lag.max = 24,
        main = paste("CCF —", label, "vs sp500_ret"), 
        ylab = "CCF", xlab = "Lag (months)")
    dev.off()
}

# for (v in predictors) {
#   x <- na.omit(train_data[[v]])
#   
#   if (v %in% non_stationary_vars) {
#     x             <- diff(x)
#     # FIX: Drop first element to align lengths instead of double differencing returns
#     sp500_aligned <- sp500_train[-1]
#     label         <- paste0("diff(", v, ")")
#   } else {
#     sp500_aligned <- sp500_train
#     label         <- v
#   }
#   
#   min_len       <- min(length(x), length(sp500_aligned))
#   x             <- tail(x, min_len)
#   sp500_aligned <- tail(sp500_aligned, min_len)
#   
#   ccf(x, sp500_aligned,
#       lag.max = 24, col = "steelblue", lwd = 2, ci.col = "firebrick",
#       main = paste("CCF —", label, "vs sp500_ret"),
#       ylab = "CCF", xlab = "Lag (months)")
# }
# 
# par(mfrow = c(1, 1))

# # ── CCF log-diff(VIXCLS) vs SP500 ────────────────────────────────────────────
# d1 <- na.omit(train_data[, c("observation_date", "VIXCLS", "sp500_ret")])
# d1 <- d1[order(d1$observation_date), ]
# 
# vix_ld <- diff(log(d1$VIXCLS))
# # FIX: Align index without double-differencing S&P 500 returns
# sp1    <- d1$sp500_ret[-1] 
# min1   <- min(length(vix_ld), length(sp1))
# 
# ccf(tail(vix_ld, min1), tail(sp1, min1),
#     lag.max = 24,
#     main = "CCF — log-diff(VIXCLS) vs SP500",
#     ylab = "CCF", xlab = "Lag (months)",
#     col = "pink", lwd = 2, ci.col = "firebrick")
# 
# # ── CCF log-diff(INDPRO) vs SP500 ────────────────────────────────────────────
# d2 <- na.omit(train_data[, c("observation_date", "INDPRO", "sp500_ret")])
# d2 <- d2[order(d2$observation_date), ]
# 
# indpro_ld <- diff(log(d2$INDPRO))
# # FIX: Align index without double-differencing S&P 500 returns
# sp2       <- d2$sp500_ret[-1]
# min2      <- min(length(indpro_ld), length(sp2))
# 
# ccf(tail(indpro_ld, min2), tail(sp2, min2),
#     lag.max = 24,
#     main = "CCF — log-diff(INDPRO) vs SP500",
#     ylab = "CCF", xlab = "Lag (months)",
#     col = "steelblue", lwd = 2, ci.col = "firebrick")
# 
# par(mfrow = c(1, 1))
