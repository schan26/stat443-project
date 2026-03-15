
# load raw data
indpro <- read.csv("data/raw/INDPRO.csv")
yield_curve <- read.csv("data/raw/T10Y2Y.csv")
cpi <- read.csv("data/raw/CPI.csv", skip = 11)
vix <- read.csv("data/raw/VIXCLS.csv")
sp500 <- read.csv("data/raw/S&P 500.csv")


# change data type of observation_date 
indpro$observation_date <- as.Date(indpro$observation_date)
yield_curve$observation_date <- as.Date(yield_curve$observation_date)
vix$observation_date <- as.Date(vix$observation_date)
sp500$Date <- as.Date(sp500$Date)

tail(cpi)
tail(indpro)
tail(yield_curve)
tail(vix)
tail(sp500)
