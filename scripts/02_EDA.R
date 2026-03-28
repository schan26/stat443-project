# EDA code 

data <- read.csv("data/processed/cleanedData.csv")

data_ts <- ts(merged_data$sp500_ret, start = c(1990,2), frequency = 12)
train = window(data_ts, start = c(1990,2), end = c(2015,4)) #train data from 1990-02 to 2015-04
holdout = window(data_ts, start = c(2015,5)) #holdout data from 2015-05 onwards
