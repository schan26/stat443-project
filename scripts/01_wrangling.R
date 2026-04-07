library(tidyr)
library(dplyr)

# ── Helper function ────────────────────────────────────────────────────────────
fill_missing_values <- function(df, date_col = "observation_date") {
  
  all_dates <- data.frame(observation_date = seq(min(df[[date_col]]), max(df[[date_col]]), by = "day"))
  
  df <- df %>%
    full_join(all_dates, by = "observation_date") %>%
    arrange(observation_date)
  
  df <- fill(df, -observation_date, .direction = "down")
  
  return(df)
}

# ── Load raw data ──────────────────────────────────────────────────────────────
indpro <- read.csv("data/raw/INDPRO.csv")
cpi <- read.csv("data/raw/CPI.csv", skip = 11)
vix <- read.csv("data/raw/VIXCLS.csv")
sp500 <- rename(read.csv("data/raw/S&P 500.csv"), observation_date = Date, sp500_ret = Change..)[,c("observation_date","sp500_ret")]
gold <- rename(read.csv("data/raw/Gold.csv"), observation_date = Date, Gold = Value)
silver <- rename(read.csv("data/raw/Silver.csv"), observation_date = Date, Silver = Value)
copper <- rename(read.csv("data/raw/Copper.csv"), observation_date = Date, Copper = Value)
BAA10Y <- read.csv("data/raw/BAA10YM.csv")

# ── Parse dates ────────────────────────────────────────────────────────────────
indpro$observation_date <- as.Date(indpro$observation_date)
vix$observation_date <- as.Date(vix$observation_date)
sp500$observation_date <- as.Date(sp500$observation_date)
gold$observation_date <- as.Date(gold$observation_date, format = "%m/%d/%Y")
silver$observation_date <- as.Date(silver$observation_date, format = "%m/%d/%Y")
copper$observation_date <- as.Date(copper$observation_date, format = "%m/%d/%Y")
BAA10Y$observation_date <- as.Date(BAA10Y$observation_date)


# ── Forward-fill daily series ──────────────────────────────────────────────────
vix <- fill_missing_values(vix)
copper <- fill_missing_values(copper)

# ── Reshape CPI data ───────────────────────────────────────────────────────────
cpi <- pivot_longer(cpi, cols = Jan:Dec, names_to = "month", values_to = "CPI") #change from wide to long format
cpi$observation_date <- paste0(cpi$Year,"-",cpi$month, "-01") #assign default date to be the first of every month
cpi <- cpi[,c("observation_date","CPI")]
cpi$observation_date = as.Date(cpi$observation_date, format = "%Y-%b-%d")

# ── Clean S&P 500 returns ──────────────────────────────────────────────────────
sp500$sp500_ret = as.numeric(sub("%","",sp500$sp500_ret))

# ── Merge data sources ─────────────────────────────────────────────────────────
merged_data <- cpi %>%
  inner_join(vix, by = "observation_date") %>%
  inner_join(indpro, by = "observation_date") %>%
  inner_join(sp500, by = "observation_date") %>%
  inner_join(gold, by = "observation_date") %>%
  inner_join(silver, by = "observation_date") %>%
  inner_join(copper, by = "observation_date") %>%
  inner_join(BAA10Y, by = "observation_date")

# ── Format output ──────────────────────────────────────────────────────────────
merged_data$observation_date = format(merged_data$observation_date, "%Y-%m")

# ── Save processed data ────────────────────────────────────────────────────────
# Two CPI values in 2025 are missing due to the 2025 lapse in appropriations.
write.csv(merged_data, file = "data/processed/cleanedData.csv")

