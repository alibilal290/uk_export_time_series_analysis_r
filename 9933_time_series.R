# Load Libraries and data set ################################################

library(magrittr)
library(ggplot2)

exports_ds <- read.csv("series-190324.csv", header = FALSE)
exports_ds <- as.numeric(c(exports_ds[-(1:138), -1])) # No row deleted

# EDA and Pre-processing ######################################################

exports_ds <- ts(exports_ds, start = 1998, frequency = 12)
exports_ds
print(summary(exports_ds))
str(exports_ds)
anyNA(exports_ds, class)

################################# Visualise Data set ###

# Convert the time series to a data frame
uk_exports_df <- data.frame(
  Date = seq(as.Date("1998-01-01"), by = "month", length = length(exports_ds)),
  Exports = as.numeric(exports_ds)
)

# Extracting Year and Month from the Date
uk_exports_df$Year <- format(uk_exports_df$Date, "%Y")
uk_exports_df$Month <- format(uk_exports_df$Date, "%m")

# Check the first few rows to ensure correctness
head(uk_exports_df)
mean_exports <- mean(uk_exports_df$Exports)

# Line plot
ggplot(data = uk_exports_df, aes(x = Date, y = Exports)) +
  geom_line(color = "darkblue") +
  geom_hline(yintercept = mean_exports, linetype = "dashed", color = "red") + # Add mean line
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "MBS 25 UK Export Data", x = "Year", y = "Exports (£M)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Box plot
ggplot(data = uk_exports_df, aes(x = Year, y = Exports)) +
  geom_boxplot() +
  labs(title = "MBS 25 UK Export Data", x = "Year", y = "Exports (£M)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Log Transformation ##########################################################

exports_ds_log <- log(exports_ds)

# Adding a log-transformed column
uk_exports_df$Log_Exports <- log(uk_exports_df$Exports)

ggplot(data = uk_exports_df, aes(x = Date, y = Log_Exports)) +
  geom_line(color = "darkblue") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Log Transformed UK Exports Data", x = "Year", y = "Log Values") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Split Train and test Data set ###############################################

months_for_test <- (12*5) + 1 # 1 month of 2024
training_data_log <- exports_ds_log[1:(length(exports_ds_log) - months_for_test)]
test_data_log <- exports_ds_log[(length(exports_ds_log) - months_for_test + 1):length(exports_ds_log)]
length(training_data_log)
length(test_data_log)

# De-trend Data ###############################################################

TIME <- seq(from = 1998, by = 1/12, length = length(exports_ds_log))
TIME.TRAIN <- seq(from = 1998, by = 1/12, length = length(training_data_log))
TIME.TEST <- TIME[(length(TIME) - months_for_test + 1):length(TIME)]

length(TIME)
length(TIME.TEST)
length(TIME.TRAIN)

trend1 <- lm(training_data_log ~., data.frame(TIME = poly(TIME.TRAIN, degree = 1, raw = TRUE)))
trend2 <- lm(training_data_log ~., data.frame(TIME = poly(TIME.TRAIN, degree = 2, raw = TRUE)))
trend3 <- lm(training_data_log ~., data.frame(TIME = poly(TIME.TRAIN, degree = 3, raw = TRUE)))

print(summary(trend1))
print(summary(trend2))
print(summary(trend3))

print(AIC(trend1))
print(AIC(trend2))
print(AIC(trend3))
# Both 2 and 3 trends have same AIC Values, Why? [1] -185.3638
plot(TIME.TRAIN, training_data_log, type = "l", main = "Trends of UK Exports - 1998 to 2023",
     xlab = "Year", ylab = "Logged Values")
lines(TIME.TRAIN, fitted(trend1), col = 'red', lwd = 2)
lines(TIME.TRAIN, fitted(trend2), col = 'green', lwd = 2)
lines(TIME.TRAIN, fitted(trend3), col = 'blue', lwd = 2)
#legend()


# Remove trend from data
data_notrend <- training_data_log - fitted(trend1)
plot(TIME.TRAIN, data_notrend, type = "l", main = "De-trended UK Exports", xlab = "Year", ylab = "Logged Values")


# Fit Simple seasonality model ################################################
# Create SIN and COS matrices
ncol <- 6
SIN <- COS <- matrix(nrow = length(TIME.TRAIN), ncol = ncol)

for (i in 1:ncol){
  SIN[,i] <- sin(2*pi*i*TIME.TRAIN)
  COS[,i] <- cos(2*pi*i*TIME.TRAIN)
}

seasonality <- lm(data_notrend ~. -1, data.frame(SIN = SIN[, 1:6], COS = COS[, 1:6]))

plot(TIME.TRAIN, training_data_log, type = "l", main = "Seasonality",
     xlab = "Year", ylab = "Logged Values")
lines(TIME.TRAIN, fitted(trend1) + fitted(seasonality), col = "red")
# lines (mean_value)

AIC(seasonality) # AIC Value = [1] -194.9349

# Fit Harmonic Seasonality Model and Compare AIC Values #######################


# Residuals ###################################################################
residuals <- as.numeric(data_notrend) - fitted(seasonality)

residual_plot <- data.frame(TIME.TRAIN, residuals)

ggplot(residual_plot, aes(x = TIME.TRAIN, y = residuals)) +
  geom_line() +
  geom_point(shape = 20) +
  labs(title = "Residuals of Seasonal Model",
       x = "Year",
       y = "Residuals") +
  theme_minimal()

acf(residuals, main = "ACF of Residuals")
pacf(residuals, main = "PACF of Residuals")

# Since the residuals are NOT stationary, Differencing applied within ARIMA model with degree =1

# Comprehensive ADF Test WO library ###########################################

# AIC Calculation for ARIMA Models ############################################
n <- length(residuals)
norder <- 6
p <- q<- c(1:norder) -1  # This defines the orders of the AR and MA components
aic <- matrix(0, norder, norder, dimnames = list(paste0("AR", p), paste0("MA", q)))

# Compute AIC for different ARIMA models
for (i in 1:norder){
  for (j in 1:norder){
    modij <- arima(residuals, c(p[i], 0, q[j]), method = "ML", optim.control = list(maxit = 1000))
    aic[i,j] <-modij$aic - 2 *(p[i] + q[j] +1) + 2 *(p[i] + q[j] +1) * n / (n - p[i] - q[j] -2)
  }
}

# Plot AIC values 
aic_vector <- as.vector(aic)
plot(aic_vector, ylab = "AIC Values", type = "l", main = "AIC Values")

# Print the AIC matrix and the best ARIMA model order
print(aic)
min_aic_value <- min(aic)
best_arima_order <- which(aic == min_aic_value, arr.ind = TRUE)
print(paste("Best ARIMA(p,d,q) Model: ARIMA(", rownames(aic)[best_arima_order[1]], 
            ",0,", colnames(aic)[best_arima_order[2]], ")", sep = ""))

#[1] "Best ARIMA(p,d,q) Model: ARIMA(AR3,0,MA4)"
# [1] "Best ARIMA(p,d,q) Model: ARIMA(AR1,0,MA1)" Latest

# Fitting ARIMA Model #########################################################

# Extract the ARIMA order with the lowest AIC value
aic_index <- which(aic == min(aic), arr.ind = TRUE)
porder <- aic_index[1,1] -1
qorder <- aic_index[1,2] -1

residuals.model <- arima(residuals, order = c(porder, 0, qorder), method = "ML",
                         optim.control = list(maxit = 10000))

summary(residuals.model)

# Producing diagnostic plots ##################################################
par(mfrow=c(2,2)) # Set up the graphics layout to display 4 plots at once
plot(residuals.model$residuals, main="Residuals of the ARIMA Model", ylab="Residuals")
acf(residuals.model$residuals, main="ACF of Residuals")
pacf(residuals.model$residuals, main="PACF of Residuals")
tsdiag(residuals.model) # Provides additional diagnostic plots
par(mfrow=c(1,1))

# Shapiro-Wilk Normality Test on ARIMA Model Residuals
shapiro.test(residuals.model$residuals)

# Out-of-Sample Testing #######################################################

n_test <- length(test_data_log)  # Define the length of your test data

# Forecasting using the ARIMA model
forecasts <- predict(residuals.model, n.ahead = n_test)

# Extract forecasted values
trend1 <- lm(training_data_log ~ TIME.TRAIN)
trend_test <- predict(trend1, newdata = data.frame(TIME.TRAIN = TIME.TEST))
seasonality_model <- lm(data_notrend ~ . -1, data.frame(SIN = SIN[, 1:ncol], COS = COS[, 1:ncol]))

SIN.TEST <- matrix(nrow = length(TIME.TEST), ncol = ncol)
COS.TEST <- matrix(nrow = length(TIME.TEST), ncol = ncol)
for (i in 1:ncol) {
  SIN.TEST[,i] <- sin(2*pi*i*TIME.TEST)
  COS.TEST[,i] <- cos(2*pi*i*TIME.TEST)
}

# Now predict the seasonality for the test period
seasonality_test <- predict(seasonality_model, newdata = data.frame(SIN = SIN.TEST, COS = COS.TEST))

# Forecasting residuals for the test period
forecasts <- predict(residuals.model, n.ahead = n_test)
forecasted_residuals <- forecasts$pred

# Assuming the trend and seasonality components are calculated as described above
forecasted_values <- exp(forecasted_residuals + trend_test + seasonality_test)

# Calculate RMSE manually
actual_test_data <- exp(test_data_log)  # Reverse the log transformation of the test data
rmse <- sqrt(mean((forecasted_values - actual_test_data)^2))

# Display the RMSE
print(paste("RMSE for out-of-sample forecasts:", rmse))

#Visual comparison using ggplot2
test_period <- length(training_data_log) + (1:length(test_data_log))

# Data frame for plotting
plot_data <- data.frame(
  Time = test_period,
  Actual = actual_test_data,
  Forecast = forecasted_values
)

ggplot(plot_data, aes(x = Time)) +
  geom_line(aes(y = Actual, colour = "Actual")) +
  geom_line(aes(y = Forecast, colour = "Forecast")) +
  labs(title = "Out-of-Sample Forecast vs Actual",
       y = "Values", colour = "Legend") +
  theme_minimal()
