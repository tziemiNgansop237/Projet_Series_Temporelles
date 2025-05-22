
###############################################################################
#                         Project: Time Series (2025)                        #
#                         Authors: TZIEMI NGANSOP & TANGOUO KUETE IVANA      #
#                         Date   : 22/05/2025                                #
#                         Purpose: ARIMA Analysis and Modeling               #
###############################################################################




#==============================================================================
#                             Loading packages                                #
#==============================================================================

require(zoo) #convenient and easy-to-use time series format
require(tseries)
#install.packages("seastests")
#install.packages("patchwork")

library(seastests)
library(patchwork)
library(dplyr)
library(readr)
library(ellipse)
require(zoo)
library(ggplot2)
library(gridExtra)
library(fUnitRoots)
library(forecast)
library(tseries)
library(lubridate)
library(urca)

# Loading and cleaning the dataset 
### Setting the working directory
###Replace this directory with your own
chemin <- "C:/Users/damso/Desktop/p/IPI_202501.csv"
base <- read_csv(chemin)

summary(base)
#==============================================================================#

# **************************************************************************#
#                                                                           #
#                                                                           #
#                         PART I: The Data                                  #
#                                                                           #
#                                                                           #
# **************************************************************************#
#==============================================================================#

# ========================================================================= #
# ============================= Question 1 ================================ #
# ========================================================================= #


###______-formatting the variables of interest-______#

## Selecting our variables of interest: time and IT data (CI)
data <- base %>% select(c(`NAF rev. 2`,`CZ`))

## Renaming columns 
colnames(data)<-c("Periode","indice")

## Formatting the time variable: "Year-Month-Day"
data$Periode <- substr(data$Periode, 2, nchar(data$Periode))
data$Periode <-  as.Date(paste0(data$Periode, "01"), format = "%Y%m%d")


#####___Defining training and test datasets__#

##_________Since ARMA is a short-term model  
## we will keep 12 observations for testing 

#     '''Test: Data From 01-01-2024 to 01-01-2025'''
##     '''Training: Data Before 01-01-2024'''


# Training data
train_data <- data %>%
  filter(Periode < as.Date("2024-01-01"))

# Test data
test_data <- data %>%
  filter(Periode >= as.Date("2024-01-01"))


#    ''' Copy of the train_data '''

train_copie <- train_data
#==============================================================================#

# ========================================================================= #
# =============================== Question 2 ============================== #
# ========================================================================= #


#==============================================================================#
### 2-a) Checking the seasonality of the time series

## Since this is a monthly series, we will check for seasonality 
## over 12 months (1 year) using the QS test (Q-statistic test)
ts_data <- ts(data$indice, frequency = 12)  
isSeasonal(ts_data, test = "qs")  # "qs" = test based on intra-period variance

###########################################################
# This test indicates that the series is not seasonal,
# so we will focus on stationarity.
#_____________________________#___________________________#


#==============================================================================#
### 2-a) Data visualization before processing

windows()

ggplot(train_data, aes(x = Periode, y = indice)) +
  geom_line(color = "black") +  # Series color
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "", y = "Index values", title = "") +
  theme(panel.background = element_rect(fill = "lightgrey"))  

#### Checking simple and partial autocorrelations of the raw series

windows()
par(mfrow=c(1,2))
acf(train_data$indice,40,main="");pacf(train_data$indice,40,main="")


# we suspect the presence of a deterministic trend, so we will 
# perform an Augmented Dickey-Fuller stationarity test
#==============================================================================#



##### 2-b) Stationarity test of the series


# In this section, we will write a function `adftest_valid` that tests, 
# for k from 1 to kmax, the stationarity of the AR using `adftest`, but 
# retains only the one for which the residuals are weak white noise.


### a- Function to check residual autocorrelation 
Qtests <- function(var, k, fitdf = 0) {
  pvals <- apply(matrix(1:k), 1, FUN = function(l) {
    pval <- if (l <= fitdf) NA else Box.test(var, lag = l, type = "Ljung-Box", 
    fitdf = fitdf)$p.value  # The Ljung-Box test checks the null hypothesis 
                            #that autocorrelations up to lag are zero
    return(c("lag" = l, "pval" = pval))
  })
  return(t(pvals))
}


### b - Augmented Dickey-Fuller Test
# Function to find the optimal lag for the Dickey-Fuller regression
adfTest_valid <- function(var, kmax, adftype){
  k <- 0
  noautocorr <- 0 # Are the residuals weak white noise?
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    test_adf <- adfTest(var, lags=k, type=adftype)
    pvals <- Qtests(test_adf@test$lm$residuals, 20, fitdf = length(
                   test_adf@test$lm$coefficients))[,2] # Perform portmanteau test
    if (sum(pvals<0.05,na.rm=T)==0) { # if no p_value falls in the rejection zone
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")


    k <- k+1
  }
  return(test_adf)
}
#==============================================================================#

## Applying the function to find the effective lag for the ADF test
best_adf <- adfTest_valid(var=train_data$indice,kmax=24,adftype="ct") 
best_adf

## --------------The number of lag selected is k=2----------------------------#


### 2-c) Eliminating non-stationarity

# We will difference the original series with a lag of 1 and name the transformed
                                                                 #series d_indice

train_data["d_indice"]=c(NA,diff(train_data$indice,1)) # First-order differencing
train_copie["d_indice"]=c(NA,diff(train_data$indice,1)) # First-order differencing



## Plotting simple autocorrelations of both original and differenced series

windows()
par(mfrow=c(1,2))
acf(train_data$indice,lag.max=24,main="",xaxt = "n")
axis(1, at = 0:24, labels = FALSE)
text(x = 0:24, y = par("usr")[3] - 0.05, labels = 0:24, srt = 90, adj = 1, xpd = TRUE)

acf(train_data$d_indice[-1],lag.max=24,main="",xaxt = "n")
axis(1, at = 0:24, labels = FALSE)
text(x = 0:24, y = par("usr")[3] - 0.05, labels = 0:24, srt = 90, adj = 1, xpd = TRUE)


## 2-d) Stationarity test of the differenced series d_indice

best_dadf <- adfTest_valid(var=train_data$d_indice,kmax=24,adftype="ct") # ADF with 3 lags: residuals OK? OK

best_dadf

### We observe that the differenced series is stationary

## Using other tests to confirm previous results
# KPSS test:
kpss.test(train_data$d_indice[-1]) # stationarity
pp.test(train_data$d_indice[-1])   # stationarity

# Stationarity is clearly confirmed


# ========================================================================= #
# =============================== Question 3 ============================== #
# ========================================================================= #

### Graphically represent the chosen series before and after transformation ----

windows()
# Creating both plots
serie_brute <- ggplot(data = train_copie[-1,], aes(x = Periode, y = indice)) +
  geom_line() +
  labs(x = "", y = "Index", title = "")+
  theme(panel.background = element_rect(fill = "lightgrey"))

serie_diff <- ggplot(data = train_copie[-1,], aes(x = Periode, y = d_indice)) +
  geom_line() +
  labs(x = "", y = "Differenced index", title = "")+
  theme(panel.background = element_rect(fill = "lightgrey"))

# Displaying both plots side by side
grid.arrange(serie_brute, serie_diff, nrow = 1)
#==============================================================================#

# **************************************************************************#
#                                                                           #
#                                                                           #
#                         PART II: ARMA Models                              #
#                                                                           #
#                                                                           #
# **************************************************************************#

#==============================================================================#

# ========================================================================= #
# =============================== Question 4 ============================== #
# ========================================================================= #

### Choosing ARMA(p,q) model and checking its validity

# To choose the p and q orders, we will use autocorrelation (ACF) and partial autocorrelation (PACF) functions

windows()
par(mfrow=c(1,2))
acf(train_data$d_indice[-1],20,main="",xaxt = "n")
axis(1, at = 0:20, labels = FALSE)
text(x = 0:20, y = par("usr")[3] - 0.05, 
     labels = 0:20, srt = 90, adj = 1, xpd = TRUE)
pacf(train_data$d_indice[-1],20,main="",xaxt = "n") # up to 20 months lag
axis(1, at = 1:20, labels = FALSE)
text(x = 1:20, y = par("usr")[3] - 0.02, 
     labels = 1:20, srt = 90, adj = 1, xpd = TRUE)

pmax=2 ; qmax=2


## Estimating all possible models

# Initialize the list of valid models: those with white noise residuals and well-specified p and q orders

# We explore all (p, q) pairs between 0 and pmax/qmax, then we test if the ARIMA(p, 0, q) model:
  
# a. Can be estimated without error.

# b. Produces residuals close to white noise.

# c. Has significant coefficients.

# d. Is reasonably parsimonious (pure AR or MA if possible).

modeles_valides <- list()
nb <-0
for (p in 0:pmax) {
  for (q in 0:qmax) {
    modele <- try(arima(train_data$d_indice[-1], order = c(p, 0, q), 
                        include.mean = FALSE)) # Attempt to estimate ARIMA
    if ((sum(is.na(sqrt(diag(modele$var.coef))))==0) & (all(Qtests(
              modele$residuals, 24, fitdf = p+q)[,2] > 0.05, na.rm = TRUE)
               )) { # test for good residuals and well-estimated standard errors
      validite <- abs(modele$coef/sqrt(diag(modele$var.coef)))>1.96 # model validity test (coefficient significance)
      ajout <- FALSE
      if ((p==0 | q==0) & (validite[length(validite)])) {ajout <- TRUE}
      else if ((validite[p]) & (validite[p+q])) {ajout <- TRUE}
      if (ajout){
        nb <- nb + 1
        modeles_valides <- append(modeles_valides, list(list(paste("model_", nb, sep = ""), p, q)))
      }
    }
  }
}


modeles_valides



# The models that passed these tests are:
#--------MA(0,2), ARMA(1,2), ARMA(2,1)-----------#
# We will project these models on the training data
# to compare and then choose the best model according to 
# AIC and BIC criteria.

#==============================================================================#
# p = 0, q = 2
ma2 <- estim <- arima(train_data$d_indice,c(0,0,2))

# p=1,q=2
ar1ma2 <- estim <- arima(train_data$d_indice,c(1,0,2))

# p=2, q=1
ar2ma1 <- estim <- arima(train_data$d_indice,c(2,0,1))
#==============================================================================#

# Choosing the best model among the 3 valid ones based on AIC and BIC criteria

models <- list(ma2 = ma2, ar1ma2 = ar1ma2, ar2ma1 = ar2ma1)
scores <- sapply(models, function(mod) {
  c(AIC = AIC(mod), BIC = BIC(mod), logLik = logLik(mod))})
print(scores)

# Find the columns with the lowest values
min_aic_model <- colnames(scores)[which.min(scores["AIC", ])]
min_bic_model <- colnames(scores)[which.min(scores["BIC", ])]

# Clean output
cat("Best model according to AIC:", min_aic_model, "\n")
cat("Best model according to BIC:", min_bic_model, "\n")



#-----_______Best models:_______-------#  
#          __ AIC: ARMA(1,2)  __
#           _   BIC: MA(2)    _



# Likelihood ratio test: MA(2) VS ARMA(1,2)
# To choose between the two models, we will perform
# a likelihood ratio test with the following hypotheses:
# H0: the MA(2) model is sufficient (the AR(1) term adds nothing).
# H1: the ARMA(1,2) model is better (the AR(1) term improves the likelihood).

#==============================================================================#

# Compute the likelihood ratio test statistic
ll_ratio <- 2 * (scores["logLik", "ar1ma2"] - scores["logLik", "ma2"])
df_diff <- (length(coef(ar1ma2)) - length(coef(ma2))) 

# Compute the test statistic and p-value
p_value <- 1 - pchisq(ll_ratio, df_diff)
cat("P-value : ", p_value, "\n")

#  -- The p-value (0.04) is less than 5%, so the
# difference between the two models is significant.
#  -- The ARMA(1,2) model better explains
# the data than the MA(2) model.
#==============================================================================#


# ____________________________Model: ARMA(1,2)_________________________________#


#==============================================================================
##---------5. Express the ARIMA(p,d,q) model for the selected series-----------#
#==============================================================================


# _____________________Base model: ARIMA(1,1,2)_____________________________#
arima112 <- arima(train_data$indice,c(1,1,2),include.mean=F)


# Residual validity test
all(Qtests(arima112$residuals, 24, fitdf=3)[,2]>0.05,na.rm = TRUE)
# - Result: TRUE i.e. all p-values are greater than 0.05
# - Conclusion: the residuals show no significant autocorrelation
#   at all tested lags, indicating that the model fits the data well.


# Significance of the coefficients
arima112$coef/sqrt(diag(arima112$var.coef))
abs(arima112$coef/sqrt(diag(arima112$var.coef)))>1.96
# - Purpose: Test if the model's coefficients are greater than 1.96
# - Result: TRUE, all coefficients are greater than 1.96
# - Conclusion: The AR1, MA1, and MA2 components provide
# significant information, so the ARIMA(1,1,2) model appears well specified


## Plot of residuals and their autocorrelations

par(mfrow=c(1,2))

plot(arima112$residuals,main="",xlab="",ylab="",col="darkblue")
acf(coredata(arima112$residuals),lag.max = 20,main="",xaxt = "n")
# Add x-axis ticks from 1 to 20 without labels
axis(1, at = 0:20, labels = FALSE)
# Add rotated labels at 90 degrees
# Use text() function to add rotated labels
text(x = 0:20, y = par("usr")[3] - 0.05, labels = 0:20, srt = 90, adj = 1, xpd = TRUE)


##_______________Check invertibility of the ARIMA(1,1,2) model______________________#


#==============================================================================#
#  To do this, we check whether the roots of the ARMA(1,2) polynomial lie outside the unit circle
# - Start by extracting the ARIMA model coefficients
# - Then select MA1 and MA2 coefficients
#  - Then compute the polynomial roots

summary(arima112)
# Extract coefficients
coefs <- arima112$coef
theta <- c(1, coefs["ma1"], coefs["ma2"])

# Compute the polynomial roots
roots <- polyroot(theta)
print(roots)

#  - The roots are c(1.409179+0i, -3.131788-0i)
#  - Their moduli are 1.409179 and 3.131788 respectively
#  - All roots lie outside the unit circle

#___Conclusion: The ARIMA(1,1,2) model is identifiable and correctly estimated____#
#___It provides reliable predictions and is robust to past errors____#


#==============================================================================#

#******************************************************************************#
#*                        Part 3: Forecasting                                 *#
#******************************************************************************#


#==============================================================================#
#6. Write the equation satisfied by the confidence region at level α
#___for the future values (XT+1, XT+2).
#**************************(see Report)*************************************#

#7. Specify the assumptions used to obtain this region.

#**************************(see Report)*************************************#

# Convert the series to data frame for ggplot
data_frame <- data.frame(values = as.numeric(arima112$residuals))


#==============================================================================#
#__8. Graphically represent this region for α = 95%. Comment.

#___---- FORECAST OVER 5 MONTHS ----___#
#__Forecast with confidence intervals
pred <- forecast(arima112, h = 5, level = 95)

#__Display of predicted values and intervals
print(pred)

# ---- EXTRACT TEST DATA ----
# First, check if test_data$d_indice contains at least 5 values
test_vals <- test_data$indice[1:5]

#________---- PLOT THE GRAPH ----_______#


#__Full series to plot
observed_values <- c(train_data$indice, rep(NA, 5))
test_values <- c(rep(NA, length(train_data$indice)), test_vals)
predicted_values <- c(rep(NA, length(train_data$indice)), as.numeric(pred$mean))
IC_up <- c(rep(NA, length(train_data$indice)), pred$upper)
IC_low <- c(rep(NA, length(train_data$indice)), pred$lower)


#__Create a sequence of monthly dates
date_index <- seq(as.Date("1990-01-01"), by = "month", length.out = length(train_data$d_indice) + 5)

#__Build a DataFrame for the plot
df_plot <- data.frame(
  Date = date_index,
  Train = c(train_data$indice, rep(NA, 5)),
  Test = c(rep(NA, length(train_data$indice)), head(test_data$indice, 5)),
  Forecasts = c(rep(NA, length(train_data$indice)), as.numeric(pred$mean)),
  IC_low = c(rep(NA, length(train_data$indice)), as.numeric(pred$lower)),
  IC_up = c(rep(NA, length(train_data$indice)), as.numeric(pred$upper))
)

#==============================================================================#

pred <- forecast(arima112, h = 2, level = 95)

# Assume a correlation rho between the two forecasts
rho <- 0.8  # assumed correlation
sigma1 <- (pred$upper[1] - pred$lower[1]) / (2 * 1.96)
sigma2 <- (pred$upper[2] - pred$lower[2]) / (2 * 1.96)

# Compute covariance from standard deviations and correlation
cov_val <- rho * sigma1 * sigma2

# Covariance matrix with correlation
Sigma <- matrix(c(sigma1^2, cov_val,
                  cov_val, sigma2^2), nrow = 2)

# Forecast center
mu <- c(as.numeric(pred$mean[1]), as.numeric(pred$mean[2]))


ellipse_pts <- ellipse(Sigma, centre = mu, level = 0.95)
ellipse_df <- as.data.frame(ellipse_pts)

#==============================================================================#
#___Time series plot with forecasts and confidence intervals


################################################
#==============================================================================#
#___Confidence interval plot: ellipse
p1 <- ggplot(ellipse_df, aes(x = x, y = y)) +
  geom_path(color = "black") +
  geom_point(aes(x = mu[1], y = mu[2]), color = "black", size = 2) +
  labs(x = "Forecast March 2022", y = "Forecast April 2022",
       title = "Joint confidence ellipse") +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

#==============================================================================#
#___Time series plot with forecasts

p2 <- ggplot(df_plot, aes(x = Date)) +
  geom_line(aes(y = Train), color = "black", size = 0.4) +
  geom_line(aes(y = Test), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(aes(y = Forecasts), color = "red", size = 0.6) +
  geom_ribbon(aes(ymin = IC_low, ymax = IC_up), fill = "gray80", alpha = 0.5) +
  annotate("rect",
           xmin = df_plot$Date[length(train_data$indice) + 1],
           xmax = max(df_plot$Date),
           ymin = -Inf, ymax = Inf,
           fill = "gray90", alpha = 0.4) +
  labs(x = "Date", y = "Manufacturing production index",
       title = "Forecasts of the index and confidence intervals") +
  theme_minimal(base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "none"
  )

# display both graphs in one window
p2 + p1 

#==============================================================================#

#__9. Conditions to improve forecasting

#**************************(see Report)************************************#

#------------------------------------------------------------------------------#
#==================================-END-=======================================#
#------------------------------------------------------------------------------#
