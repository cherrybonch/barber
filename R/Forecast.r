#decompose time series, predict trend by arima, random by regression + seasonal
#
forecast <- function(inputVector, period){
   numberOfOfficces <- c(3:5,7:35)

  if (!require("forecast")) {
    install.packages("forecast")
    library("forecast")
  }
  else library("forecast")
  if (!require("caret")) {
    install.packages("caret")
    library("caret")
  }
  else library("caret")
  if (!require("e1071")) {
    install.packages("e1071")
    library("e1071")
  }
  else library("e1071")

  if (!require("plotly")) {
     install.packages("plotly")
     library("plotly")
   }

  office <- data.frame()

  if((length(inputVector) > 1)&(length(inputVector) < 31)){
    p <- 1
    for(i in inputVector){
      if(p == 1) office <- data[(match(inputVector, numberOfOfficces)*1204+1):((match(inputVector, numberOfOfficces)+1)*1204), c(1, 5, 9:11, 16, 18:20)]
      else {
        temporaryData <- data[(match(inputVector, numberOfOfficces)*1204+1):((match(inputVector, numberOfOfficces)+1)*1204) ,c(1, 5, 9:11, 16, 18:20)]
        office[, 3:9] <- office[, 3:9] + temporaryData[, 3:9]
      }
      p <- p + 1
      t <- data.frame(t(apply(office[, 7:9], 1, function(x){return(as.numeric(as.character(x))/sum(as.numeric(as.character(x))))})))
      for(i in 1:3) t[is.na(t[, i]), i] <- 0
      colnames(t) <- c("Class1", "Class2", "Class3")
      predictedVector <- office$Trans_sum
      office <- cbind(office, t)
    }
    office <- office[min(which(office$Trans_sum!=0)):length(office$Trans_sum), ]
    if(period > nrow(office)/2) print("Warning: too long prediction period")
    predictedVector <- office$Trans_sum
      } else {
        if(length(inputVector) == 31){
          office <- data[1:1204 ,c(1, 5, 9:11, 16, 18:23)]
          if(period > nrow(office)/2) print("Warning: too long prediction period")
          predictedVector <- office$Trans_sum
        } else {
          if(length(inputVector) == 1){
            office <- data[(match(inputVector, numberOfOfficces)*1204+1):((match(inputVector, numberOfOfficces)+1)*1204) ,c(1, 5, 9:11, 16, 18:20)]
            t <- data.frame(t(apply(office[, 7:9], 1, function(x){return(as.numeric(as.character(x))/sum(as.numeric(as.character(x))))})))
            for(i in 1:3) t[is.na(t[, i]), i] <- 0
            colnames(t) <- c("Class1", "Class2", "Class3")
            office <- cbind(office, t)
            office <- office[min(which(office$Trans_sum!=0)):length(office$Trans_sum), ]
            if(period > nrow(office)/2) print("Warning: too long prediction period")
            predictedVector <- office$Trans_sum
          }
        }
      }

  TimeSeries <- msts(predictedVector, seasonal.periods = seasonalPeriod)
  ValidationSeries <- msts(predictedVector[1:(length(predictedVector) - period)], seasonal.periods = seasonalPeriod)

  #####################################################
  #holt-winters
  modelHW <- HoltWinters(TimeSeries, seasonal = "multiplicative", alpha = 0.1, beta = FALSE, gamma = FALSE)
  predictHW <- predict(modelHW, n.ahead = period)
  validationModelHW <- HoltWinters(TimeSeries[1:(length(TimeSeries) - period)], seasonal = "multiplicative", alpha = 0.1, beta = FALSE, gamma = FALSE)
  validationHW <- predict(validationModelHW, n.ahead = period)
  #exponential smoothing
  predictHolt <- holt(TimeSeries, h = period, alpha = 0.1)
  validationHolt <- holt(TimeSeries[1:(length(TimeSeries) - period)], h = period, alpha = 0.1)
  #arima from package forecast
  trendModel <- auto.arima(TimeSeries)
  predictedARIMA <- forecast.Arima(object = trendModel, h = period)
  trendModelValidation <- auto.arima(TimeSeries[1:(length(TimeSeries) - period)])
  validationARIMA <- forecast.Arima(object = trendModelValidation , h = period)
  #adaptive composition of models
  Validations <- matrix(c(validationHW, validationHolt$mean, validationARIMA$mean), ncol = 3, nrow = period)
  err <- matrix(0, nrow = period, ncol = 3)
  esm <- matrix(0, nrow = period, ncol = 3)
  esmAbs <- matrix(0, nrow = period, ncol = 3)
  weight <- rep(0, 3)
  gamma <- 0.1
  weightedPredict <- rep(0, period)
  #
    for(j in 1:3){
      err[, j] <- TimeSeries[(length(TimeSeries) - period + 1):length(TimeSeries)] - Validations[, j] #error
      for(t in 2:period){
        esmAbs[1, j] <- abs(err[1, j])
        esmAbs[t, j] <- gamma*abs(err[t, j]) + (1 - gamma)*esmAbs[t - 1, j] #exponential smoothing mean of absolute of error
      }
    }
    #choosing weights
    for(j in 1:3) weight[j] <- esmAbs[t, j]^(-1) / sum((esmAbs[t, ])^(-1))

    #making predictions
    for(j in 1:3) weightedPredict <- weight[j]*Predictions[, j] + weightedPredict

  pl <- data.frame(data = c(TimeSeries,weightedPredict))
  plot <- plot_ly(data = pl, x = seq(1:1384), y = data)
  plot
}
