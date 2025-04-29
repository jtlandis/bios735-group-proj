##-------------------------------------------##
# Source ####
library(forecast)
library(ggplot2)
library(tscount)
##-------------------------------------------##

##-------------------------------------------##
# Arma model ####
model_fit_TS = function(dat, item, brand,
                        p, q){
  
  fit = Arima(dat[dat$item == item&dat$brand == brand,]$QTY, order = c(p,0,q),
              xreg = dat[dat$item==item&dat$brand==brand,]$PROMO)
  
  return(fit)
  
}

predict_TS = function(fit, ahead, xreg){
  
  fc = forecast(fit, h = ahead, xreg = xreg)
  
  return(fc)
}
##-------------------------------------------##


##-------------------------------------------##
# Data load ####
train = readRDS("data_set_tidy_train.rds")
test = readRDS("data_set_tidy_test.rds")
##-------------------------------------------##


##-------------------------------------------##
# Fit procedure ####
fit = model_fit_TS(train, 1, "B1", 3, 2)
train.sub = subset(train, item==1 & brand=="B1")
fit.val = fit$fitted
fit.val[fit.val<0] = 0
fit.val = round(fit.val)
# in-sample performance
ts.plot(
  train.sub$QTY,
  fit.val,
  col = c("black","red"),
  lty = c(1,2),
  lwd = 2,
  xlab   = "Time",
  ylab   = "Quantity",
  main   = sprintf("ARMA(%d,%d) In‐Sample Fit", 3, 2),
  cex.lab = 1.5,
  cex.axis = 1.5
)
legend(
  "topleft",
  legend = c("Actual","Fitted"),
  col = c("black","red"),
  lty = c(1,2),
  lwd = 2,
  bty    = "n"
)
##-------------------------------------------##


##-------------------------------------------##
# Prediction procedure ####

# 1-step prediction
test.sub = subset(test, item==1 & brand=="B1")
test.ts = ts(test.sub$QTY)
train.ts = ts(train.sub$QTY)
H = length(test.ts)
pred.one = numeric(H)      # to store each 1-step forecast
fith = fit       # local copy

for(i in seq_len(H)) {
  # a) forecast exactly 1 step ahead
  fc = forecast(fith, h = 1, xreg = test.sub$PROMO[i])
  pred.one[i] = as.numeric(fc$mean)
  
  # b) append the *true* test value, refit via model recycling
  y.hist = ts(c(as.numeric(train.ts), as.numeric(test.ts[1:i])),
                   start = start(train.ts),
                   frequency = frequency(train.ts))
  fith = Arima(y.hist, model = fith, xreg = c(train.sub$PROMO,test.sub$PROMO[1:i]))
}

pred.one[pred.one<0] = 0
pred.one = round(pred.one)
# out-sample performance
ts.plot(
  test.sub$QTY,
  ts(pred.one),
  col = c("black","red"),
  lty = c(1,2),
  lwd = 2,
  xlab   = "Time",
  ylab   = "Quantity",
  main   = sprintf("ARMA(%d,%d) Out‐Sample One-step Prediction", 3, 2),
  cex.lab = 1.5,
  cex.axis = 1.5
)
legend(
  "topleft",
  legend = c("Actual","Predicted"),
  col = c("black","red"),
  lty = c(1,2),
  lwd = 2,
  bty    = "n"
)


# apply R package
pred = predict_TS(fit, ahead = nrow(test.sub), xreg = test.sub$PROMO)
pred = pred$mean
pred[pred<0] = 0
pred = round(pred)
# out-sample performance
ts.plot(
  test.sub$QTY,
  ts(pred),
  col = c("black","red"),
  lty = c(1,2),
  lwd = 2,
  xlab   = "Time",
  ylab   = "Quantity",
  main   = sprintf("ARMA(%d,%d) Out‐Sample Prediction", 3, 2),
  cex.lab = 1.5,
  cex.axis = 1.5
)
legend(
  "topleft",
  legend = c("Actual","Predicted"),
  col = c("black","red"),
  lty = c(1,2),
  lwd = 2,
  bty    = "n"
)
##-------------------------------------------##
