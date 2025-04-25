# par_model: Intercept only

    Code
      print(res)
    Output
      $beta
      named numeric(0)
      
      $gamma
      (Intercept) 
                3 
      

# par_model: lags only

    Code
      print(res)
    Output
      $beta
           lag1      lag2      lag3 
      0.3302994 0.2160618 0.2842066 
      
      $gamma
      named numeric(0)
      

# par_model: lags and covar

    Code
      print(res)
    Output
      $beta
               lag1          lag2          lag3 
       0.3393938358 -0.0004108964  0.3352749123 
      
      $gamma
      (Intercept)       PROMO 
       0.80049809  0.08179525 
      

# par_model: estimate covariates by ave items

    Code
      print(res)
    Output
      $beta
              lag1         lag2 
       0.574515755 -0.004101278 
      
      $gamma
      (Intercept)       PROMO     brandB2     brandB3     brandB4 
       0.18229385  0.14951757  0.01065082 -0.01781969  0.20425868 
      

