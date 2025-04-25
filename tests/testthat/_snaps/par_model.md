# par_model: Intercept only

    Code
      print(res)
    Output
      $beta
      named numeric(0)
      
      $gamma
      (Intercept) 
         1.945908 
      
      $iter
      [1] 8
      
      $objective
      [1] -590271.8
      

# par_model: lags only

    Code
      print(res)
    Output
      $beta
           lag1      lag2      lag3 
      0.4171832 0.2074944 0.1876104 
      
      $gamma
      named numeric(0)
      
      $iter
      [1] 13
      
      $objective
      [1] -243874.2
      

# par_model: lags and covar

    Code
      print(res)
    Output
      $beta
           lag1      lag2      lag3 
      0.3952482 0.1516265 0.1224757 
      
      $gamma
      (Intercept)       PROMO 
         1.021149    1.381067 
      
      $iter
      [1] 19
      
      $objective
      [1] -221708.4
      

# par_model: estimate covariates by ave items

    Code
      print(res)
    Output
      $beta
            lag1       lag2 
      0.33331122 0.07663363 
      
      $gamma
      (Intercept)       PROMO     brandB2     brandB3     brandB4 
        0.8911661   1.8618050  -1.3918708  -0.7322984   0.4059405 
      
      $iter
      [1] 22
      
      $objective
      [1] -44374.99
      

