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
      0.4171512 0.2074937 0.1876541 
      
      $gamma
      named numeric(0)
      
      $iter
      [1] 12
      
      $objective
      [1] -243874.2
      

# par_model: lags and covar

    Code
      print(res)
    Output
      $beta
           lag1      lag2      lag3 
      0.3952092 0.1516124 0.1225360 
      
      $gamma
      (Intercept)       PROMO 
         1.021168    1.381041 
      
      $iter
      [1] 20
      
      $objective
      [1] -221708.4
      

# par_model: estimate covariates by ave items

    Code
      print(res)
    Output
      $beta
            lag1       lag2 
      0.33338226 0.07670208 
      
      $gamma
      (Intercept)       PROMO     brandB2     brandB3     brandB4 
        0.8916504   1.8616106  -1.3919850  -0.7326483   0.4056676 
      
      $iter
      [1] 23
      
      $objective
      [1] -44374.99
      

