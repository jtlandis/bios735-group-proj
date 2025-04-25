# par_model: Intercept only

    Code
      print(res)
    Output
      $beta
      named numeric(0)
      
      $gamma
      (Intercept) 
                3 
      
      $iter
      [1] 2
      
      $objective
      [1] -590277.5
      

# par_model: lags only

    Code
      print(res)
    Output
      $beta
           lag1      lag2      lag3 
      0.3582281 0.1617883 0.1623495 
      
      $gamma
      named numeric(0)
      
      $iter
      [1] 4
      
      $objective
      [1] -246967.5
      

# par_model: lags and covar

    Code
      print(res)
    Output
      $beta
             lag1        lag2        lag3 
       0.38472482 -0.04287368  0.26482992 
      
      $gamma
      (Intercept)       PROMO 
        1.0959704   0.5965945 
      
      $iter
      [1] 6
      
      $objective
      [1] -232420.9
      

# par_model: estimate covariates by ave items

    Code
      print(res)
    Output
      $beta
            lag1       lag2 
      0.67023309 0.08149927 
      
      $gamma
      (Intercept)       PROMO     brandB2     brandB3     brandB4 
       1.84511733  1.11705547 -0.35467703 -0.07485625  0.46713457 
      
      $iter
      [1] 8
      
      $objective
      [1] -49846.45
      

