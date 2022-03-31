# whole formula can be generated

    Code
      formula_parser(xc, quo(Xt.pr.13C), quo(Xt.pr.12C))
    Output
      
      Call:
      lm(formula = Xt.pr.13C ~ Xt.pr.12C, data = data)
      
      Coefficients:
      (Intercept)    Xt.pr.12C  
         14.32184      0.01051  
      

---

    Code
      formula_parser(xc, quo(Xt.pr.13C), quo(Xt.pr.12C), type = "Rm")
    Output
      
      Call:
      lm(formula = Xt.pr.13C ~ -1 + Xt.pr.12C, data = data, weights = 1/Xt.pr.12C)
      
      Coefficients:
      Xt.pr.12C  
        0.01098  
      

---

    Code
      formula_parser(xc, quo(Xt.pr.13C), quo(Xt.pr.12C), quo(flag), type = "Rm")
    Output
      
      Call:
      lm(formula = Xt.pr.13C ~ -1 + Xt.pr.12C * flag, data = data, 
          weights = 1/Xt.pr.12C)
      
      Coefficients:
                    Xt.pr.12C            flagconfluent            flagdivergent  
                      0.01054                 13.49926                374.36711  
      Xt.pr.12C:flagdivergent  
                     -0.01053  
      

---

    Code
      formula_parser(xc, quo(Xt.pr.13C), quo(Xt.pr.12C), type = "GLS",
      transformation = "ppt")
    Output
      Generalized least squares fit by REML
        Model: Xt.pr.13C ~ -1 + I(Xt.pr.12C/1000) 
        Data: data 
        Log-restricted-likelihood: -18982.78
      
      Coefficients:
      I(Xt.pr.12C/1000) 
                10.9742 
      
      Variance function:
       Structure: fixed weights
       Formula: ~1/I(Xt.pr.12C/1000) 
      Degrees of freedom: 3900 total; 3899 residual
      Residual standard error: 173.5608 

---

    Code
      formula_parser(tibble::add_column(xc, execution = 1), quo(Xt.pr.13C), quo(
        Xt.pr.12C), type = "LME", nest = quo(file.nm), transformation = "ppt")
    Output
      Linear mixed-effects model fit by REML
        Data: data 
        Log-restricted-likelihood: -18982.78
        Fixed: Xt.pr.13C ~ -1 + I(Xt.pr.12C/1000) 
      I(Xt.pr.12C/1000) 
                10.9742 
      
      Random effects:
       Formula: ~-1 + I(Xt.pr.12C/1000) | execution
              I(Xt.pr.12C/1000)
      StdDev:         0.2411059
      
       Formula: ~-1 + I(Xt.pr.12C/1000) | file.nm %in% execution
              I(Xt.pr.12C/1000) Residual
      StdDev:         0.2411059 173.5608
      
      Variance function:
       Structure: fixed weights
       Formula: ~1/I(Xt.pr.12C/1000) 
      Number of Observations: 3900
      Number of Groups: 
                   execution file.nm %in% execution 
                           1                      1 

