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
      

