## code to prepare `DATASET` dataset goes here

names_model <- tibble::tibble(
  name = c(
    "M",
    "F",
    "p",
    "hat_M",
    # "GM_R",
    # "dNorm",
    "hat_RS",
    # "RS_R_se_inter",
    "dAIC",
    "p" #,
    # "RS_R_intra",
    # "RS_R_se_intra",
    # "p_intra"
    ),
  derived = c(rep("R", 3), rep("M_R", 4)),
  # latex = c(
  #   "$\\bar{R}$",
  #   "$F_{q,n–p–1}$",
  #   "$p_{R'}$",
  #   "$\\hat{\\bar{R}}$",
  #   "$\\hat{\\prod{R}}$",
  #   "$\\delta_{\\hat{\\prod{R}} - \\hat{\\bar{R}}} (\\text{\\textperthousand})$",
  #   "$\\epsilon_{\\hat{R}} (\\text{\\textperthousand})$",
  #   "$s_{\\epsilon_{\\bar{R}}} (\\text{\\textperthousand})$",
  #   "$\\Delta AIC$",
  #   "$p_{\\hat{R}}$",
  #   "$\\epsilon_{\\hat{R'}} (\\text{\\textperthousand})$",
  #   "$s_{\\epsilon_{\\hat{R'}}} (\\text{\\textperthousand})$",
  #   "$p_{\\hat{R'}}$"
  #   ),
  type = c(
    rep("Ratio method", 3),
    rep("Restricted Maximum Likelihood optimization", 4)
    ),
  label = c(
    "mean R of analysis",
    "F test statistic of joint model hypothesis test of single analysis intra-variation",
    "p value of joint model hypothesis test of single analysis intra-variation",
    "mean R of group of analyses",
    # "geometric mean R of group of analysis",
    # "relative difference in per mille between geometric and arithmetic mean",
    "predicted relative standard deviation of a set of analyses (inter-variation ~ external reproducibility)",
    # "standard error of the predicted relative standard deviation of a set of analyses (inter-variation ~ external reproducibility)",
    "difference in AIC between full and reduced model in likelihood ratio test for inter-variation",
    "p value of likelihood ratio test for inter-variation"#,
    # "predicted relative standard deviation of a set of augmentation cycles (longitudinal test for intra-variation)",
    # "standard error of the predicted relative standard deviation of a set of augmentation cycles (longitudinal test for intra-variation)",
    # "p value of likelihood ratio test for intra-variation"
    )
)

usethis::use_data(names_model, overwrite = TRUE)
