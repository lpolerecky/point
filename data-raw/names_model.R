## code to prepare `DATASET` dataset goes here

names_model <- tibble::tibble(
  name = c(
    "M_R_*",
    "f",
    "f_cl",
    "f_cu",
    "F_vl",
    "p_F",
    "M_R_*",
    "GM_R",
    "dNorm",
    "RS_R_inter",
    "RS_R_se_inter",
    "dAIC_inter",
    "p_inter",
    "RS_R_intra",
    "RS_R_se_intra",
    "p_intra"
    ),
  latex = c(
    "$\\bar{R}$",
    "$\\text{Cohen's} f$",
    "$\\[$",
    "$\\]$",
    "$F_{q,n–k–1}$",
    "$p_{R'}$",
    "$\\bar{\\sum{R}}$",
    "$\\bar{\\prod{R}}$",
    "$\\delta_{\\prod{R} - \\sum{R}}$ (\u2030)",
    "$\\hat{\\epsilon_{\\sum{R}}$ (\u2030)",
    "$s_{\\hat{\\epsilon_{\\sum{R}}}$ (\u2030)",
    "$\\Delta AIC$",
    "$p_{\\sum{R}}$",
    "$\\hat{\\epsilon_{\\sum{R'}}$ (\u2030)",
    "$s_{\\hat{\\epsilon_{\\sum{R}}}$ (\u2030)",
    "$p_{\\sum{R'}}$"
    ),
  type = c(
    rep("Ordinary Least Squares", 6),
    rep("Restricted Maximum Likelihood optimization", 10)
    ),
  label = c(
    "mean R of analysis",
    "cohen's f",
    "cohen's f lower bound CI",
    "cohen's f upper bound CI",
    "F test statistic of joint model hypothesis test of single analysis intra-variation",
    "p value of joint model hypothesis test of single analysis intra-variation",
    "mean R of group of analyses",
    "geometric mean R of group of analysis",
    "relative difference in per mille between geometric and arithmetic mean",
    "predicted relative standard deviation of a set of analyses (inter-variation ~ external reproducibility)",
    "standard error of the predicted relative standard deviation of a set of analyses (inter-variation ~ external reproducibility)",
    "difference in AIC between full and reduced model in likelihood ratio test for inter-variation",
    "p value of likelihood ratio test for inter-variation",
    "predicted relative standard deviation of a set of augmentation cycles (longitudinal test for intra-variation)",
    "standard error of the predicted relative standard deviation of a set of augmentation cycles (longitudinal test for intra-variation)",
    "p value of likelihood ratio test for intra-variation"
    )
)

usethis::use_data(names_model, overwrite = TRUE)
