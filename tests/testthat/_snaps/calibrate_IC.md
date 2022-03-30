# isotope values can be converted

    Code
      calib_R(0.0111, reference = 0.011237, type = "composition", input = "R",
        output = "delta")
    Output
      [1] -12.19187

---

    Code
      calib_R(0.0111, reference = "VPDB", isotope = "13C", type = "composition",
        input = "R", output = "delta")
    Output
      [1] -12.19187

---

    Code
      calib_R(0.0111, reference = 0.011237, type = "composition", input = "R",
        output = "F")
    Output
      [1] 0.01097814

---

    Code
      calib_R(-25, reference = "VPDB", isotope = "13C", type = "enrichment", input = "delta",
        output = "alpha", y = -105)
    Output
      [1] 0.9179487

