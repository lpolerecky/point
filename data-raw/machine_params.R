## Names for the machine parameter settings

cam_params <- tibble(
  cam = c(
    "Raw data file",
    "Block number",
    "Meas. per block",
    "Rejection(sigma)",
    "Slit Preset",
    "Lens Preset",
    "Pre Sputtering Time (s)",
    "Raster (um)",
    "Scanning On",
    "Blanking",
    "Raster",
    "Tuning Mode",
    "Regulation Mode",
    "Check every (fr.)",
    "Secondary Ion Beam Centering",
    "during Acq",
    "Width Horizontal(V)",
     "Vertical(V)",
    "E0S Centering",
    "Width(V)",
    "E0P Offset (V)",
    "Primary Current before acq",
    "after acq",
    "FC Background before acq",
    "FC Background after acq"
     ),
# New names for internal processing
  point = c(
    "file_raw.mt",
    "bl_num.mt",
    "meas_bl.mt",
    "rejection.mt",
    "slit.mt",
    "lens.mt",
    "presput.mt",
    "rast_com.mt",
    "scan.mt",
    "blank_rast.mt",
    "raster.mt",
    "tune.mt",
    "reg_mode.mt",
    "chk_frm.mt",
    "sec_ion_cent.mt",
    "acq.mt",
    "width_hor.mt",
    "width_ver.mt",
    "E0S_cent.mt",
    "width_V.mt",
    "E0P_off.mt",
    "prim_cur_start.mt",
    "prim_cur_after.mt",
    "FC_start.mt",
    "FC_after.mt"
     ),
  use = c("ID", rep("meta", 24))
  )


usethis::use_data(cam_params, overwrite = TRUE, internal = TRUE)
