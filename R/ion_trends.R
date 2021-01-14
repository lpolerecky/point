#' Predicting trends in ionization efficiency
#'
#' Fluctuations in electronics, development of the sputter pit geometry and the
#' analysed substrate can all cause trends and fluctuations in the secondary ion
#' current. This function attempts to accommodate the global trend in the
#' ionization trend by application of a GAM model. The nested variant can then
#' be used to gauge whether a set of analyses are likely to originate from a
#' homogenoues substrate at the level of individual analysis.
#'
#' @param .df A tibble containing raw ion count data.
#' @param ... Variables for grouping.
#' @param .Xt A variable constituting the ion count rate.
#' @param .N A variable constituting the ion counts.
#' @param .species A variable constituting the species analysed.
#' @param .t A variable constituting the time increments.
#' @param .nest A logical indicating whether a nested mixed GAM model is
#' applied.
#' @param .group The grouping variable for the mixed GAM model.
#' @param .plot Logical indicating whether to plot ion trends
#' @param .method Method for calculating de-trended single ion counts
#' @param .hide A logical indicating whether only processed data should be
#' returned. If \code{TRUE} The model parameters are contained as an attribute
#' named \code{"modeldata"}.
#' @export
#' @examples
#' # raw ion counts on carbonate
#' tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb_pr <- cor_IC(tb_rw)
#'
#' # remove zero count analysis
#' tb_0 <- zeroCt(tb_pr, "12C", "40Ca 16O", file.nm, .warn = FALSE)
#'
#' # predict ionization trends
#' predict_ionize(tb_0, file.nm)
#'
#' # predict ionization trends mixed type
#' tb_bel <- filter(tb_0, str_detect(sample.nm, "Belemnite"))
#'
#' predict_ionize(tb_bel, sample.nm, file.nm, .nest = TRUE, .group = file.nm)
predict_ionize <- function(.df, ..., .Xt = Xt.pr, .N = N.pr,
                           .species = species.nm, .t = t.nm,
                           .nest = FALSE, .group = NULL,
                           .plot = TRUE, .method = "median", .hide = TRUE
                           ){

  species <- enquo(.species)
  gr_by <- enquos(..., .species, .named = TRUE)
  group <- enquo(.group)
  Xt <- enquo(.Xt)
  N <- enquo(.N)
  t <- enquo(.t)

  if (.nest) nst_by <- gr_by[!sapply(gr_by, as_name) %in% as_name(group)]

  # new quosures
  if (.nest) {
    Xt.mdl <- quo_updt(Xt, post = "grp", update_post = TRUE) # group-wise
    M_Xt.l0 <- quo_updt(Xt, pre = "M", post = "grp", update_post = TRUE)
    } else {
      Xt.mdl <- quo_updt(Xt, post = "nlt", update_post = TRUE) # per analysis
      M_Xt.l0 <- quo_updt(Xt, pre = "M", post = "nlt", update_post = TRUE)
  }

  Xt.l0 <- quo_updt(Xt, post = "l0", update_post = TRUE)
  N.l0 <-  quo_updt(N, post = "l0", update_post = TRUE)

  mth_switch <- function(mth, arg) {
    switch(
      mth,
      median = call2("median", arg),
      mean  = call2("mean", arg),
      stop("unknown method", call. = FALSE)
      ) %>%
      eval()
    }

  # metadata
  if (ncol(select(.df, ends_with(".mt"))) == 0) df <- unfold(.df)


# Groupwise
  if (.nest) {
    df <- nest(df, data = -c(!!! nst_by)) %>%
      mutate(gam_out = purrr::map(data, ~gam_fun(.x, Xt, t, group, Xt.mdl))) %>%
      unnest(cols = c(data, gam_out)) %>%
      group_by(!!! nst_by) %>%
      mutate(
        # groupwise mean
        !! M_Xt.l0 := mth_switch(.method,!! Xt.mdl),
        # mean plus random intercept correction
        !! Xt.l0 := !! M_Xt.l0 + (!! Xt - !! Xt.mdl) + ran_in.ml,
        !! N.l0 := !! Xt.l0 * (min(!! t) - .data$tc.mt)
        )
# Single analysis
    } else {
      df <- nest(df, data = -c(!!! gr_by)) %>%
        mutate(
          !! Xt.mdl := purrr::map(data, ~as.vector(gam_form(.x, Xt, t)))
          ) %>%
        unnest(cols = c(data, !! Xt.mdl)) %>%
        group_by(!!! gr_by) %>%
        mutate(
          !! M_Xt.l0 := mth_switch(.method,!! Xt.mdl),
          # mean correction
          !! Xt.l0 := !! M_Xt.l0 + (!! Xt - !! Xt.mdl),
          !! N.l0:= !! Xt.l0 * (min(!! t) - .data$tc.mt)
        )
      }

  df <- ungroup(df)

  if (.plot) {
    if(.nest) {
      facets_gr <- quos(!!! nst_by)
    } else {
      facets_gr <- quos(!!! gr_by)
    }

    # De-trended timeseries data
    df_pl <- mutate(
      df,
      hat_Y = !! Xt.mdl,
      hat_min = !! Xt.mdl,
      hat_max = !! Xt.mdl,
      flag = "good"
      )

    plot_args <- list2(
      df = df_pl,
      stat = NULL,
      y = Xt,
      x = t,
      ion1 = NA,
      ion2 = NA,
      plot_title = "timeseries",
      plot_type = "static",
      args = NULL,
      iso = FALSE,
      isoscale = NULL,
      !!! facets_gr
    )

    # control with environment
    data_env <- env(data = df)
    eval(
      expr(
        gg_default(!!! plot_args) +
          geom_hline(aes(yintercept = !! M_Xt.l0), color = "red", size = 1.1)
        ),
         data_env
      ) %>%
      print()


    # hide model data
    if (.hide) df <- fold(df, c(".mt",".rw", ".ml"))
    return(df)
    }
  # hide model data
  if (.hide) df <- fold(df, c(".mt",".rw", ".ml"))
  return(df)
}


