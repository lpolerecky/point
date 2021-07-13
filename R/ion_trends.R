#' Predicting trends in ionization efficiency
#'
#' Fluctuations in electronics, development of the sputter pit geometry and the
#' analysed substrate can all cause trends and fluctuations in the secondary ion
#' current. This function attempts to accommodate the global trend in the
#' ionization trend by application of a GAM model. The nested variant can then
#' be used to gauge whether a set of analyses are likely to originate from a
#' homogeneous substrate at the level of individual analysis.
#'
#' @param .IC A tibble containing raw ion count data.
#' @param ... Variables for grouping.
#' @param .X A variable constituting the ion count rate.
#' @param .N A variable constituting the ion counts.
#' @param .species A variable constituting the species analysed.
#' @param .t A variable constituting the time increments.
#' @param .bl_t A variable constituting the blanking time.
#' @param .nest A variable identifying a groups of analyses which indicates
#' whether a nested mixed GAM model is applied.
#' @param .plot Logical indicating whether to plot ion trends
#' @param .method Method for calculating de-trended single ion counts
#' @param .hide A logical indicating whether only processed data should be
#' returned. If \code{TRUE} The model parameters are contained as an attribute
#' named \code{"modeldata"}.
#' @export
#' @examples
#'
#' # remove zero count analysis
#' tb_0 <- zeroCt(real_IC, "12C", "40Ca 16O", sample.nm, file.nm, .warn = FALSE)
#'
#' # predict ionization trends
#' predict_ionize(tb_0, file.nm)
predict_ionize <- function(.IC, ..., .nest = NULL, .X = NULL, .N = NULL,
                           .species = NULL, .t = NULL, .bl_t = NULL,
                           .plot = TRUE, .method = "median", .hide = TRUE){

  # Metadata
  if (ncol(select(.IC, ends_with(".mt"))) == 0) .IC <- unfold(.IC)

  # Quoting the call (user-supplied expressions)
  # Additional arguments
  args <- inject_args(
    .IC,
    enquos(.X = .X, .N = .N, .species = .species, .t = .t, .bl_t = .bl_t),
    type = c("processed", "group", "meta")
    )

  # Grouping and nesting
  gr_by <- enquos(...) %>% append(args[[".species"]])
  nest <- enquo(.nest)

  # New quosures
  if (is_symbol(get_expr(nest))) {
    # group-wise
    nest_gr <- gr_by[!sapply(gr_by, as_name) %in% as_name(nest)]
    X.mdl <-quo_updt(args[[".X"]], post = "grp", update_post =  TRUE)
    M_X.l0 <- quo_updt(args[[".X"]], pre = "M", post = "grp", TRUE)
  } else {
    # per analysis
    X.mdl <- quo_updt(args[[".X"]], post = "nlt", update_post = TRUE)
    M_X.l0 <- quo_updt(args[[".X"]], pre = "M", post = "nlt", TRUE)
  }

  X.l0 <- quo_updt(args[[".X"]], post = "l0", update_post = TRUE)
  N.l0 <- quo_updt(args[[".N"]], post = "l0", update_post = TRUE)

  # Group-wise
  if (is_symbol(get_expr(nest))) {
    IC <- tidyr::nest(.IC, data = -c(!!! nest_gr)) %>%
      mutate(gam_out = purrr::map(.data$data, ~gam_fun(.x, args, X.mdl))) %>%
      tidyr::unnest(cols = c(.data$data, .data$gam_out)) %>%
      group_by(!!! nest_gr) %>%
      mutate(
        # Group wise mean
        !! M_X.l0 := mth_switch(.method,!! X.mdl),
        # mean plus random intercept correction
        !! X.l0 := !! M_X.l0 + (!! args[[".X"]] - !! X.mdl) + .data$ran_in.ml,
        !! N.l0 := !! X.l0 * (min(!! args[[".t"]]) - (!! args[[".bl_t"]] / 1e3))
        )
    # Single analysis
  } else {
    IC <- tidyr::nest(.IC, data = -c(!!! gr_by)) %>%
      mutate(
        !! X.mdl := purrr::map(.data$data, ~as.vector(gam_form(.x, args)))
        ) %>%
      tidyr::unnest(cols = c(.data$data, !! X.mdl)) %>%
      group_by(!!! gr_by) %>%
      mutate(
        !! M_X.l0 := mth_switch(.method, !! X.mdl),
        # mean correction
        !! X.l0 := !! M_X.l0 + (!! args[[".X"]]- !! X.mdl),
        !! N.l0 := !! X.l0 * (min(!! args[[".t"]]) - (!! args[[".bl_t"]] / 1e3))
        )
    }

  IC <- ungroup(IC)

  if (.plot) {
    if (is_symbol(get_expr(nest))) {
      facets_gr <- quos(!!! nest_gr)
    } else {
      facets_gr <- quos(!!! gr_by)
    }

    plot_args <- list2(
        .IC = IC,
        .x = args[[".t"]],
        .y = args[[".X"]],
        .flag = NULL,
        .method = "timeseries",
        .plot_type = "static",
        !!! facets_gr,
        .hat = X.mdl,
        .alpha_level = NULL
        )
    ggh <- geom_hline(aes(yintercept = !! M_X.l0), color = "blue", size = 1.1)

    # Plot
    expr(gg_base(!!! plot_args) + ggh) %>%
      eval() %>%
      print()

    # Hide model data
    if (.hide) IC <- fold(IC, c(".mt",".rw", ".ml"))
    IC
    }
  # Hide model data
  if (.hide) IC <- fold(IC, c(".mt",".rw", ".ml"))
  IC
}


mth_switch <- function(mth, arg) {
  switch(
    mth,
    median = call2("median", arg),
    mean  = call2("mean", arg),
    stop("unknown method", call. = FALSE)
    ) %>%
    eval()
}
