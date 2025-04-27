#' PAR mean_t estimates
#' @param spec a par_model_spec object
#' @return numeric vector
#' @export
get_par_mt <- function(spec) {
  assert_valid_par_model_spec(spec)
  get_mt_cpp(spec$Y, spec$X, spec$beta, spec$gamma)
}

is_par_model_spec <- function(x) {
  inherits(x, "par_model_spec")
}

collapse <- function(x) {
  paste(x, collapse = ", ")
}

assert_valid_par_model_spec <- function(spec, .call = parent.frame()) {
  err <- NULL
  if (!is_par_model_spec(spec)) {
    abort("object is not a par_model_spec", call = .call)
  }

  if (!is.numeric(spec$Y)) {
    err <- c(err, "Y must be numeric")
  }

  if (!is.matrix(spec$X)) {
    err <- c(
      err,
      sprintf("X must be a matrix, not %s", collapse(class(spec$X)))
    )
  }

  if (!is.numeric(spec$X)) {
    err <- c(err, sprintf("X must be numeric, not %s", typeof(spec$X)))
  }

  if (!is.numeric(spec$beta)) {
    err <- c(err, sprintf("beta must be numeric, not %s", typeof(spec$beta)))
  }

  if (!is.numeric(spec$gamma)) {
    err <- c(err, sprintf("gamma must be numeric, not %s", typeof(spec$gamma)))
  }

  if (length(spec$Y) != nrow(spec$X)) {
    err <- c(err, "Y must have the same number of rows as X")
  }

  if ((length(spec$beta) + length(spec$gamma)) != ncol(spec$X)) {
    err <- c(
      err,
      "beta and gamma must have the same number of elements as columns in X"
    )
  }

  if (length(err)) {
    rlang::abort(
      c("There are errors with the par_model_spec object", err),
      call = .call,
      class = "par_model_spec_error"
    )
  }
  invisible(spec)
}

get_data_pseudo_complete <- function(spec) {
  y_sym <- rlang::f_lhs(spec$formula)
  time <- spec$model_call$time
  if (any(grepl("^lag[0-9]+$", colnames(spec$.data)))) {
    data <- spec$.data |>
      slice(1L) |>
      dplyr::reframe(
        "{y_sym}" := rev(dplyr::c_across(dplyr::matches("^lag[0-9]+$"))),
        "{time}" := (!!time) - rev(seq_len(length(!!y_sym))),
        ...pseudo_col = 1L
      ) |>
      dplyr::bind_rows(spec$.data) |>
      select(-dplyr::matches("^lag[0-9]+$")) |>
      dplyr::group_by(!!!dplyr::groups(spec$.data)) |>
      dplyr::arrange(!!time, .by_group = TRUE)
    psuedo_col <- pull(data, ...pseudo_col)
    data <- select(data, -...pseudo_col)
    attr(data, "original_data") <- which(is.na(psuedo_col))
    data
  } else {
    data <- spec$.data
    attr(data, "original_data") <- seq_len(nrow(data))
    data
  }
}
