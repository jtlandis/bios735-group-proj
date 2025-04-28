#' Perform EM - Silence a covariate
#' @param spec A model specification object.
#' @param col The column name of the covariate to silence.
#' @param subset A logical expression to subset the data.
#' @param tol The tolerance for convergence.
#' @param max_iter The maximum number of iterations.
#' @param show_plots Whether to print diagnostic information.
#' @examples
#' model_spec <- dplyr::filter(
#'      data_set_tidy,
#'      brand %in% c("B1", "B3"),
#'      item %in% c("1", "3")) |>
#' par_model_mat(
#'  formula = QTY ~ PROMO*brand*item,
#'  time = DATE,
#'  nlag = 4,
#'  groups = c(brand, item)
#' )
#' model_spec
#' em_out <- par_em_effective(model_spec,
#'    col = PROMO,
#'    subset = PROMO == 1,
#'    tol = 1e-6,
#'    show_plots = TRUE)
#' if (require("ggplot2")) {
#'   plot_data <- dplyr::ungroup(em_out$spec$.data) |>
#'      dplyr::mutate(em_mt = get_par_mt(em_out$spec),
#'             mt_original = get_par_mt(em_out$original_spec))
#'   ggplot(plot_data, aes(DATE)) +
#'     geom_line(aes(y = QTY, color = factor("Original Data"))) +
#'     geom_line(aes(y = mt_original, color = factor("Original Model")),
#'               linetype = "dotted", alpha = 0.8) +
#'     geom_rug(data = ~subset(.x, PROMO == 1)) +
#'     geom_line(aes(y = em_mt, color = factor("EM Model")),
#'               linetype = "dotted", alpha = 0.5) +
#'     facet_grid(cols = vars(brand), rows = vars(item)) +
#'     scale_color_manual(values = c("Original Data" = "grey",
#'                                    "Original Model" = "steelblue",
#'                                    "EM Model" = "red")) +
#'     theme_classic() +
#'     labs(color = "Model") +
#'     theme(legend.position = c(.8, .8))
#' }
#' @return a list of em model (spec), orignal model, z embeddings, and loglikelihood
#' @export
par_em_effective <- function(
  spec,
  col,
  subset = NULL,
  tol = 1e-3,
  max_iter = 100,
  show_plots = FALSE
) {
  pi <- 0.5

  data <- get_data_pseudo_complete(spec)
  original_slice <- attr(data, "original_data")
  data <- ungroup(data)
  #data_model <- ungroup(spec$.data)
  res <- fit_par_bfgs(spec)
  beta <- spec$beta <- res$beta
  q <- length(beta)
  gamma <- spec$gamma <- res$gamma
  original_spec <- spec
  new_call <- rlang::call_modify(spec$model_call, data = quote(data))
  new_model <- rlang::new_function(
    args = rlang::exprs(data = ),
    body = rlang::expr({
      model <- !!new_call
      model$beta <- beta
      model$gamma <- gamma
      model
    })
  )

  ll_original <- ll <- res$ll
  N <- nrow(data)
  # probability that promotion is effective
  z <- matrix(0, nrow = N, ncol = 2)
  colnames(z) <- c("EFFECTIVE", "NOT_EFFECTIVE")
  sym_col <- enexpr(col)
  sym_sub <- enexpr(subset)
  which_effective <- if (!is.null(sym_sub)) {
    data |>
      mutate(..i = seq_len(dplyr::n())) |>
      filter(!!sym_sub) |>
      pull(..i)
  } else {
    seq_len(nrow(data))
  }

  # probability that promotion is effective
  z[which_effective, 1] <- 0.9
  z[, 2] <- 1 - z[, 1]

  model_z <- rlang::new_function(
    args = rlang::exprs(data = , z = , i = z),
    body = expr({
      data |>
        mutate(
          ...z = z,
          ...i = i,
          !!sym_col := !!sym_col * ...i
        ) |>
        new_model()
    })
  )

  estimate <- rlang::new_function(
    args = rlang::exprs(data = , z = , i = ),
    body = expr({
      model <- model_z(data, z, i)
      mt <- get_par_mt(model)
      dpois(model$Y, lambda = mt) * model$.data$...z
    })
  )
  density_effective <- function() {
    x <- z[which_effective, 1]
    x <- x[!is.na(x)]
    if (length(x) > 1) {
      plot(
        density(x),
        main = paste(round(100 * pi, 2), "% Promotion Effectiveness")
      )
      rug(x)
    }
  }
  eps <- Inf
  iter <- 0L
  while (eps > tol & iter < max_iter) {
    iter <- iter + 1L
    ll0 <- ll
    #E step
    z[original_slice, 1] <- estimate(data, z[, 1], rep(1, N))
    z[original_slice, 2] <- estimate(data, z[, 2], rep(0, N))
    rsum <- rowSums(z)
    if (any(is_zero <- rsum == 0)) {
      warning("Zero row sums detected")
      rsum[is_zero] <- 1
      z[is_zero, ] <- 0
    }
    z <- z / rsum

    pi <- mean(z[which_effective, 1], na.rm = TRUE)

    #m step
    m_model <- model_z(data, z[, 1])
    res <- bfgs_cpp(
      m_model$Y,
      m_model$X,
      m_model$beta,
      m_model$gamma
    )
    m_model$beta <- beta <- res$beta
    m_model$gamma <- gamma <- res$gamma
    ll <- res$objective
    eps <- abs(ll - ll0) / abs(ll0)
    if (show_plots) {
      density_effective()
    }
  }

  walk_nodes <- function(nodes) {
    N <- length(nodes)
    for (i in rev(seq_len(N))) {
      node <- nodes[[i]]
      if (is.call(node)) {
        nodes[[i]] <- walk_nodes(node)
      } else if (is.name(node)) {
        if (identical(node, sym_col)) {
          nodes[[i]] <- expr((!!sym_col:z))
        }
      }
    }
    nodes
  }

  new_call$formula <- walk_nodes(new_call$formula)
  data <- mutate(data, z = z[, 1])

  model_m <- eval(new_call)
  res <- fit_par_bfgs(model_m)
  model_m$beta <- res$beta
  model_m$gamma <- res$gamma
  list(
    spec = model_m,
    original_spec = original_spec,
    z = z,
    pi = pi,
    ll = ll,
    ll_original = ll_original
  )
}
